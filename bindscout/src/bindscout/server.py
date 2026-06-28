"""Local web app for BindScout.

A thin FastAPI server that serves a single-page UI and runs the pipeline by
calling the MCP server's `prepare_target` tool through an in-process FastMCP
client (deterministic — no LLM, no agent). The browser cannot run gemmi /
network calls, so this server is the only backend: input target -> run pipeline
-> return the chosen RCSB assembly, the trimmed target, and the summary.

    uv run python -m bindscout.server      # then open http://127.0.0.1:8000
"""
from __future__ import annotations

import re
from pathlib import Path

import requests as _requests
from fastapi import FastAPI
from fastapi.responses import FileResponse, HTMLResponse, JSONResponse, Response
from fastmcp import Client

from .mcp_server import mcp
from .uniprot import search_proteins

PROTTER_URL = "https://protter.ethz.ch/create"

ROOT = Path(__file__).resolve().parents[2]
OUTPUTS = ROOT / "outputs"
INDEX = ROOT / "frontend" / "index.html"

# UniProt accession pattern (official regex). Anything else is treated as a name.
_ACCESSION = re.compile(
    r"^([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})$")

app = FastAPI(title="BindScout")


@app.get("/", response_class=HTMLResponse)
def index() -> HTMLResponse:
    # Explicit encoding: Path.read_text() defaults to the platform locale codepage
    # (cp1252 on Windows), which mis-decodes the file's UTF-8 punctuation (·, Å, →)
    # into mojibake once HTMLResponse re-encodes it as UTF-8.
    return HTMLResponse(INDEX.read_text(encoding="utf-8"))


@app.get("/api/search")
def search(q: str):
    """Typeahead: free-text gene/protein name (or accession) -> UniProt candidates."""
    try:
        return search_proteins(q)
    except Exception as exc:
        return JSONResponse({"error": str(exc)}, status_code=502)


@app.get("/api/run")
async def run(target: str):
    """Resolve the target, run the pipeline via the MCP tool, return summary + file URLs."""
    target = (target or "").strip()
    if not target:
        return JSONResponse({"error": "empty target"}, status_code=400)

    args: dict = {"prefer_antibody": True, "outdir": str(OUTPUTS)}
    if _ACCESSION.match(target.upper()):
        args["accession"] = target.upper()
    else:
        args["name"] = target

    try:
        async with Client(mcp) as client:
            result = await client.call_tool("prepare_target", args)
    except Exception as exc:  # surface a clean error to the UI
        msg = str(exc)
        # Strip the MCP wrapper ("Error calling tool 'prepare_target': ...").
        msg = re.sub(r"^Error calling tool '[^']*':\s*", "", msg)
        if "no PDB structures found" in msg:
            msg = ("No structures are available for this target — PDBe/RCSB have "
                   "no experimental or predicted models indexed for it, so there is "
                   "nothing to trim. Try a different target.")
        return JSONResponse({"error": msg}, status_code=400)

    data = result.data or {}
    summary = data.get("summary", {})
    name = summary.get("target") or target
    return {
        "summary": summary,
        "files": {
            "original": f"/files/{name}/original.cif",
            "trimmed": f"/files/{name}/trimmed.pdb",
        },
    }


@app.get("/api/protter/{accession}")
def protter(accession: str):
    """Proxy the Protter topology SVG for a UniProt accession."""
    if not _ACCESSION.match(accession.upper()):
        return JSONResponse({"error": "invalid accession"}, status_code=400)
    try:
        resp = _requests.get(
            PROTTER_URL,
            params={"up": accession.upper(), "format": "svg"},
            timeout=15,
        )
        resp.raise_for_status()
        return Response(content=resp.content, media_type="image/svg+xml")
    except _requests.HTTPError as exc:
        return JSONResponse({"error": f"Protter returned {exc.response.status_code}"},
                            status_code=502)
    except Exception as exc:
        return JSONResponse({"error": str(exc)}, status_code=502)


@app.get("/files/{target}/{name}")
def files(target: str, name: str):
    """Serve an emitted artifact (trimmed.pdb / original.cif) for the viewers."""
    # prevent path traversal
    if "/" in name or ".." in name or "/" in target or ".." in target:
        return JSONResponse({"error": "bad path"}, status_code=400)
    path = OUTPUTS / target / name
    if not path.is_file():
        return JSONResponse({"error": "not found"}, status_code=404)
    return FileResponse(path)


def main() -> None:
    import uvicorn
    uvicorn.run(app, host="127.0.0.1", port=8000)


if __name__ == "__main__":
    main()
