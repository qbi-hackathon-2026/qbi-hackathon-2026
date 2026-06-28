"""Local web app for BindScout.

A thin FastAPI server that serves a single-page UI and runs the pipeline by
calling the MCP server's `prepare_target` tool through an in-process FastMCP
client (deterministic — no LLM, no agent). The browser cannot run gemmi /
network calls, so this server is the only backend: input target -> run pipeline
-> return the chosen RCSB assembly, the trimmed target, and the summary.

    uv run python -m bindscout.server      # then open http://127.0.0.1:8000
"""
from __future__ import annotations

import os
import re
from pathlib import Path
from xml.etree import ElementTree as ET

import requests as _requests
from dotenv import load_dotenv
from fastapi import FastAPI
from fastapi.responses import FileResponse, HTMLResponse, JSONResponse, Response
from fastmcp import Client
from pydantic import BaseModel

PROTTER_URL = "https://protter.ethz.ch/create"

# Matches the web app's default styles (signal peptide, disulfide bonds, variants, PTMs).
# Style descriptors are the URL parameter keys; UniProt feature codes are the values.
# Special chars are left unencoded — Protter's server parses them raw.
_SVG_PAD_TOP    = 5
_SVG_PAD_BOTTOM = 5
_SVG_PAD_LEFT   = 5
_SVG_PAD_RIGHT  = 45  # legend text anchors at rightmost circle; text overflows right


def _crop_svg(svg_bytes: bytes) -> bytes:
    """Rewrite the SVG viewBox to the tight content bounding box + padding.

    Skips <defs> children (pattern/glyph definitions) since they are not
    directly rendered and their coordinates are relative to their own viewport.
    """
    try:
        ET.register_namespace("", "http://www.w3.org/2000/svg")
        ET.register_namespace("xlink", "http://www.w3.org/1999/xlink")
        root = ET.fromstring(svg_bytes)
        SVG_NS = "http://www.w3.org/2000/svg"
        defs_tag = f"{{{SVG_NS}}}defs"

        mins_x, mins_y, maxs_x, maxs_y = [], [], [], []

        def _walk(el: ET.Element) -> None:
            for child in el:
                if child.tag == defs_tag:
                    continue  # skip glyph/pattern definitions
                tag = child.tag.split("}")[-1] if "}" in child.tag else child.tag
                try:
                    if tag == "circle":
                        cx, cy, r = float(child.get("cx", 0)), float(child.get("cy", 0)), float(child.get("r", 0))
                        mins_x.append(cx - r); maxs_x.append(cx + r)
                        mins_y.append(cy - r); maxs_y.append(cy + r)
                    elif tag == "rect":
                        x, y = float(child.get("x", 0)), float(child.get("y", 0))
                        w, h = float(child.get("width", 0)), float(child.get("height", 0))
                        mins_x.append(x); maxs_x.append(x + w)
                        mins_y.append(y); maxs_y.append(y + h)
                    elif tag in ("text", "use"):
                        if child.get("x") is None or child.get("y") is None:
                            continue  # skip hidden internal data elements
                        x, y = float(child.get("x")), float(child.get("y"))
                        mins_x.append(x); maxs_x.append(x)
                        mins_y.append(y); maxs_y.append(y)
                except (TypeError, ValueError):
                    pass
                _walk(child)

        _walk(root)

        if not mins_x:
            return svg_bytes

        vx = min(mins_x) - _SVG_PAD_LEFT
        vy = min(mins_y) - _SVG_PAD_TOP
        vw = max(maxs_x) - vx + _SVG_PAD_RIGHT
        vh = max(maxs_y) - vy + _SVG_PAD_BOTTOM

        root.set("viewBox", f"{vx:.2f} {vy:.2f} {vw:.2f} {vh:.2f}")
        # Remove fixed pt dimensions so the browser scales from the viewBox aspect ratio.
        root.attrib.pop("width", None)
        root.attrib.pop("height", None)
        return ET.tostring(root, encoding="unicode", xml_declaration=False).encode()
    except Exception:
        return svg_bytes  # fall back to unmodified SVG on any parse error


_PROTTER_DEFAULT_STYLES = (
    "&n:signal peptide,fc:red,bc:red=UP.SIGNAL"
    "&n:disulfide bonds,s:box,fc:greenyellow,bc:greenyellow=UP.DISULFID"
    "&n:variants,s:diamond,fc:orange,bc:orange=UP.VARIANT"
    "&n:PTMs,s:box,fc:forestgreen,bc:forestgreen=UP.CARBOHYD,UP.MOD_RES"
)

ROOT = Path(__file__).resolve().parents[2]
OUTPUTS = ROOT / "outputs"
INDEX = ROOT / "frontend" / "index.html"

# Load each user's own ANTHROPIC_API_KEY (and any other secrets) from a local,
# gitignored `bindscout/.env` — so nobody has to export it per-terminal and no
# key is ever committed. Anchored to the package dir so it works regardless of
# the working directory the server is launched from. Must run BEFORE chat.py
# constructs the Anthropic client (which reads ANTHROPIC_API_KEY from the env).
# An already-set environment variable wins (override=False), so an explicit
# `export` still takes precedence over the file.
load_dotenv(ROOT / ".env")

from .chat import chat as run_chat        # noqa: E402  (import after load_dotenv)
from .mcp_server import mcp                # noqa: E402
from .uniprot import search_proteins       # noqa: E402

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
        url = (f"{PROTTER_URL}?up={accession.upper()}&tm=auto"
               f"{_PROTTER_DEFAULT_STYLES}&legend&format=svg")
        resp = _requests.get(url, timeout=15)
        resp.raise_for_status()
        return Response(content=_crop_svg(resp.content), media_type="image/svg+xml")
    except _requests.HTTPError as exc:
        return JSONResponse({"error": f"Protter returned {exc.response.status_code}"},
                            status_code=502)
    except Exception as exc:
        return JSONResponse({"error": str(exc)}, status_code=502)


class ChatTurn(BaseModel):
    role: str       # "user" | "assistant"
    content: str


class ChatRequest(BaseModel):
    messages: list[ChatTurn]
    context: dict | None = None   # summary of the target loaded in the UI, if any


@app.post("/api/chat")
async def chat(req: ChatRequest):
    """Conversational assistant: Claude orchestrates the deterministic pipeline tools."""
    if not os.getenv("ANTHROPIC_API_KEY"):
        return JSONResponse(
            {"error": "The assistant is unavailable — set ANTHROPIC_API_KEY in the "
                      "server environment and restart."}, status_code=503)

    msgs = [{"role": m.role, "content": m.content} for m in req.messages]
    if not msgs:
        return JSONResponse({"error": "empty conversation"}, status_code=400)

    try:
        result = await run_chat(msgs, outdir=str(OUTPUTS), context=req.context)
    except Exception as exc:
        msg = re.sub(r"^Error calling tool '[^']*':\s*", "", str(exc))
        return JSONResponse({"error": msg}, status_code=400)

    # If the assistant re-ran the full pipeline, hand back file URLs + summary so
    # the UI can refresh the 3D viewers in place.
    files = None
    summary = None
    prepared = result.get("prepared")
    if prepared:
        summary = prepared.get("summary") or {}
        name = summary.get("target")
        if name:
            files = {
                "original": f"/files/{name}/original.cif",
                "trimmed": f"/files/{name}/trimmed.pdb",
            }
    return {
        "reply": result["reply"],
        "tools_used": result["tools_used"],
        "files": files,
        "summary": summary,
    }


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
