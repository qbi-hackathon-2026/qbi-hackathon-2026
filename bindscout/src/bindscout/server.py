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

from dotenv import load_dotenv
from fastapi import FastAPI
from fastapi.responses import FileResponse, HTMLResponse, JSONResponse
from fastmcp import Client
from pydantic import BaseModel

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
