"""Conversational assistant over the BindScout pipeline.

A thin agent loop: Claude orchestrates the SAME deterministic FastMCP tools the
rest of the app uses (resolve_target, structures, prepare_target, ...). The model
decides which tool to call and with what parameters, and explains the result; it
never computes structural biology itself. Every fact the assistant reports comes
from a deterministic tool call.

The loop is intentionally manual (not the SDK tool runner) so it can bridge the
in-process FastMCP `Client(mcp)` directly — the same client `server.py` already
uses for `/api/run` — without a second MCP transport.
"""
from __future__ import annotations

import json
from typing import Any, Optional

from fastmcp import Client

from .mcp_server import mcp

MODEL = "claude-opus-4-8"
MAX_TOOL_ROUNDS = 8          # cap the agent loop so a confused turn can't run away
MAX_TOOL_RESULT_CHARS = 20000  # keep large summaries from blowing up the context

SYSTEM = """\
You are BindScout's assistant. BindScout is a deterministic pipeline that prepares
a design-ready protein target for de novo binder design: it resolves a UniProt
entry, ranks experimental/AlphaFold structures, trims to the relevant
extracellular domain, picks hotspot residues and an epitope patch, builds an
"avoid" set, and emits a BindCraft config and a 3D view.

You have tools that run the real pipeline stages. Use them — never invent
structural facts. When a user asks something you can answer by running a stage
(which structure was picked and why, the ECD ranges, glycosylation sites, the
avoid set, re-trimming with different parameters), call the tool and report what
it returned.

Key tools:
- resolve_target / extracellular_ranges / membrane_proximal / glycosylation:
  facts about the protein from UniProt.
- structures: the ranked candidate structures and the selection reasoning.
- prepare_target: run the FULL pipeline and (re)emit artifacts. It accepts
  parameters you can vary on request: prefer_antibody, membrane_buffer,
  patch_radius, patch_size, no_patch. Use it to re-run with different settings
  (e.g. a larger membrane buffer, or a no-patch run).

Honesty rules — this matters:
- BindScout does NOT compute fold stability, folding free energy (ΔΔG),
  expressibility, or binding affinity. If asked whether an isolated domain "will
  fold" or "stay stable", say plainly that BindScout doesn't model that, give
  only what the structure data supports (e.g. whether the region is a discrete
  Pfam/topology domain), and suggest a dedicated tool if they need a real answer.
- Distinguish what a tool returned (fact) from your own reasoning (clearly
  flagged). Never present a guess as a computed result.
- Keep answers short and concrete. Lead with the answer. Cite specific residues,
  PDB IDs, and numbers from tool output.
"""


def _anthropic_tools(mcp_tools: list) -> list[dict]:
    """Convert FastMCP tool definitions into Anthropic tool schemas."""
    out: list[dict] = []
    for t in mcp_tools:
        out.append({
            "name": t.name,
            "description": (t.description or "").strip(),
            "input_schema": t.inputSchema or {"type": "object", "properties": {}},
        })
    return out


async def chat(messages: list[dict], *, outdir: str,
               context: Optional[dict] = None) -> dict[str, Any]:
    """Run one assistant turn over the conversation, executing tools as needed.

    `context` is the summary of the target currently loaded in the UI (if any),
    so the assistant knows what the user means by "this protein"/"this PDB"
    without having to re-resolve it.

    Returns {reply, tools_used, prepared} where `prepared` is the prepare_target
    payload if the assistant (re)ran the full pipeline this turn (so the UI can
    refresh the viewers), else None.
    """
    try:
        import anthropic
    except ImportError:
        raise RuntimeError(
            "The assistant needs the 'anthropic' package — run `uv sync`.")

    client = anthropic.AsyncAnthropic()  # reads ANTHROPIC_API_KEY from the env

    system = SYSTEM
    if context:
        ctx = json.dumps(context, default=str)[:8000]
        system += (
            "\n\n--- CURRENTLY LOADED TARGET ---\n"
            "The user has this target open in the UI right now. When they say "
            '"this protein", "this PDB", "the chosen structure", "the hotspots", '
            "etc., they mean THIS one — do not ask them which protein. Answer from "
            "this summary directly (it already includes the structure-selection "
            "reasoning); only call a tool when you need data it doesn't contain or "
            "to re-run the pipeline with new parameters:\n" + ctx)

    convo: list[dict] = list(messages)
    tools_used: list[dict] = []
    prepared: Optional[dict] = None

    async with Client(mcp) as mcp_client:
        tools = _anthropic_tools(await mcp_client.list_tools())

        resp = None
        for _ in range(MAX_TOOL_ROUNDS):
            resp = await client.messages.create(
                model=MODEL, max_tokens=4096,
                system=system, tools=tools, messages=convo,
            )
            if resp.stop_reason != "tool_use":
                break

            # Echo the assistant turn (incl. tool_use blocks) back into history.
            convo.append({"role": "assistant", "content": resp.content})

            results = []
            for block in resp.content:
                if block.type != "tool_use":
                    continue
                args = dict(block.input or {})
                # Pin emitted artifacts to the app's outputs dir so the viewers
                # and /files route can find them.
                if block.name == "prepare_target":
                    args.setdefault("outdir", outdir)
                try:
                    r = await mcp_client.call_tool(block.name, args)
                    payload = r.data
                    if block.name == "prepare_target":
                        prepared = payload
                    content = json.dumps(payload, default=str)[:MAX_TOOL_RESULT_CHARS]
                    is_error = False
                except Exception as exc:  # surface the failure to the model
                    content = f"error: {exc}"
                    is_error = True
                tools_used.append({"tool": block.name, "input": args})
                results.append({
                    "type": "tool_result",
                    "tool_use_id": block.id,
                    "content": content,
                    "is_error": is_error,
                })
            convo.append({"role": "user", "content": results})

        reply = ""
        if resp is not None:
            reply = "".join(b.text for b in resp.content if b.type == "text").strip()

    return {"reply": reply, "tools_used": tools_used, "prepared": prepared}
