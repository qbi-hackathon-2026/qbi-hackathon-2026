# TrimProt

**Deterministic protein target-preparation for de novo binder design.**

Give TrimProt a protein — a gene name, protein name, or UniProt accession — and it
prepares a design-ready target: it resolves the UniProt entry, finds and ranks
structures, trims to the relevant domain, identifies candidate **hotspot**
residues and an epitope **patch**, builds an **avoid** set, and renders it all in
an interactive 3D viewer. When no suitable experimental structure exists, it falls
back to the **AlphaFold** model automatically.

All structural-biology logic is **pure and deterministic** — no LLM calls in the
pipeline. The same code paths are available three ways: a web app, a CLI, and an
**MCP server** for AI tooling.

```bash
cd trimprot
uv sync
uv run python -m trimprot.server      # then open http://127.0.0.1:8000
```

## Three ways to use it

| Interface | Command | For |
|---|---|---|
| **Web app** | `uv run python -m trimprot.server` | Interactive: search a target, view the chosen assembly beside the trimmed target, inspect hotspots/patch/warnings. |
| **CLI** | `uv run trimprot --protein EGFR` | Scripting / batch; outputs land in `trimprot/outputs/<TARGET>/`. |
| **MCP** | `uv run python -m trimprot.mcp_server` | Agents/AI tooling — one tool per pipeline stage (see [MCP server](#mcp-server-for-ai-tooling)). |

## What it produces

For a target it emits, under `trimprot/outputs/<TARGET>/`:

- `trimmed.pdb` — the design-ready trimmed target (author numbering preserved)
- `original.cif` — the chosen assembly (or AlphaFold model)
- `summary.json` — structured result: chosen structure + reasoning, topology, ECD
  ranges, hotspots, epitope patch, avoid set, warnings
- `hotspots.csv` — ranked hotspots with an `in_patch` flag
- a BindCraft config for the epitope patch

The web app renders two panels: the **chosen RCSB assembly** (colored by chain)
beside the **trimmed target** with hotspots in green and the patch in orange, plus
chain/domain labels, per-residue hover, and "focus patch."

## Correctness guarantees

1. Everything operates on **author numbering** (`gemmi residue.seqid.num` + insertion code). Never renumbered.
2. UniProt↔PDB(auth) mapping comes from **PDBe SIFTS** (segment-based, offsets + gaps), cross-checked against `_pdbx_sifts_xref_db`.
3. Unobserved UniProt residues → `avoid`.
4. Topology from UniProt **TOPO_DOM** only (never inferred from position). Handles type I and type II correctly.
5. Homo-oligomers download the **biological assembly** by default.
6. **Membrane-proximal exclusion:** the ECD terminus adjacent to TRANSMEM is excluded by a buffer (`--membrane-buffer`, default 12).
7. **AlphaFold fallback:** when no experimental structure exists (or none covers
   ≥ 50% of the ECD), the AlphaFold model is used, with pLDDT masking — residues
   ≤ 50 treated as unobserved and 50–70 as low-confidence, both folded into the
   avoid set.

## Web API

The server (`trimprot/src/trimprot/server.py`) exposes:

- `GET /api/search?q=<text>` — typeahead: gene/protein name or accession →
  ranked UniProt candidates.
- `GET /api/run?target=<gene|accession>` — runs the pipeline (via the MCP
  `prepare_target` tool, in-process) and returns the summary + file URLs.
- `GET /files/<TARGET>/<name>` — serves emitted artifacts (`trimmed.pdb`,
  `original.cif`, `summary.json`, `hotspots.csv`).
- `GET /` — the single-page viewer.

## MCP server (for AI tooling)

`trimprot.mcp_server` is a [FastMCP](https://github.com/jlowin/fastmcp) server that
exposes each deterministic pipeline stage as a tool, so an agent can run the whole
pipeline or orchestrate individual stages:

| Tool | Purpose |
|---|---|
| `resolve_target` | gene/protein/accession → UniProt record + features |
| `extracellular_ranges` | ECD + transmembrane ranges (UniProt numbering) |
| `membrane_proximal` | membrane-proximal terminus + excluded buffer |
| `glycosylation` | predicted glycosylation sites (CARBOHYD ∪ N-X-[S/T]) |
| `structures` | ranked candidate structures + chosen target/partner chains |
| `sifts_mapping` | per-residue UniProt↔PDB(auth) mapping with observed flags |
| `prepare_target` | run the full pipeline and emit the artifact set |

```bash
cd trimprot
uv run python -m trimprot.mcp_server
```

The web app itself calls `prepare_target` through an in-process FastMCP client —
the exact same deterministic path, no LLM.

## Development

```bash
cd trimprot
uv sync --extra dev
uv run pytest -m "not network"     # offline unit tests (run in CI)
uv run pytest                       # full suite (hits UniProt/PDBe/RCSB)
```

See [`trimprot/README.md`](trimprot/README.md) for the pipeline architecture and
[`CONTRIBUTING.md`](CONTRIBUTING.md) for how to contribute.

`validate_structures.py` is a read-only auditor: after generating outputs it
re-queries RCSB to verify each chosen structure's antibody/partner labelling, and
exits nonzero if any target is mislabeled.

## License

[MIT](LICENSE) © TrimProt contributors.
