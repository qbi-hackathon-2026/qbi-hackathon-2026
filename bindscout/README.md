# BindScout

Deterministic **target-prep and binding-epitope discovery pipeline** for de novo binder design.

Input: a gene/protein name (or UniProt accession).
Output: a design-ready trimmed target structure, ranked hotspot/epitope-patch
residues, an avoid set, a BindCraft config, and a self-contained 3Dmol.js viewer
that highlights the hotspot residues on the trimmed target.

All structural-biology logic is **pure and deterministic** (no LLM calls). The
exact same code paths are exposed two ways:

- **CLI:** `python -m bindscout.cli --protein EGFR`
- **MCP:** `python -m bindscout.mcp_server` (FastMCP server; one tool per pipeline stage)

## Hard correctness guarantees

1. Everything operates on **author numbering** (`gemmi residue.seqid.num` + insertion code). Never renumbered.
2. UniProt↔PDB(auth) mapping comes from **PDBe SIFTS** (segment-based, offsets + gaps), cross-checked against `_pdbx_sifts_xref_db`.
3. Unobserved UniProt residues → `avoid`.
4. Topology from UniProt **TOPO_DOM** only (never inferred from position). Handles type I and type II correctly.
5. Homo-oligomers download the **biological assembly** by default.
6. **Membrane-proximal exclusion:** the ECD terminus adjacent to TRANSMEM is excluded by a buffer (`--membrane-buffer`, default 12).

## Quickstart

```bash
uv sync --extra dev
uv run bindscout --protein EGFR
uv run pytest -m "not network"     # offline unit tests
uv run pytest                       # full suite (hits UniProt/PDBe/RCSB)
```

Outputs land in `./outputs/<TARGET>/`.

## Web app (interactive viewer)

A local single-page app: type a target (gene name or UniProt accession), run the
pipeline, and view the chosen **RCSB assembly** (left, colored by chain) beside the
**trimmed target** (right) — hotspots in green, the BindCraft patch in orange,
Pfam **domains** (Domain I, II, …) and **chains** labeled, and per-residue **hover**.
The backend runs the pipeline through the MCP server (in-process FastMCP client —
deterministic, no LLM).

```bash
uv sync
uv run python -m bindscout.server    # then open http://127.0.0.1:8000
```

`server.py` exposes `/api/run?target=EGFR` (runs `prepare_target` via MCP),
`/files/<TARGET>/<name>` (serves `trimmed.pdb` / `original.cif`), and serves
`frontend/index.html`.
