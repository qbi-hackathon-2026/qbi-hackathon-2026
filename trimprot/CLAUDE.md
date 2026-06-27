# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

TrimProt: a target-prep assistant for de novo binder design. Given a protein
name/gene/UniProt accession, it resolves the target, decides whether it needs
extracellular-domain trimming, finds and ranks PDB structures, trims to the
extracellular domain via SIFTS residue mapping, annotates avoid-residues and
hotspots, and serves all of it through a small web UI (FastAPI backend +
single-page vanilla-JS/NGL.js frontend). See `README.md` for the full pipeline
description, API table, and known limitations — read it before making
behavioral changes.

## Running it

```bash
cd trimprot/backend/app
python -m pip install fastapi uvicorn requests gemmi
python -m uvicorn api:app --host 127.0.0.1 --port 8000
```

Then open http://127.0.0.1:8000/. There is no test suite and no linter
configured — verification is done by running the server and checking the
`/api/run` JSON output and/or the rendered page (see "Verifying changes"
below).

Every backend module (`uniprot.py`, `pdb_search.py`, `sifts.py`, `trim.py`,
`annotate.py`, `run_egfr.py`) has a `__main__` block that exercises just that
module against EGFR (`P00533`) / PDB `1YY9` — run any of them directly to
debug one pipeline stage in isolation without the web server.

## Architecture

The pipeline is a strict linear chain; each stage's output is the next
stage's input. `run_egfr.py:run(accession)` is the orchestrator — read it
first when tracing how data flows end to end.

```
uniprot.py        -> fetch UniProt entry, extract features, decide if trim is
                      needed, find extracellular domain boundaries (UniProt
                      numbering)
pdb_search.py      -> find all PDB entries referencing that UniProt accession
                      (RCSB search API), fetch resolution/entity/sequence
                      details in batched GraphQL, rank by ECD coverage +
                      resolution + partner-chain bonus
sifts.py           -> map UniProt residue positions <-> PDB label_seq_id via
                      PDBe SIFTS, then label_seq_id <-> author residue number
                      via the mmCIF `_pdbx_poly_seq_scheme` table
trim.py            -> download the chosen PDB's full mmCIF, use sifts.py to
                      figure out which author-numbered residues fall inside
                      the extracellular domain, build a new gemmi.Structure
                      containing only those
annotate.py        -> map avoid-residue/PTM positions (from uniprot.py
                      features) into the structure's author numbering;
                      detect hotspots either via real partner-chain heavy-atom
                      contacts (gemmi.NeighborSearch) or, if no partner chain
                      exists, a surface-exposure fallback heuristic
summary.py         -> assemble the final structured JSON (target info, domain
                      boundaries, structure choice + reasoning, trim stats,
                      avoid-residue lists, hotspot list + source)
api.py             -> FastAPI app: /api/search, /api/run, /api/files/*, and
                      serves frontend/index.html at /
frontend/index.html -> single-file UI: NGL.js dual 3D viewer, search box with
                      autosuggest, summary panel, download links. All wiring
                      (color schemes, residue highlighting) lives in one
                      <script> block.
```

### Residue numbering — the part most likely to bite you

There are at least three different residue-numbering schemes in play, and
mixing them up silently produces wrong-looking-but-not-crashing results:

- **UniProt position** — 1-indexed position in the canonical UniProt
  sequence. This is what `uniprot.py`'s feature lists (glycosylation sites,
  disulfide bonds, PTMs, topological domains) are expressed in.
- **PDB `label_seq_id`** — sequential, gap-free position within one PDB
  polymer entity's full construct (including unobserved/missing residues).
  SIFTS (`sifts.fetch_unp_segments`) maps UniProt positions to this via a
  per-segment linear offset (`label_seq_start - unp_start`).
- **PDB author residue number (`auth_seq_num` / `seqid.num`)** — the number
  actually printed in the PDB/mmCIF file, recovered per-`label_seq_id` from
  `_pdbx_poly_seq_scheme` (`sifts.parse_poly_seq_scheme`). This is the only
  scheme gemmi's `Structure`/`Chain`/`Residue` objects expose directly
  (`res.seqid.num`), and the only one valid for selecting atoms in a
  downloaded structure or for NGL viewer selections in the frontend.

`annotate.unp_to_auth_map` builds the UniProt -> auth_seq_num dict that
everything downstream (avoid-residue lookup, hotspot exclusion, the
viewer's auth-residue lists in `run_egfr.py`'s `result_summary["viewer"]`)
depends on. If you add a new UniProt-derived annotation, route it through
this map rather than assuming any numbering scheme lines up 1:1.

### Partner-chain vs. surface-exposure hotspot detection

`annotate.find_partner_chains` treats any non-target polymer chain in the
downloaded structure as a "partner" (antibody, ligand protein, etc.) with
no further filtering. If a partner chain exists, hotspots come from real
heavy-atom contacts (`annotate.interface_hotspots`); otherwise from
`annotate.surface_exposed_hotspots`, which is a CB-contact-count heuristic,
**not** true SASA — gemmi has no built-in solvent-accessible-surface-area
calculation, and `freesasa` does not currently have an installable wheel for
Python 3.14 in this environment. Don't treat the surface-exposure numbers as
biophysically rigorous; they're a rough fallback.

### The `_cache` in api.py

`/api/run` caches one result per accession in an in-memory dict (`api.py`'s
`_cache`). It does **not** check whether the underlying output files in
`trimprot/backend/output/` still exist — if those files are deleted (they're
gitignored, regenerated per run) while the server process keeps running, the
cache will keep returning stale file references and `/api/files/...`
will 404. Pass `?refresh=true` to force re-running the pipeline, or restart
the server.

### Frontend

`frontend/index.html` is a single static file with no build step. It uses
relative `/api/...` paths when served from the FastAPI app's own origin, but
falls back to an absolute `http://127.0.0.1:8000` base (`API_BASE` constant
near the top of the `<script>` block) when opened from a non-http origin
(e.g. a local file or a webview), since relative fetches can't resolve
without an http(s) origin.

NGL viewer gotcha: call `stage.handleResize()` before `stage.autoView()`
when loading a structure — `autoView()` frames the camera based on the
canvas's current size, and if that's called before the container has its
real laid-out size, the structure renders far too small/off-screen to see.

## Verifying changes

There's no automated test suite. To check a change actually works:

1. Start the server (see "Running it").
2. Hit `curl http://127.0.0.1:8000/api/run?accession=<acc>&refresh=true` and
   inspect the JSON (avoid-residue counts, hotspot list, chosen PDB) for
   sanity.
3. Load the page in a browser and confirm the 3D viewer panes actually
   render (not just that the canvas exists — NGL can silently fail to frame
   content; see the `handleResize`/`autoView` ordering note above) and that
   the summary/avoid/hotspot panels populate.

Good secondary test targets beyond EGFR (`P00533` / PDB `1YY9`, has a bound
antibody Fab): TNF-alpha (`P01375` / PDB `8ZUI`) exercises a different
transmembrane topology and partner-chain hotspot path end to end.
