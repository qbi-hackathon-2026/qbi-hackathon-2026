# TrimProt

A target-prep assistant for de novo binder design (BindCraft, RFdiffusion, ProteinMPNN,
AlphaFold-Multimer validation, etc.). Given a protein name, gene name, or UniProt
accession, TrimProt:

1. Resolves the target in UniProt (preferring the human-tagged entry).
2. Decides whether extracellular-domain trimming is needed (i.e. does the target have
   a transmembrane region) and, if so, determines the extracellular domain boundaries.
3. Searches RCSB PDB for all structures covering that target and ranks them by
   extracellular-domain sequence coverage, resolution, and whether a bound partner
   (antibody/Fab/nanobody) is present.
4. Maps UniProt sequence positions to the chosen PDB's residue numbering via SIFTS,
   then trims the structure down to just the extracellular domain.
5. Annotates "avoid" residues (glycosylation sites, disulfide-bonded cysteines, other
   PTMs, missing/unresolved residues) that shouldn't be picked as design hotspots.
6. Identifies candidate hotspot residues — either from a real antibody/partner
   interface (heavy-atom contacts) if one exists in the chosen structure, or from an
   inferred-surface-exposure fallback if not.
7. Serves all of this through a web UI: a side-by-side 3D viewer (original vs.
   trimmed), a structured summary explaining every decision, and a downloadable
   trimmed PDB.

## Requirements

- Python 3.10+ (developed/tested on 3.14)
- Outbound internet access to `rest.uniprot.org`, `search.rcsb.org`,
  `data.rcsb.org`, `files.rcsb.org`, and `www.ebi.ac.uk` (PDBe SIFTS) — there is no
  local database, every run hits these APIs live.

## Installation

```bash
cd trimprot/backend/app
python -m pip install fastapi uvicorn requests gemmi
```

## Running it

Start the API + web server from `trimprot/backend/app`:

```bash
cd trimprot/backend/app
python -m uvicorn api:app --host 127.0.0.1 --port 8000
```

Then open **http://127.0.0.1:8000/** in a browser. The page auto-runs the pipeline for
EGFR (`P00533`) on load; use the search box (top of the page) to look up a different
protein/gene name — human-tagged UniProt entries are surfaced first in the autosuggest
dropdown.

`trimprot/frontend/index.html` can also be opened directly as a file (e.g.
double-clicked) instead of through the server's `/` route — it detects a non-http
origin and falls back to talking to `http://127.0.0.1:8000` for its API calls. The
server must still be running for this to work.

### Running pipeline modules standalone

Every backend module has a `if __name__ == "__main__"` block useful for debugging a
single stage without the web server, e.g.:

```bash
cd trimprot/backend/app
python uniprot.py       # fetch + parse EGFR's UniProt entry, print trim decision
python pdb_search.py    # find + rank all PDB structures for EGFR
python sifts.py         # fetch SIFTS mapping for 1YY9, report missing residues
python trim.py          # download 1YY9, trim to ECD, write 1YY9_trimmed.pdb
python annotate.py      # avoid-residues + interface hotspots for 1YY9
python run_egfr.py      # full end-to-end pipeline for EGFR, writes output/ + prints summary JSON
```

## API

| Endpoint | Description |
|---|---|
| `GET /` | Serves the web UI (`frontend/index.html`). |
| `GET /api/search?q=<text>` | Resolves a free-text protein/gene name or UniProt accession to a ranked list of candidates (human-tagged first). |
| `GET /api/run?accession=<acc>&refresh=<bool>` | Runs the full pipeline for a UniProt accession and returns the structured summary JSON (see below). Results are cached in-memory per accession; pass `refresh=true` to force a re-run. |
| `GET /api/files/{filename}` | Serves a generated structure/summary file from `trimprot/backend/output/` (downloads). |

## Output

Each pipeline run writes to `trimprot/backend/output/` (gitignored, regenerated on
every run):

- `{PDB_ID}_full.cif` / `{PDB_ID}_full.pdb` — the original chosen structure.
- `{PDB_ID}_ECD_trimmed.pdb` — the trimmed extracellular-domain-only structure.
- `{PDB_ID}_summary.json` — the structured summary (same JSON returned by `/api/run`).

The summary JSON includes: which UniProt accession/isoform was used, the extracellular
domain boundaries, which PDB was chosen and why (vs. alternatives), how many residues
were trimmed away, the avoid-residue lists (glycosylation/disulfide/other-PTM/missing)
with both UniProt and structure residue numbers, and the candidate hotspot list with
its source (known partner interface vs. inferred surface exposure).

## Known limitations

- Single chain only: the pipeline picks one UniProt-mapped chain in the chosen PDB
  entry; it doesn't currently handle targets that need multiple chains stitched
  together.
- The no-partner hotspot fallback (`annotate.surface_exposed_hotspots`) is a crude
  CB-contact-count heuristic, not true solvent-accessible surface area — `gemmi`
  (the structure library in use) has no built-in SASA calculation, and `freesasa`
  currently has no installable wheel for Python 3.14.
- No predicted-model (AlphaFold) fallback yet when no experimental PDB structure
  exists for a target — `pdb_search.rank_structures` will raise if the RCSB search
  returns zero hits.
