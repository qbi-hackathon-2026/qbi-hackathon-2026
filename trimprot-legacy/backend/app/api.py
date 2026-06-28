"""FastAPI server: runs the TrimProt pipeline and serves the viewer webpage."""
import os
from enum import Enum

from fastapi import FastAPI, HTTPException
from fastapi.responses import FileResponse, HTMLResponse
from fastapi.staticfiles import StaticFiles

import run_egfr
import uniprot

APP_DIR = os.path.dirname(os.path.abspath(__file__))
OUT_DIR = os.path.normpath(os.path.join(APP_DIR, "..", "output"))
STATIC_DIR = os.path.normpath(os.path.join(APP_DIR, "..", "..", "frontend"))
os.makedirs(OUT_DIR, exist_ok=True)

app = FastAPI(title="TrimProt")

_cache: dict[str, dict] = {}

# MIME types for structure file formats
_MEDIA_TYPES = {
    ".pdb": "chemical/x-pdb",
    ".cif": "chemical/x-mmcif",
}


class StructureType(str, Enum):
    trimmed = "trimmed"   # ECD-only trimmed structure
    full    = "full"      # complete original structure (PDB or AlphaFold)
    cif     = "cif"       # full structure in mmCIF format


def _require_cached(accession: str) -> dict:
    """Return cached result or raise 404 with a helpful message."""
    if accession not in _cache:
        raise HTTPException(
            status_code=404,
            detail=f"Pipeline not yet run for {accession!r}. Call /api/run?accession={accession} first.",
        )
    return _cache[accession]


def _require_file(filename: str) -> str:
    """Return absolute path to an output file, or raise 404 if missing from disk."""
    path = os.path.join(OUT_DIR, filename)
    if not os.path.isfile(path):
        raise HTTPException(
            status_code=404,
            detail=f"File {filename!r} not found on disk. Re-run with ?refresh=true to regenerate.",
        )
    return path


@app.get("/api/search")
def search(q: str):
    try:
        return uniprot.search_proteins(q)
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


@app.get("/api/run")
def run_pipeline(accession: str = "P00533", refresh: bool = False):
    if accession not in _cache or refresh:
        try:
            _cache[accession] = run_egfr.run(accession)
        except Exception as e:
            raise HTTPException(status_code=500, detail=str(e))
    return _cache[accession]


# ── Download endpoints ────────────────────────────────────────────────────────
# These are the endpoints the UI should use — they're accession-aware and
# return consistent responses regardless of whether the result came from a
# crystal structure or an AlphaFold model.

@app.get("/api/downloads/{accession}")
def list_downloads(accession: str):
    """
    Return metadata describing what structure files are available for download.

    The UI can call this after /api/run to get the exact URLs and human-readable
    labels for each file, without needing to know any filename conventions.

    Response shape:
      {
        "accession": "B6A8C7",
        "alphafold_fallback": true,
        "source": "AlphaFold",
        "files": {
          "trimmed": { "url": "/api/download/B6A8C7/trimmed", "filename": "...", "description": "...", "format": "PDB" },
          "full":    { ... },
          "cif":     { ... }
        }
      }
    """
    result = _require_cached(accession)
    viewer = result["viewer"]
    is_af = result.get("alphafold_fallback", False)
    source = "AlphaFold" if is_af else f"PDB {viewer['pdb_id']}"

    base = f"/api/download/{accession}"
    return {
        "accession": accession,
        "alphafold_fallback": is_af,
        "source": source,
        "files": {
            "trimmed": {
                "url": f"{base}/trimmed",
                "filename": viewer["trimmed_file"],
                "description": f"Trimmed extracellular domain ({source})",
                "format": "PDB",
            },
            "full": {
                "url": f"{base}/full",
                "filename": viewer["original_file"],
                "description": f"Full structure ({source})",
                "format": "PDB",
            },
            "cif": {
                "url": f"{base}/cif",
                "filename": viewer["original_cif_file"],
                "description": f"Full structure ({source}, mmCIF)",
                "format": "mmCIF",
            },
        },
    }


@app.get("/api/download/{accession}/{file_type}")
def download_structure(accession: str, file_type: StructureType):
    """
    Download a structure file by accession and type.

    file_type:
      trimmed — ECD-only trimmed PDB (the main design input)
      full    — complete original structure as PDB
      cif     — complete original structure in mmCIF format

    Works identically whether the result came from a crystal structure or
    AlphaFold — the viewer block in the pipeline result always carries the
    correct filenames for both paths.

    Returns 404 if /api/run hasn't been called yet for this accession, or if
    the output files have been deleted since the last run (re-run with
    ?refresh=true to regenerate).
    """
    result = _require_cached(accession)
    viewer = result["viewer"]

    filename = {
        StructureType.trimmed: viewer["trimmed_file"],
        StructureType.full:    viewer["original_file"],
        StructureType.cif:     viewer["original_cif_file"],
    }[file_type]

    path = _require_file(filename)
    ext = os.path.splitext(filename)[1].lower()
    media_type = _MEDIA_TYPES.get(ext, "application/octet-stream")
    return FileResponse(path, filename=filename, media_type=media_type)


# ── Legacy generic file endpoint (kept for NGL viewer fetch calls) ────────────

@app.get("/api/files/{filename}")
def get_file(filename: str):
    path = os.path.join(OUT_DIR, filename)
    if not os.path.isfile(path):
        raise HTTPException(status_code=404, detail="file not found")
    ext = os.path.splitext(filename)[1].lower()
    return FileResponse(path, filename=filename, media_type=_MEDIA_TYPES.get(ext, "application/octet-stream"))


@app.get("/", response_class=HTMLResponse)
def index():
    with open(os.path.join(STATIC_DIR, "index.html"), encoding="utf-8") as f:
        return f.read()


app.mount("/static", StaticFiles(directory=STATIC_DIR), name="static")


if __name__ == "__main__":
    import uvicorn
    uvicorn.run(app, host="127.0.0.1", port=8000)
