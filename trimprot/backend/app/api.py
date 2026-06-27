"""FastAPI server: runs the EGFR target-prep pipeline and serves the viewer webpage."""
import os
import json

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


@app.get("/api/search")
def search(q: str):
    try:
        return uniprot.search_proteins(q)
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


@app.get("/api/run")
def run_pipeline(accession: str = "P00533", refresh: bool = False, verbose: bool = False):
    cache_key = (accession, verbose)
    if cache_key not in _cache or refresh:
        try:
            _cache[cache_key] = run_egfr.run(accession, verbose=verbose)
        except Exception as e:
            raise HTTPException(status_code=500, detail=str(e))
    return _cache[cache_key]


@app.get("/api/files/{filename}")
def get_file(filename: str):
    path = os.path.join(OUT_DIR, filename)
    if not os.path.isfile(path):
        raise HTTPException(status_code=404, detail="file not found")
    return FileResponse(path, filename=filename)


@app.get("/", response_class=HTMLResponse)
def index():
    with open(os.path.join(STATIC_DIR, "index.html"), encoding="utf-8") as f:
        return f.read()


app.mount("/static", StaticFiles(directory=STATIC_DIR), name="static")


if __name__ == "__main__":
    import uvicorn
    uvicorn.run(app, host="127.0.0.1", port=8000)
