"""HTTP + structure-download layer.

All network access funnels through here so it can be cached on disk (deterministic
re-runs, friendlier to rate limits) and so every downloaded artifact lands inside
the project sandbox (./cache by default), never in ~ or /tmp.
"""
from __future__ import annotations

import hashlib
import json
import os
import time
from pathlib import Path
from typing import Any, Optional

import requests

# Sandbox-local cache root. Overridable via env for tests, but defaults inside
# the project directory per the isolation rules.
_DEFAULT_CACHE = Path(__file__).resolve().parents[2] / "cache"
CACHE_DIR = Path(os.environ.get("TRIMPROT_CACHE", _DEFAULT_CACHE))

_SESSION: Optional[requests.Session] = None
USER_AGENT = "trimprot/0.1 (deterministic target-prep pipeline)"
TIMEOUT = 60


def cache_dir() -> Path:
    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    return CACHE_DIR


def _session() -> requests.Session:
    global _SESSION
    if _SESSION is None:
        s = requests.Session()
        s.headers.update({"User-Agent": USER_AGENT})
        _SESSION = s
    return _SESSION


def _key(*parts: str) -> str:
    h = hashlib.sha256("||".join(parts).encode()).hexdigest()[:24]
    return h


def _request(method: str, url: str, *, params=None, json_body=None,
             headers=None, max_retries: int = 4) -> requests.Response:
    last_exc = None
    for attempt in range(max_retries):
        try:
            resp = _session().request(
                method, url, params=params, json=json_body,
                headers=headers, timeout=TIMEOUT,
            )
            # Retry transient server errors / rate limits.
            if resp.status_code in (429, 500, 502, 503, 504):
                raise requests.HTTPError(f"{resp.status_code} for {url}")
            return resp
        except (requests.RequestException,) as exc:  # noqa: PERF203
            last_exc = exc
            sleep = min(2 ** attempt, 8)
            time.sleep(sleep)
    raise RuntimeError(f"request failed after {max_retries} retries: {url}: {last_exc}")


def get_json(url: str, params: Optional[dict] = None,
             *, cache: bool = True, allow_404: bool = False) -> Any:
    """GET returning parsed JSON, cached on disk by url+params."""
    keystr = url + "?" + json.dumps(params or {}, sort_keys=True)
    path = cache_dir() / "http" / f"{_key('GET', keystr)}.json"
    if cache and path.exists():
        return json.loads(path.read_text())
    resp = _request("GET", url, params=params)
    if allow_404 and resp.status_code == 404:
        data = None
    else:
        resp.raise_for_status()
        data = resp.json()
    if cache:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(data))
    return data


def post_json(url: str, payload: dict, *, cache: bool = True) -> Any:
    """POST JSON (used for RCSB GraphQL), cached by url+payload."""
    keystr = url + "::" + json.dumps(payload, sort_keys=True)
    path = cache_dir() / "http" / f"{_key('POST', keystr)}.json"
    if cache and path.exists():
        return json.loads(path.read_text())
    resp = _request("POST", url, json_body=payload,
                    headers={"Content-Type": "application/json"})
    resp.raise_for_status()
    data = resp.json()
    if cache:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(data))
    return data


def download_text(url: str, dest: Path, *, cache: bool = True) -> Path:
    """Download a text file (e.g. mmCIF) to dest, cached."""
    dest.parent.mkdir(parents=True, exist_ok=True)
    if cache and dest.exists() and dest.stat().st_size > 0:
        return dest
    resp = _request("GET", url)
    resp.raise_for_status()
    dest.write_bytes(resp.content)
    return dest


def structure_path(pdb_id: str) -> Path:
    return cache_dir() / "structures" / f"{pdb_id.lower()}.cif"


def download_cif(pdb_id: str) -> Path:
    """Download the asymmetric-unit mmCIF from RCSB into the sandbox cache."""
    pdb_id = pdb_id.lower()
    dest = structure_path(pdb_id)
    url = f"https://files.rcsb.org/download/{pdb_id}.cif"
    return download_text(url, dest)
