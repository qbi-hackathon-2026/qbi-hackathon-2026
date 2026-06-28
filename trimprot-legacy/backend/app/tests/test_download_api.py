"""
tests/test_download_api.py
==========================
Tests for the structured download endpoints added to api.py:

  GET /api/downloads/{accession}          — metadata listing available files
  GET /api/download/{accession}/{type}    — actual file download (type: trimmed | full | cif)

All tests use a pre-built fake pipeline result injected into the server's
_cache — no network access, no real files on disk (file content is written
to a tmp_path for the download tests that need the file to exist).

Run with:
  pytest -v tests/test_download_api.py       # no network required
"""

from __future__ import annotations

import os
import sys
from pathlib import Path
from unittest.mock import patch

import pytest
from fastapi.testclient import TestClient

sys.path.insert(0, str(Path(__file__).parent.parent))

import api as api_module
from api import app

client = TestClient(app, raise_server_exceptions=False)


# ── Shared fake pipeline results ──────────────────────────────────────────────

def _af_result(accession: str = "B6A8C7") -> dict:
    """Minimal AlphaFold-path pipeline result (mirrors what run_egfr.run produces)."""
    return {
        "alphafold_fallback": True,
        "target": {"uniprot_accession": accession, "uniprot_id": "TARM1_HUMAN", "isoform": "canonical"},
        "viewer": {
            "pdb_id": f"AF-{accession}",
            "auth_chain": "A",
            "partner_chains": [],
            "original_file":     f"AF-{accession}_full.pdb",
            "original_cif_file": f"AF-{accession}_full.cif",
            "trimmed_file":      f"AF-{accession}_ECD_trimmed.pdb",
            "hotspot_auth_residues": [50, 60, 70],
            "avoid_glycosylation_auth_residues": [],
            "avoid_disulfide_auth_residues": [],
            "avoid_other_ptm_auth_residues": [],
            "avoid_low_plddt_auth_residues": [20, 21],
            "missing_auth_residues_note": "pLDDT ≤ 50 residues treated as unobserved",
        },
    }


def _pdb_result(accession: str = "P00533", pdb_id: str = "1IVO") -> dict:
    """Minimal crystal-structure-path pipeline result."""
    return {
        "alphafold_fallback": False,
        "target": {"uniprot_accession": accession, "uniprot_id": "EGFR_HUMAN", "isoform": "canonical"},
        "viewer": {
            "pdb_id": pdb_id,
            "auth_chain": "A",
            "partner_chains": ["B"],
            "original_file":     f"{pdb_id}_full.pdb",
            "original_cif_file": f"{pdb_id}_full.cif",
            "trimmed_file":      f"{pdb_id}_ECD_trimmed.pdb",
            "hotspot_auth_residues": [384, 408],
            "avoid_glycosylation_auth_residues": [128],
            "avoid_disulfide_auth_residues": [],
            "avoid_other_ptm_auth_residues": [],
            "missing_auth_residues_note": "missing residues have no atoms",
        },
    }


# ═════════════════════════════════════════════════════════════════════════════
# /api/downloads/{accession}  —  metadata endpoint
# ═════════════════════════════════════════════════════════════════════════════

class TestListDownloads:

    def test_404_when_not_run(self):
        r = client.get("/api/downloads/NOTRUN")
        assert r.status_code == 404
        assert "not yet run" in r.json()["detail"].lower()

    def test_af_metadata_shape(self):
        with patch.dict(api_module._cache, {"B6A8C7": _af_result()}):
            r = client.get("/api/downloads/B6A8C7")
        assert r.status_code == 200
        data = r.json()
        assert data["accession"] == "B6A8C7"
        assert data["alphafold_fallback"] is True
        assert data["source"] == "AlphaFold"
        assert set(data["files"].keys()) == {"trimmed", "full", "cif"}

    def test_pdb_metadata_shape(self):
        with patch.dict(api_module._cache, {"P00533": _pdb_result()}):
            r = client.get("/api/downloads/P00533")
        assert r.status_code == 200
        data = r.json()
        assert data["alphafold_fallback"] is False
        assert "1IVO" in data["source"]

    def test_af_urls_reference_accession(self):
        with patch.dict(api_module._cache, {"B6A8C7": _af_result()}):
            r = client.get("/api/downloads/B6A8C7")
        files = r.json()["files"]
        for key, entry in files.items():
            assert "B6A8C7" in entry["url"], f"{key} URL doesn't include accession"
            assert entry["url"].startswith("/api/download/")

    def test_af_filenames_contain_accession(self):
        with patch.dict(api_module._cache, {"B6A8C7": _af_result()}):
            r = client.get("/api/downloads/B6A8C7")
        files = r.json()["files"]
        assert "B6A8C7" in files["trimmed"]["filename"]
        assert "B6A8C7" in files["full"]["filename"]
        assert "B6A8C7" in files["cif"]["filename"]

    def test_trimmed_description_mentions_ecd(self):
        with patch.dict(api_module._cache, {"B6A8C7": _af_result()}):
            r = client.get("/api/downloads/B6A8C7")
        desc = r.json()["files"]["trimmed"]["description"].lower()
        assert "trimmed" in desc or "extracellular" in desc

    def test_formats_are_labeled(self):
        with patch.dict(api_module._cache, {"B6A8C7": _af_result()}):
            r = client.get("/api/downloads/B6A8C7")
        files = r.json()["files"]
        assert files["trimmed"]["format"] == "PDB"
        assert files["full"]["format"] == "PDB"
        assert files["cif"]["format"] == "mmCIF"


# ═════════════════════════════════════════════════════════════════════════════
# /api/download/{accession}/{type}  —  file download endpoint
# ═════════════════════════════════════════════════════════════════════════════

class TestDownloadStructure:

    def test_404_when_accession_not_run(self):
        r = client.get("/api/download/NOTRUN/trimmed")
        assert r.status_code == 404
        assert "not yet run" in r.json()["detail"].lower()

    def test_422_for_invalid_file_type(self):
        with patch.dict(api_module._cache, {"B6A8C7": _af_result()}):
            r = client.get("/api/download/B6A8C7/invalid")
        assert r.status_code == 422

    def test_404_when_file_missing_from_disk(self):
        with patch.dict(api_module._cache, {"B6A8C7": _af_result()}):
            with patch("api.OUT_DIR", "/nonexistent/path"):
                r = client.get("/api/download/B6A8C7/trimmed")
        assert r.status_code == 404
        assert "not found on disk" in r.json()["detail"].lower()

    def _write_fake_files(self, tmp_path: Path, accession: str, is_af: bool) -> dict:
        """Write placeholder files and patch OUT_DIR; return the result dict."""
        result = _af_result(accession) if is_af else _pdb_result()
        viewer = result["viewer"]
        for key in ("original_file", "original_cif_file", "trimmed_file"):
            (tmp_path / viewer[key]).write_bytes(b"FAKE PDB/CIF CONTENT")
        return result

    def test_af_trimmed_download_returns_200(self, tmp_path):
        result = self._write_fake_files(tmp_path, "B6A8C7", is_af=True)
        with patch.dict(api_module._cache, {"B6A8C7": result}):
            with patch("api.OUT_DIR", str(tmp_path)):
                r = client.get("/api/download/B6A8C7/trimmed")
        assert r.status_code == 200

    def test_af_full_download_returns_200(self, tmp_path):
        result = self._write_fake_files(tmp_path, "B6A8C7", is_af=True)
        with patch.dict(api_module._cache, {"B6A8C7": result}):
            with patch("api.OUT_DIR", str(tmp_path)):
                r = client.get("/api/download/B6A8C7/full")
        assert r.status_code == 200

    def test_af_cif_download_returns_200(self, tmp_path):
        result = self._write_fake_files(tmp_path, "B6A8C7", is_af=True)
        with patch.dict(api_module._cache, {"B6A8C7": result}):
            with patch("api.OUT_DIR", str(tmp_path)):
                r = client.get("/api/download/B6A8C7/cif")
        assert r.status_code == 200

    def test_pdb_trimmed_download_returns_200(self, tmp_path):
        result = self._write_fake_files(tmp_path, "P00533", is_af=False)
        with patch.dict(api_module._cache, {"P00533": result}):
            with patch("api.OUT_DIR", str(tmp_path)):
                r = client.get("/api/download/P00533/trimmed")
        assert r.status_code == 200

    def test_trimmed_pdb_content_type(self, tmp_path):
        result = self._write_fake_files(tmp_path, "B6A8C7", is_af=True)
        with patch.dict(api_module._cache, {"B6A8C7": result}):
            with patch("api.OUT_DIR", str(tmp_path)):
                r = client.get("/api/download/B6A8C7/trimmed")
        assert "chemical/x-pdb" in r.headers["content-type"]

    def test_cif_content_type(self, tmp_path):
        result = self._write_fake_files(tmp_path, "B6A8C7", is_af=True)
        with patch.dict(api_module._cache, {"B6A8C7": result}):
            with patch("api.OUT_DIR", str(tmp_path)):
                r = client.get("/api/download/B6A8C7/cif")
        assert "chemical/x-mmcif" in r.headers["content-type"]

    def test_content_disposition_filename_is_set(self, tmp_path):
        result = self._write_fake_files(tmp_path, "B6A8C7", is_af=True)
        with patch.dict(api_module._cache, {"B6A8C7": result}):
            with patch("api.OUT_DIR", str(tmp_path)):
                r = client.get("/api/download/B6A8C7/trimmed")
        disposition = r.headers.get("content-disposition", "")
        assert "AF-B6A8C7" in disposition
        assert "ECD_trimmed.pdb" in disposition

    def test_trimmed_url_from_metadata_resolves(self, tmp_path):
        """The URL returned by /api/downloads must actually work."""
        result = self._write_fake_files(tmp_path, "B6A8C7", is_af=True)
        with patch.dict(api_module._cache, {"B6A8C7": result}):
            with patch("api.OUT_DIR", str(tmp_path)):
                meta = client.get("/api/downloads/B6A8C7").json()
                trimmed_url = meta["files"]["trimmed"]["url"]
                r = client.get(trimmed_url)
        assert r.status_code == 200

    def test_full_url_from_metadata_resolves(self, tmp_path):
        result = self._write_fake_files(tmp_path, "B6A8C7", is_af=True)
        with patch.dict(api_module._cache, {"B6A8C7": result}):
            with patch("api.OUT_DIR", str(tmp_path)):
                meta = client.get("/api/downloads/B6A8C7").json()
                full_url = meta["files"]["full"]["url"]
                r = client.get(full_url)
        assert r.status_code == 200
