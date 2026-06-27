"""
tests/test_alphafold_pipeline.py
=================================
End-to-end pipeline integration tests for proteins with no crystal structure.

All three targets are human, reviewed Swiss-Prot entries with a transmembrane
region and annotated extracellular domain, but no PDB cross-reference — they
exercise the AlphaFold fallback path in run_egfr.run().

Run with:
    pytest -v -m network tests/test_alphafold_pipeline.py
"""

from __future__ import annotations

import json
import sys
from pathlib import Path
from unittest.mock import patch

import pytest

sys.path.insert(0, str(Path(__file__).parent.parent))

import run_egfr


# ── Target table ─────────────────────────────────────────────────────────────
# (accession, gene, ecd_start, ecd_end, n_tm, description)
TARGETS = [
    ("A0A1B0GTW7", "CIROP",  21,  735, 1,
     "Ciliated left-right organizer metallopeptidase; large single-pass ECD"),
    ("B6A8C7",     "TARM1",  17,  236, 1,
     "T-cell-interacting activating receptor on myeloid cells; immune target"),
    ("A6BM72",     "MEGF11", 20,  848, 1,
     "Multiple EGF-like domains protein 11; large ECD, multiple EGF repeats"),
]


@pytest.mark.network
@pytest.mark.parametrize("accession,gene,ecd_start,ecd_end,n_tm,description", TARGETS,
                         ids=[t[1] for t in TARGETS])
def test_alphafold_fallback_used(accession, gene, ecd_start, ecd_end, n_tm, description, tmp_path):
    """Pipeline must select AlphaFold because no PDB structure exists."""
    with patch("run_egfr.OUT_DIR", str(tmp_path)):
        result = run_egfr.run(accession)

    assert result.get("alphafold_fallback") is True, (
        f"{gene} ({accession}): expected AlphaFold fallback but got a crystal structure path"
    )


@pytest.mark.network
@pytest.mark.parametrize("accession,gene,ecd_start,ecd_end,n_tm,description", TARGETS,
                         ids=[t[1] for t in TARGETS])
def test_ecd_boundaries_match_uniprot(accession, gene, ecd_start, ecd_end, n_tm, description, tmp_path):
    """Extracellular domain boundaries must match UniProt topological annotation."""
    with patch("run_egfr.OUT_DIR", str(tmp_path)):
        result = run_egfr.run(accession)

    ecd = result["extracellular_domain"]
    assert ecd["start"] == ecd_start, f"{gene}: ECD start {ecd['start']} != expected {ecd_start}"
    assert ecd["end"]   == ecd_end,   f"{gene}: ECD end {ecd['end']} != expected {ecd_end}"
    assert ecd["trimming_decision"]["trimmed"] is True


@pytest.mark.network
@pytest.mark.parametrize("accession,gene,ecd_start,ecd_end,n_tm,description", TARGETS,
                         ids=[t[1] for t in TARGETS])
def test_output_files_written(accession, gene, ecd_start, ecd_end, n_tm, description, tmp_path):
    """Trimmed PDB, full PDB, CIF, and summary JSON must all be written."""
    with patch("run_egfr.OUT_DIR", str(tmp_path)):
        run_egfr.run(accession)

    assert (tmp_path / f"AF-{accession}_ECD_trimmed.pdb").exists(), f"{gene}: missing trimmed PDB"
    assert (tmp_path / f"AF-{accession}_full.pdb").exists(),        f"{gene}: missing full PDB"
    assert (tmp_path / f"AF-{accession}_full.cif").exists(),        f"{gene}: missing full CIF"
    assert (tmp_path / f"AF-{accession}_summary.json").exists(),    f"{gene}: missing summary JSON"


@pytest.mark.network
@pytest.mark.parametrize("accession,gene,ecd_start,ecd_end,n_tm,description", TARGETS,
                         ids=[t[1] for t in TARGETS])
def test_trimming_stats_are_sane(accession, gene, ecd_start, ecd_end, n_tm, description, tmp_path):
    """Residues kept must be > 0 and ≤ ECD length; trimming must remove some residues."""
    with patch("run_egfr.OUT_DIR", str(tmp_path)):
        result = run_egfr.run(accession)

    trimming = result["trimming"]
    ecd_length = ecd_end - ecd_start + 1
    assert trimming["residues_kept"] > 0, f"{gene}: no residues kept after trim"
    assert trimming["residues_kept"] <= ecd_length, (
        f"{gene}: kept {trimming['residues_kept']} residues but ECD is only {ecd_length} aa"
    )
    assert trimming["residues_trimmed_away"] > 0, f"{gene}: expected TM/intracellular residues to be trimmed"


@pytest.mark.network
@pytest.mark.parametrize("accession,gene,ecd_start,ecd_end,n_tm,description", TARGETS,
                         ids=[t[1] for t in TARGETS])
def test_plddt_summary_fields_present(accession, gene, ecd_start, ecd_end, n_tm, description, tmp_path):
    """Summary must include all AlphaFold pLDDT metadata fields."""
    with patch("run_egfr.OUT_DIR", str(tmp_path)):
        result = run_egfr.run(accession)

    required = {
        "alphafold_version",
        "alphafold_model_url",
        "alphafold_ecd_mean_plddt",
        "alphafold_ecd_n_low_plddt",
        "alphafold_ecd_n_unobserved",
        "alphafold_low_plddt_ranges",
        "alphafold_unobs_ranges",
    }
    missing = required - result.keys()
    assert not missing, f"{gene}: summary missing keys: {missing}"

    assert isinstance(result["alphafold_ecd_mean_plddt"], float)
    assert 0 < result["alphafold_ecd_mean_plddt"] <= 100
    assert isinstance(result["alphafold_low_plddt_ranges"], list)
    assert isinstance(result["alphafold_unobs_ranges"], list)


@pytest.mark.network
@pytest.mark.parametrize("accession,gene,ecd_start,ecd_end,n_tm,description", TARGETS,
                         ids=[t[1] for t in TARGETS])
def test_hotspots_found(accession, gene, ecd_start, ecd_end, n_tm, description, tmp_path):
    """Surface-exposure heuristic must find at least some hotspot candidates."""
    with patch("run_egfr.OUT_DIR", str(tmp_path)):
        result = run_egfr.run(accession)

    hotspots = result["hotspots"]
    assert hotspots["count"] > 0, f"{gene}: no hotspot candidates found"
    assert "alphafold" in hotspots["source"].lower(), (
        f"{gene}: hotspot source should mention AlphaFold, got: {hotspots['source']}"
    )
    assert hotspots["partner_chains"] == [], f"{gene}: AF monomers should have no partner chains"

    # All hotspot UniProt positions must fall inside the ECD
    for pos in hotspots["candidate_residues_unp"]:
        assert ecd_start <= pos <= ecd_end, (
            f"{gene}: hotspot at UniProt pos {pos} is outside ECD {ecd_start}-{ecd_end}"
        )


@pytest.mark.network
@pytest.mark.parametrize("accession,gene,ecd_start,ecd_end,n_tm,description", TARGETS,
                         ids=[t[1] for t in TARGETS])
def test_viewer_block_references_af_files(accession, gene, ecd_start, ecd_end, n_tm, description, tmp_path):
    """viewer dict must reference AF-prefixed files and use chain A."""
    with patch("run_egfr.OUT_DIR", str(tmp_path)):
        result = run_egfr.run(accession)

    viewer = result["viewer"]
    assert viewer["pdb_id"] == f"AF-{accession}"
    assert viewer["auth_chain"] == "A"
    assert viewer["partner_chains"] == []
    assert viewer["original_file"].startswith("AF-")
    assert viewer["trimmed_file"].startswith("AF-")
    assert viewer["trimmed_file"].endswith("_ECD_trimmed.pdb")


@pytest.mark.network
@pytest.mark.parametrize("accession,gene,ecd_start,ecd_end,n_tm,description", TARGETS,
                         ids=[t[1] for t in TARGETS])
def test_summary_json_is_valid(accession, gene, ecd_start, ecd_end, n_tm, description, tmp_path):
    """The written summary.json must be valid JSON and round-trip cleanly."""
    with patch("run_egfr.OUT_DIR", str(tmp_path)):
        run_egfr.run(accession)

    summary_path = tmp_path / f"AF-{accession}_summary.json"
    with open(summary_path) as f:
        data = json.load(f)

    assert data["target"]["uniprot_accession"] == accession
    assert data["alphafold_fallback"] is True


@pytest.mark.network
@pytest.mark.parametrize("accession,gene,ecd_start,ecd_end,n_tm,description", TARGETS,
                         ids=[t[1] for t in TARGETS])
def test_avoid_residues_populated(accession, gene, ecd_start, ecd_end, n_tm, description, tmp_path):
    """avoid_residues section must be present; low_plddt key included for AF path."""
    with patch("run_egfr.OUT_DIR", str(tmp_path)):
        result = run_egfr.run(accession)

    avoid = result["avoid_residues"]
    assert "glycosylation_sites" in avoid
    assert "disulfide_cysteines" in avoid
    assert "other_ptms" in avoid
    assert "missing_unresolved" in avoid
    assert "low_plddt" in avoid, f"{gene}: AlphaFold path must include low_plddt avoid category"
    assert "counts" in avoid

    # Counts dict must match actual list lengths
    for key, lst in avoid.items():
        if key == "counts":
            continue
        assert avoid["counts"][key] == len(lst), (
            f"{gene}: counts[{key!r}]={avoid['counts'][key]} but len={len(lst)}"
        )
