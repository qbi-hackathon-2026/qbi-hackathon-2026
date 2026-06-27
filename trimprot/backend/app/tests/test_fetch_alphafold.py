"""
tests/test_fetch_alphafold.py
=============================
Test suite for fetch_alphafold.py.

Unit tests use synthetic in-memory CIFs — no network access required.
Integration tests (marked `network`) hit the EBI AlphaFold API and are
skipped when -m "not network" is passed or the API is unreachable.

pytest -v tests/test_fetch_alphafold.py                   # unit tests only
pytest -v tests/test_fetch_alphafold.py -m network        # integration tests
pytest -v tests/test_fetch_alphafold.py -m "not network"  # skip network tests
"""

from __future__ import annotations

import os
import sys
import tempfile
from pathlib import Path
from typing import Optional
from unittest.mock import patch

import gemmi
import pytest

# Allow direct imports of app-level modules without installing them
sys.path.insert(0, str(Path(__file__).parent.parent))

from fetch_alphafold import (
    PLDDT_LOW_THRESHOLD,
    PLDDT_UNOBS_THRESHOLD,
    AvoidEntry,
    AlphaFoldResult,
    ResidueConfidence,
    _avoid_reason,
    _build_avoid_contributions,
    _confidence_band,
    _ecd_set_from_ranges,
    _extract_residue_confidences,
    _merge_ranges,
    _residue_plddt,
    get_alphafold_structure,
)


# ═════════════════════════════════════════════════════════════════════════════
# Helpers — synthetic CIF builder
# ═════════════════════════════════════════════════════════════════════════════

def _make_cif(residues: list[tuple[int, str, float, str]]) -> gemmi.Structure:
    """
    Build a gemmi Structure from a list of (auth_seq_id, res_name, plddt, chain).
    Only Cα atoms are added (sufficient for pLDDT extraction and tests).
    """
    cif_lines = [
        "data_TEST",
        "_entry.id TEST",
        "loop_",
        "_atom_site.group_PDB",
        "_atom_site.id",
        "_atom_site.type_symbol",
        "_atom_site.label_atom_id",
        "_atom_site.label_alt_id",
        "_atom_site.label_comp_id",
        "_atom_site.label_asym_id",
        "_atom_site.label_entity_id",
        "_atom_site.label_seq_id",
        "_atom_site.pdbx_PDB_ins_code",
        "_atom_site.Cartn_x",
        "_atom_site.Cartn_y",
        "_atom_site.Cartn_z",
        "_atom_site.occupancy",
        "_atom_site.B_iso_or_equiv",
        "_atom_site.auth_seq_id",
        "_atom_site.auth_comp_id",
        "_atom_site.auth_asym_id",
        "_atom_site.auth_atom_id",
        "_atom_site.pdbx_PDB_model_num",
    ]
    for atom_id, (seq_id, res_name, plddt, chain) in enumerate(residues, 1):
        label_seq = atom_id
        cif_lines.append(
            f"ATOM {atom_id} C CA . {res_name} {chain} 1 {label_seq} ? "
            f"{float(atom_id):.3f} {float(atom_id*2):.3f} {float(atom_id*3):.3f} "
            f"1.00 {plddt:.2f} {seq_id} {res_name} {chain} CA 1"
        )
    cif_lines.append("#")

    with tempfile.NamedTemporaryFile(suffix=".cif", mode="w", delete=False) as f:
        f.write("\n".join(cif_lines))
        path = f.name

    try:
        return gemmi.read_structure(path)
    finally:
        os.unlink(path)


# ── Shared fixtures ───────────────────────────────────────────────────────────

EGFR_RESIDUES = [
    # (auth_seq_id, res_name, plddt, chain)
    (25,  "GLN", 88.2, "A"),   # confident — ECD start (mature chain)
    (100, "ARG", 91.5, "A"),   # very high
    (200, "TYR", 85.0, "A"),   # confident
    (384, "ARG", 87.3, "A"),   # confident — known EGFR hotspot
    (408, "ASP", 90.1, "A"),   # very high — known hotspot
    (620, "SER", 68.4, "A"),   # LOW  50 < pLDDT < 70  → avoid_low
    (635, "PRO", 45.2, "A"),   # VERY LOW ≤ 50          → avoid_unobs
    (645, "VAL", 48.8, "A"),   # VERY LOW               → avoid_unobs
]

EGFR_ECD_RANGES = [(25, 645)]


# ═════════════════════════════════════════════════════════════════════════════
# Unit tests — pure functions (no network)
# ═════════════════════════════════════════════════════════════════════════════

class TestConfidenceBand:
    def test_very_high(self):
        assert _confidence_band(91.0) == "very_high"

    def test_confident(self):
        assert _confidence_band(75.0) == "confident"
        assert _confidence_band(70.01) == "confident"

    def test_low(self):
        assert _confidence_band(69.9) == "low"
        assert _confidence_band(50.01) == "low"

    def test_very_low(self):
        assert _confidence_band(50.0) == "very_low"
        assert _confidence_band(10.0) == "very_low"

    def test_boundary_70_is_confident(self):
        # pLDDT == 70.0 exactly → "confident" (threshold is strict <)
        assert _confidence_band(70.0) == "confident"

    def test_boundary_50_is_very_low(self):
        # pLDDT == 50.0 exactly → "very_low" (threshold is ≤)
        assert _confidence_band(50.0) == "very_low"


class TestAvoidReason:
    def test_no_reason_above_threshold(self):
        assert _avoid_reason(70.0) is None
        assert _avoid_reason(90.0) is None

    def test_low_plddt_reason(self):
        reason = _avoid_reason(65.0)
        assert reason is not None
        assert "low_plddt" in reason
        assert "70" in reason

    def test_unobserved_reason(self):
        reason = _avoid_reason(45.0)
        assert reason is not None
        assert "unobserved" in reason
        assert "50" in reason

    def test_exactly_at_unobs_threshold(self):
        reason = _avoid_reason(50.0)
        assert reason is not None
        assert "unobserved" in reason


class TestMergeRanges:
    def test_empty(self):
        assert _merge_ranges([]) == []

    def test_single(self):
        assert _merge_ranges([42]) == [(42, 42)]

    def test_contiguous(self):
        assert _merge_ranges([1, 2, 3]) == [(1, 3)]

    def test_two_gaps(self):
        assert _merge_ranges([10, 11, 20, 21, 22]) == [(10, 11), (20, 22)]

    def test_duplicates_ignored(self):
        assert _merge_ranges([5, 5, 6, 6]) == [(5, 6)]

    def test_unsorted_input(self):
        assert _merge_ranges([3, 1, 2]) == [(1, 3)]


class TestEcdSetFromRanges:
    def test_single_range(self):
        assert _ecd_set_from_ranges([(10, 13)]) == frozenset({10, 11, 12, 13})

    def test_multiple_ranges(self):
        assert _ecd_set_from_ranges([(1, 3), (10, 11)]) == frozenset({1, 2, 3, 10, 11})

    def test_empty(self):
        assert _ecd_set_from_ranges([]) == frozenset()

    def test_single_residue_range(self):
        assert _ecd_set_from_ranges([(42, 42)]) == frozenset({42})


class TestResiduePlddt:
    def test_ca_plddt_extracted(self):
        struct = _make_cif([(25, "GLN", 88.2, "A")])
        residue = list(struct[0]["A"])[0]
        assert abs(_residue_plddt(residue) - 88.2) < 0.01

    def test_very_low_plddt(self):
        struct = _make_cif([(635, "PRO", 45.2, "A")])
        residue = list(struct[0]["A"])[0]
        assert _residue_plddt(residue) < PLDDT_UNOBS_THRESHOLD


class TestExtractResidueConfidences:
    def setup_method(self):
        self.struct = _make_cif(EGFR_RESIDUES)
        self.ecd_set = _ecd_set_from_ranges(EGFR_ECD_RANGES)
        self.records = _extract_residue_confidences(self.struct, self.ecd_set, chain_id="A")

    def test_record_count_matches_residues(self):
        assert len(self.records) == len(EGFR_RESIDUES)

    def test_auth_numbering_preserved(self):
        """Rule 1 — auth_seq_id must match the input auth numbers exactly."""
        input_ids = [r[0] for r in EGFR_RESIDUES]
        output_ids = [r.auth_seq_id for r in self.records]
        assert output_ids == input_ids

    def test_no_insertion_codes(self):
        for rec in self.records:
            assert rec.icode == ""

    def test_all_ecd_residues_flagged(self):
        for rec in self.records:
            assert rec.in_ecd is True

    def test_high_confidence_residues_have_no_avoid_reason(self):
        confident_ids = {25, 100, 200, 384, 408}
        for rec in self.records:
            if rec.auth_seq_id in confident_ids:
                assert rec.avoid_reason is None

    def test_low_plddt_residue_has_avoid_reason(self):
        rec_620 = next(r for r in self.records if r.auth_seq_id == 620)
        assert rec_620.plddt == pytest.approx(68.4, abs=0.1)
        assert rec_620.avoid_reason is not None
        assert "low_plddt" in rec_620.avoid_reason
        assert "unobserved" not in rec_620.avoid_reason

    def test_very_low_plddt_residues_treated_as_unobserved(self):
        for seq_id in (635, 645):
            rec = next(r for r in self.records if r.auth_seq_id == seq_id)
            assert rec.plddt < PLDDT_UNOBS_THRESHOLD
            assert rec.avoid_reason is not None
            assert "unobserved" in rec.avoid_reason

    def test_plddt_bands(self):
        band_map = {r.auth_seq_id: r.confidence_band for r in self.records}
        assert band_map[25]  == "confident"
        assert band_map[100] == "very_high"
        assert band_map[620] == "low"
        assert band_map[635] == "very_low"

    def test_wrong_chain_raises(self):
        with pytest.raises(RuntimeError, match="Chain 'B' not found"):
            _extract_residue_confidences(self.struct, self.ecd_set, chain_id="B")


class TestBuildAvoidContributions:
    def setup_method(self):
        struct = _make_cif(EGFR_RESIDUES)
        ecd_set = _ecd_set_from_ranges(EGFR_ECD_RANGES)
        self.records = _extract_residue_confidences(struct, ecd_set)
        self.entries = _build_avoid_contributions(self.records, chain_id="A", accession="P00533")
        self.avoid_ids = {e.auth_seq_id for e in self.entries}

    def test_confident_residues_not_in_avoid_set(self):
        for seq_id in (25, 100, 200, 384, 408):
            assert seq_id not in self.avoid_ids

    def test_low_plddt_in_avoid_set(self):
        assert 620 in self.avoid_ids

    def test_very_low_plddt_in_avoid_set(self):
        assert 635 in self.avoid_ids
        assert 645 in self.avoid_ids

    def test_reason_strings_are_correct(self):
        for entry in self.entries:
            if entry.auth_seq_id == 620:
                assert "low_plddt" in entry.reason
            elif entry.auth_seq_id in (635, 645):
                assert "unobserved" in entry.reason

    def test_uniprot_residue_equals_auth_for_alphafold(self):
        """KEYSTONE: AlphaFold auth_seq_id == UniProt residue number."""
        for entry in self.entries:
            assert entry.auth_seq_id == entry.uniprot_residue

    def test_avoid_entries_have_plddt_populated(self):
        for entry in self.entries:
            assert entry.plddt is not None
            assert 0.0 <= entry.plddt <= 100.0

    def test_non_ecd_residues_never_appear(self):
        for entry in self.entries:
            rec = next(r for r in self.records if r.auth_seq_id == entry.auth_seq_id)
            assert rec.in_ecd


class TestAvoidEntryKey:
    def test_key_tuple_structure(self):
        e = AvoidEntry(chain="A", auth_seq_id=620, icode="", reason="test")
        assert e.key() == ("A", 620, "")

    def test_keys_usable_in_set(self):
        e1 = AvoidEntry(chain="A", auth_seq_id=620, icode="", reason="r1")
        e2 = AvoidEntry(chain="A", auth_seq_id=620, icode="", reason="r2")
        e3 = AvoidEntry(chain="A", auth_seq_id=621, icode="", reason="r3")
        keys = {e1.key(), e2.key(), e3.key()}
        assert len(keys) == 2


class TestPartialEcdCoverage:
    def test_non_ecd_low_plddt_not_added(self):
        residues = [
            (50,  "ALA", 30.0, "A"),   # very low pLDDT, NOT in ECD
            (100, "GLY", 88.0, "A"),   # good pLDDT, in ECD
            (150, "SER", 65.0, "A"),   # low pLDDT, in ECD → should appear
            (200, "ARG", 92.0, "A"),   # very high, in ECD
        ]
        struct  = _make_cif(residues)
        ecd_set = _ecd_set_from_ranges([(100, 200)])
        records = _extract_residue_confidences(struct, ecd_set)
        entries = _build_avoid_contributions(records, "A", "TEST")

        avoid_ids = {e.auth_seq_id for e in entries}
        assert 50  not in avoid_ids
        assert 100 not in avoid_ids
        assert 150 in avoid_ids
        assert 200 not in avoid_ids

    def test_type_ii_ecd_c_terminal(self):
        """Type II proteins (e.g. CD38) have their ECD at the C-terminus."""
        residues = [
            (1,   "MET", 88.0, "A"),
            (10,  "ALA", 82.0, "A"),
            (43,  "GLU", 91.0, "A"),   # ECD start (type II)
            (100, "LYS", 48.0, "A"),   # ECD, very low pLDDT → avoid
            (200, "TYR", 85.0, "A"),
        ]
        struct  = _make_cif(residues)
        ecd_set = _ecd_set_from_ranges([(43, 300)])
        records = _extract_residue_confidences(struct, ecd_set)
        entries = _build_avoid_contributions(records, "A", "P28907")

        avoid_ids = {e.auth_seq_id for e in entries}
        assert 1   not in avoid_ids
        assert 10  not in avoid_ids
        assert 43  not in avoid_ids
        assert 100 in avoid_ids
        assert 200 not in avoid_ids


class TestSummaryFields:
    def _make_result(self) -> AlphaFoldResult:
        struct = _make_cif(EGFR_RESIDUES)
        ecd_set = _ecd_set_from_ranges(EGFR_ECD_RANGES)
        records = _extract_residue_confidences(struct, ecd_set)
        entries = _build_avoid_contributions(records, "A", "P00533")

        ecd_plddts = [r.plddt for r in records if r.in_ecd]
        mean = sum(ecd_plddts) / len(ecd_plddts)

        low_pos  = [r.auth_seq_id for r in records
                    if r.in_ecd and PLDDT_UNOBS_THRESHOLD < r.plddt < PLDDT_LOW_THRESHOLD]
        unob_pos = [r.auth_seq_id for r in records
                    if r.in_ecd and r.plddt <= PLDDT_UNOBS_THRESHOLD]

        return AlphaFoldResult(
            accession="P00533",
            version=4,
            structure=struct,
            model_url="https://example.com/AF-P00533.cif",
            chain_id="A",
            residues=records,
            avoid_contributions=entries,
            ecd_mean_plddt=mean,
            ecd_n_low_plddt=len(low_pos),
            ecd_n_unobserved=len(unob_pos),
            low_plddt_ranges=_merge_ranges(low_pos),
            unobs_ranges=_merge_ranges(unob_pos),
        )

    def test_required_keys_present(self):
        fields = self._make_result().summary_fields()
        required = {
            "alphafold_fallback",
            "alphafold_version",
            "alphafold_model_url",
            "alphafold_ecd_mean_plddt",
            "alphafold_ecd_n_low_plddt",
            "alphafold_ecd_n_unobserved",
            "alphafold_low_plddt_ranges",
            "alphafold_unobs_ranges",
        }
        assert required.issubset(fields.keys())

    def test_fallback_flag_always_true(self):
        assert self._make_result().summary_fields()["alphafold_fallback"] is True

    def test_ranges_are_lists_of_lists(self):
        fields = self._make_result().summary_fields()
        for key in ("alphafold_low_plddt_ranges", "alphafold_unobs_ranges"):
            for item in fields[key]:
                assert isinstance(item, list)
                assert len(item) == 2

    def test_mean_plddt_is_rounded(self):
        val = self._make_result().summary_fields()["alphafold_ecd_mean_plddt"]
        assert isinstance(val, float)
        assert round(val, 2) == val


class TestGetAlphaFoldStructureValidation:
    def test_empty_ecd_ranges_raises(self):
        with pytest.raises(RuntimeError, match="ecd_ranges must not be empty"):
            get_alphafold_structure("P00533", ecd_ranges=[])

    def test_inverted_thresholds_raise(self):
        with pytest.raises(ValueError, match="plddt_unobs_threshold"):
            get_alphafold_structure(
                "P00533",
                ecd_ranges=[(25, 645)],
                plddt_warn_threshold=40.0,
                plddt_unobs_threshold=60.0,
            )


class TestGetAlphaFoldStructureWithMockedNetwork:
    def test_full_pipeline_egfr_mock(self, tmp_path):
        cif_lines = [
            "data_AF-P00533-F1",
            "_entry.id AF-P00533-F1",
            "loop_",
            "_atom_site.group_PDB",
            "_atom_site.id",
            "_atom_site.type_symbol",
            "_atom_site.label_atom_id",
            "_atom_site.label_alt_id",
            "_atom_site.label_comp_id",
            "_atom_site.label_asym_id",
            "_atom_site.label_entity_id",
            "_atom_site.label_seq_id",
            "_atom_site.pdbx_PDB_ins_code",
            "_atom_site.Cartn_x",
            "_atom_site.Cartn_y",
            "_atom_site.Cartn_z",
            "_atom_site.occupancy",
            "_atom_site.B_iso_or_equiv",
            "_atom_site.auth_seq_id",
            "_atom_site.auth_comp_id",
            "_atom_site.auth_asym_id",
            "_atom_site.auth_atom_id",
            "_atom_site.pdbx_PDB_model_num",
        ]
        for i, (seq_id, res_name, plddt, chain) in enumerate(EGFR_RESIDUES, 1):
            cif_lines.append(
                f"ATOM {i} C CA . {res_name} {chain} 1 {i} ? "
                f"{i}.0 {i*2}.0 {i*3}.0 1.00 {plddt:.2f} "
                f"{seq_id} {res_name} {chain} CA 1"
            )
        cif_lines.append("#")

        cache = tmp_path / "alphafold"
        cache.mkdir()
        cif_path = cache / "AF-P00533-F1-model_v4.cif"
        cif_path.write_text("\n".join(cif_lines))

        fake_meta = {"latestVersion": 4, "cifUrl": str(cif_path)}

        with patch("fetch_alphafold._fetch_metadata", return_value=fake_meta):
            result = get_alphafold_structure(
                "P00533",
                ecd_ranges=EGFR_ECD_RANGES,
                cache_dir=cache,
            )

        assert result.alphafold_fallback is True
        assert result.accession == "P00533"
        assert result.chain_id == "A"

        avoid_ids = {e.auth_seq_id for e in result.avoid_contributions}
        assert 620 in avoid_ids
        assert 635 in avoid_ids
        assert 645 in avoid_ids
        assert 25  not in avoid_ids
        assert 384 not in avoid_ids
        assert 408 not in avoid_ids

        for entry in result.avoid_contributions:
            assert entry.auth_seq_id == entry.uniprot_residue

        assert (620, 620) in result.low_plddt_ranges

        unob_ids = {s for s, e in result.unobs_ranges for _ in range(s, e + 1)}
        assert 635 in unob_ids
        assert 645 in unob_ids

        import json
        json.dumps(result.summary_fields())


# ═════════════════════════════════════════════════════════════════════════════
# Integration tests (real EBI AlphaFold API)
# ═════════════════════════════════════════════════════════════════════════════

@pytest.mark.network
class TestAlphaFoldIntegration:
    """
    Live tests against the EBI AlphaFold API.  Run with:

        pytest -v -m network tests/test_fetch_alphafold.py
    """

    @pytest.fixture(scope="class")
    def egfr_result(self, tmp_path_factory):
        cache = tmp_path_factory.mktemp("af_cache")
        return get_alphafold_structure("P00533", ecd_ranges=[(25, 645)], cache_dir=cache)

    def test_structure_parses(self, egfr_result):
        assert isinstance(egfr_result.structure, gemmi.Structure)

    def test_chain_a_present(self, egfr_result):
        assert "A" in [c.name for c in egfr_result.structure[0]]

    def test_version_is_integer(self, egfr_result):
        assert isinstance(egfr_result.version, int)
        assert egfr_result.version >= 4

    def test_avoid_contributions_is_list(self, egfr_result):
        assert isinstance(egfr_result.avoid_contributions, list)

    def test_all_avoid_entries_in_ecd(self, egfr_result):
        ecd_set = _ecd_set_from_ranges([(25, 645)])
        for entry in egfr_result.avoid_contributions:
            assert entry.auth_seq_id in ecd_set

    def test_auth_equals_uniprot_for_all_avoid_entries(self, egfr_result):
        """KEYSTONE — AlphaFold numbering must match UniProt directly."""
        for entry in egfr_result.avoid_contributions:
            assert entry.auth_seq_id == entry.uniprot_residue

    def test_cif_cached_on_second_call(self, tmp_path_factory):
        cache = tmp_path_factory.mktemp("af_cache2")
        get_alphafold_structure("P00533", [(25, 645)], cache_dir=cache)
        cif_files = list(cache.glob("AF-P00533-*.cif"))
        assert len(cif_files) == 1

    @pytest.mark.parametrize("accession, ecd_ranges", [
        ("P00533", [(25, 645)]),    # EGFR — type I
        ("P28907", [(43, 300)]),    # CD38 — type II
    ])
    def test_parametrized_proteins(self, accession, ecd_ranges, tmp_path):
        result = get_alphafold_structure(accession, ecd_ranges=ecd_ranges, cache_dir=tmp_path)
        assert result.alphafold_fallback is True
        assert result.ecd_mean_plddt > 0
        ecd_set = _ecd_set_from_ranges(ecd_ranges)
        for entry in result.avoid_contributions:
            assert entry.auth_seq_id in ecd_set
