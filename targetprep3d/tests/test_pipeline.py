"""Network-backed pipeline tests over the parametrized target fixtures.

Run offline subset:  uv run pytest -m "not network"
Run everything:      uv run pytest
"""
from __future__ import annotations

import json
import warnings

import pytest

import gemmi

from targetprep3d.hotspots import _cb_positions
from targetprep3d.interface import is_amino_acid, partner_contact_counts
from targetprep3d.outputs import emit_outputs, format_hotspot_string
from targetprep3d.sifts import get_sifts_mapping
from targetprep3d.structio import load_structure
from targetprep3d.topology import get_extracellular_ranges, get_transmem_ranges

pytestmark = pytest.mark.network

PATCH_RADIUS = 11.0
PATCH_SIZE = 8

# EGFR ECD domain spans in PDB author (mature) numbering.
EGFR_DOMAINS = {"I": (1, 165), "II": (166, 309), "III": (310, 480), "IV": (481, 613)}


def _domain_of(num: int):
    for name, (lo, hi) in EGFR_DOMAINS.items():
        if lo <= num <= hi:
            return name
    return None


# ---- generic per-fixture assertions ----------------------------------------
def test_sifts_mapping_non_empty(prepared):
    meta, r = prepared
    chain = r.chains.target_chain
    assert len(r.mapping.for_chain(chain)) > 0, "SIFTS mapping empty"


def test_trimmed_numbering_equals_original(prepared):
    meta, r = prepared
    display = r.loaded.display[0]
    inv = {new: old for old, new in r.trim.rename.items()}
    for ch in r.trim.structure[0]:
        orig_name = inv.get(ch.name, ch.name)
        orig = next((c for c in display if c.name == orig_name), None)
        assert orig is not None
        orig_keys = {(res.seqid.num, res.seqid.icode) for res in orig}
        for res in ch:
            assert (res.seqid.num, res.seqid.icode) in orig_keys, (
                f"{meta['id']}: trimmed residue {res.seqid.num} not in original")


def test_no_hotspot_in_avoid_set(prepared):
    meta, r = prepared
    for h in r.hotspots:
        # hotspots are relabelled to the trimmed chain; avoid is keyed on the AU
        # target chain, but residue (num, icode) identity is what matters.
        assert not r.avoid.contains(r.chains.target_chain, h.num, h.icode), (
            f"{meta['id']}: hotspot {h.num} is in the avoid set")


def test_no_hotspot_in_membrane_buffer(prepared):
    meta, r = prepared
    mb = r.membrane_buffer_auth
    if mb is None:
        pytest.skip("no membrane buffer mapped")
    for h in r.hotspots:
        assert not (mb.start <= h.num <= mb.end), (
            f"{meta['id']}: hotspot {h.num} in membrane-proximal buffer "
            f"{mb.start}-{mb.end}")


def test_bindcraft_config_and_hotspot_string(prepared, tmp_path):
    meta, r = prepared
    out = emit_outputs(r, tmp_path)
    cfg = json.loads((out / "bindcraft_config.json").read_text())
    # exact schema keys
    assert set(cfg) == {"design_path", "binder_name", "starting_pdb", "chains",
                        "target_hotspot_residues", "lengths",
                        "number_of_final_designs"}
    assert cfg["lengths"] == [65, 150]
    assert cfg["number_of_final_designs"] == 100
    # bindcraft carries only the focused patch; hotspots.csv carries the full list
    patch = format_hotspot_string(r.patch)
    if not patch:
        # empty is allowed only with a logged reason
        assert any("no hotspots" in w for w in r.warnings)
    else:
        assert cfg["target_hotspot_residues"] == patch
        assert len(r.patch) <= len(r.hotspots)
        # full list is preserved in hotspots.csv
        import csv as _csv
        rows = list(_csv.DictReader((out / "hotspots.csv").open()))
        assert len(rows) == len(r.hotspots)


def test_expected_pdb_in_candidate_list_soft(prepared):
    meta, r = prepared
    if not meta["pdb"]:
        pytest.skip("no expected pdb for this fixture")
    ids = {c.pdb_id for c in r.choice.candidates}
    if meta["pdb"] not in ids:
        warnings.warn(f"{meta['id']}: expected pdb {meta['pdb']} not in candidates")


# ---- HARD regression: the strict ladder must reproduce these exact picks ----
def test_ladder_reproduces_expected_pick(prepared):
    """If any pick flips, the ladder is being applied as a weighted sum somewhere —
    fix the ladder, do NOT update the expected pick."""
    meta, r = prepared
    assert r.choice.chosen.pdb_id == meta["pick"], (
        f"{meta['id']}: ladder chose {r.choice.chosen.pdb_id}, expected {meta['pick']}")


def test_egfr_high_res_low_coverage_fragment_is_gated_out():
    """EGFR 3p0y (1.8Å but only domain III, ECD coverage ~0.33) must be filtered
    at the coverage gate so resolution never rescues a fragment."""
    from targetprep3d.structures import (ecd_coverage, search_structures)
    from targetprep3d.uniprot import resolve_uniprot
    from targetprep3d.topology import get_extracellular_ranges
    rec = resolve_uniprot(accession="P00533")
    ecd = get_extracellular_ranges(rec.features)
    choice = search_structures("P00533", ecd, prefer_antibody=True)
    frag = next((c for c in choice.candidates if c.pdb_id == "3p0y"), None)
    if frag is None:
        pytest.skip("3p0y not in candidate list")
    assert frag.ecd_coverage < 0.40            # below the gate
    assert choice.chosen.pdb_id != "3p0y"      # never selected
    assert choice.chosen.ecd_coverage >= 0.40  # the pick clears the gate


def test_topology_terminus_matches_type(prepared):
    meta, r = prepared
    mp = r.membrane_proximal
    assert mp is not None
    assert mp.topo_type == meta["topo"], f"{meta['id']} topology mismatch"
    assert mp.terminus == meta["mp"], f"{meta['id']} terminus mismatch"


def test_extracellular_location_relative_to_tm(prepared):
    meta, r = prepared
    ecd = get_extracellular_ranges(r.record.features)
    tm = get_transmem_ranges(r.record.features)
    assert ecd and tm
    tm_start = min(t.start for t in tm)
    tm_end = max(t.end for t in tm)
    if meta["topo"] == "type I":
        # ECD is N-terminal: lies before the TM
        assert min(e.start for e in ecd) < tm_start
        assert max(e.end for e in ecd) <= tm_end or min(e.start for e in ecd) < tm_start
    else:
        # type II ECD is C-terminal: lies after the TM
        assert max(e.end for e in ecd) > tm_end


# ---- EGFR keystone: numbering layer canary ---------------------------------
def test_egfr_keystone_nonzero_consistent_offset():
    """1YY9 UniProt->PDB(auth) offset must be non-zero and consistent.

    This is the canary that the numbering layer never assumes 1:1 identity. The
    PDBe author_residue_number is null for 1YY9, so this also exercises the
    label-seq + coordinate fallback path.
    """
    ld = load_structure("1yy9", assembly="protomer")
    m = get_sifts_mapping("1yy9", "P00533", structure=ld.analysis, target_chain="A")
    assert len(m.for_chain("A")) > 0
    off = m.consistent_offset("A")
    assert off is not None, "offset not consistent across the mapped segment"
    assert off != 0, "1:1 identity mapping assumed — numbering layer is broken"
    assert off == -24  # mature (signal peptide 1-24 cleaved)


# ---- FIX 1: partner-chain purity (hard regression) -------------------------
def test_partner_contacts_are_amino_acids_only(egfr):
    """Every counted partner atom must be in a designated antibody chain AND a
    standard amino acid — no glycan/ligand/ion/water may contribute a contact."""
    r = egfr
    assert not (r.choice.apo_fallback or not r.interface_chains)
    struct = r.loaded.analysis
    tchain = r.chains.target_chain
    ichains = set(r.interface_chains)
    cutoff = 5.0

    # Independent neighbour scan (does NOT reuse production counting).
    model = struct[0]
    ns = gemmi.NeighborSearch(model, struct.cell, cutoff).populate()
    aa_only: dict[tuple[int, str], int] = {}
    nonaa_seen = 0
    for chain in model:
        if chain.name != tchain:
            continue
        for res in chain:
            if not is_amino_acid(res):
                continue
            if not any(res.seqid.num in rng for rng in r.ecd_auth_ranges):
                continue
            aa = 0
            for atom in res:
                if atom.is_hydrogen():
                    continue
                for mk in ns.find_atoms(atom.pos, "\0", radius=cutoff):
                    cra = mk.to_cra(model)
                    if cra.chain.name not in ichains or cra.atom.is_hydrogen():
                        continue
                    if is_amino_acid(cra.residue):
                        aa += 1
                    else:
                        nonaa_seen += 1
            if aa:
                aa_only[(res.seqid.num, res.seqid.icode.strip())] = aa

    # Production interface must equal the independent amino-acid-only scan.
    prod = {(ir.num, ir.icode): ir.contact_count for ir in r.interface}
    assert prod == aa_only, "interface counts include non-amino-acid contacts"

    # Regression is meaningful: non-AA partner atoms DO sit near the target
    # (waters/ions/PEG share the antibody chain ids in 9Z9E) yet none are counted.
    assert nonaa_seen > 0, "no non-AA partner atoms present — test not exercising the filter"


def test_partner_counts_drop_when_excluding_non_amino_acids(egfr):
    r = egfr
    args = (r.loaded.analysis, r.chains.target_chain, r.interface_chains,
            5.0, r.ecd_auth_ranges)
    aa = partner_contact_counts(*args, amino_acid_only=True)
    raw = partner_contact_counts(*args, amino_acid_only=False)
    assert sum(raw.values()) > sum(aa.values())  # non-AA contacts were really removed


# ---- FIX 2: patch coherence / shape ----------------------------------------
def test_egfr_patch_is_single_domain(egfr):
    r = egfr
    doms = {_domain_of(h.num) for h in r.patch}
    assert None not in doms, "patch residue outside defined EGFR domains"
    assert len(doms) == 1, f"patch spans multiple domains: {doms}"


def test_patch_shape_and_contiguity(prepared):
    meta, r = prepared
    if not r.patch:
        pytest.skip("empty patch (logged reason)")
    assert len(r.patch) <= PATCH_SIZE
    pos = _cb_positions(r.loaded.analysis, r.chains.target_chain)
    keys = [(h.num, h.icode.strip()) for h in r.patch]
    if len(keys) == 1:
        return
    for k in keys:
        assert any(k != j and pos[k].dist(pos[j]) <= PATCH_RADIUS
                   for j in keys), f"{meta['id']}: patch residue {k} not contiguous"


def test_patch_is_subset_and_matches_bindcraft(prepared, tmp_path):
    meta, r = prepared
    out = emit_outputs(r, tmp_path)
    import csv as _csv
    rows = list(_csv.DictReader((out / "hotspots.csv").open()))
    assert "in_patch" in rows[0]
    full = {row["pdb_auth_residue"] for row in rows}
    patch = {f"{h.num}{h.icode}" if h.icode else str(h.num) for h in r.patch}
    assert patch <= full, "patch is not a subset of the full hotspot list"
    n_in_patch = sum(1 for row in rows if row["in_patch"] == "True")
    assert n_in_patch == len(r.patch)
    cfg = json.loads((out / "bindcraft_config.json").read_text())
    assert cfg["target_hotspot_residues"] == format_hotspot_string(r.patch)


def test_egfr_known_hub_in_patch_soft(egfr):
    patch_nums = {h.num for h in egfr.patch}
    known = {409, 408, 384}
    if not (patch_nums & known):
        warnings.warn(f"EGFR patch {sorted(patch_nums)} misses known hub 409/408/384")
