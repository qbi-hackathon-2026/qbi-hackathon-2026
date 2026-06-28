"""Offline unit tests for the pure deterministic logic (no network)."""
from __future__ import annotations

import gemmi

from targetprep3d.glyco import predict_glycosylation, scan_nglyc_sequons
from targetprep3d.hotspots import Hotspot
from targetprep3d.outputs import format_hotspot_string
from targetprep3d.topology import (
    Range,
    get_extracellular_ranges,
    get_membrane_proximal_terminus,
    infer_topology_type,
)
from targetprep3d.trim import trim_structure
from targetprep3d.uniprot import Feature


# ---- topology / membrane-proximal ------------------------------------------
def test_type_ii_membrane_proximal_is_n_terminal():
    # CD38-like: TM 22-42, ECD 43-300
    mp = get_membrane_proximal_terminus([Range(43, 300)], [Range(22, 42)], buffer=12)
    assert mp.terminus == "N"
    assert mp.topo_type == "type II"
    assert mp.agrees
    assert (mp.buffer.start, mp.buffer.end) == (43, 54)


def test_type_i_membrane_proximal_is_c_terminal():
    # EGFR-like: ECD 25-645, TM 646-668
    mp = get_membrane_proximal_terminus([Range(25, 645)], [Range(646, 668)], buffer=12)
    assert mp.terminus == "C"
    assert mp.topo_type == "type I"
    assert mp.agrees
    assert (mp.buffer.start, mp.buffer.end) == (634, 645)


def test_infer_topology_type():
    assert infer_topology_type(Range(25, 645), Range(646, 668)) == "type I"
    assert infer_topology_type(Range(43, 300), Range(22, 42)) == "type II"


def test_membrane_buffer_respects_buffer_size():
    mp = get_membrane_proximal_terminus([Range(43, 300)], [Range(22, 42)], buffer=20)
    assert (mp.buffer.start, mp.buffer.end) == (43, 62)


def test_extracellular_ranges_from_features():
    feats = [
        Feature("TOPO_DOM", 1, 21, "Cytoplasmic"),
        Feature("TRANSMEM", 22, 42, "Helical"),
        Feature("TOPO_DOM", 43, 300, "Extracellular"),
    ]
    ecd = get_extracellular_ranges(feats)
    assert [(r.start, r.end) for r in ecd] == [(43, 300)]


# ---- glycosylation ----------------------------------------------------------
def test_sequon_scan_excludes_proline_and_requires_st():
    # N-I-T (yes), N-P-T (no, X=P), N-G-S (yes), N-A-A (no)
    assert scan_nglyc_sequons("ANITxNPTxNGSxNAA") == {2, 10}


def test_predict_glycosylation_unions_annotated_sites():
    seq = "ANITAA"
    carbo = [Feature("CARBOHYD", 5, 5, "O-linked")]
    sites = predict_glycosylation(seq, carbo)
    assert 2 in sites and 5 in sites


# ---- hotspot string ---------------------------------------------------------
def test_hotspot_string_compression():
    hs = [Hotspot("A", n, "", None, 1, "i")
          for n in (180, 182, 190, 191, 192, 193)]
    assert format_hotspot_string(hs) == "A180,A182,A190-A193"


def test_hotspot_string_handles_insertion_codes():
    hs = [Hotspot("A", 100, "", None, 1, "i"), Hotspot("A", 100, "A", None, 1, "i")]
    assert format_hotspot_string(hs) == "A100,A100A"


# ---- trim preserves author numbering + renames chains -----------------------
def _toy_structure() -> gemmi.Structure:
    st = gemmi.Structure()
    st.cell = gemmi.UnitCell(1, 1, 1, 90, 90, 90)
    st.spacegroup_hm = "P 1"
    model = gemmi.Model("1")
    # two chains 'H' (target) and 'L', residues 10..14 with one insertion code
    for cname in ("H", "L"):
        ch = gemmi.Chain(cname)
        for num, ic in [(10, " "), (11, " "), (12, "A"), (13, " "), (14, " ")]:
            r = gemmi.Residue()
            r.name = "ALA"
            r.seqid = gemmi.SeqId(num, ic)
            a = gemmi.Atom()
            a.name = "CA"
            a.element = gemmi.Element("C")
            a.pos = gemmi.Position(num * 1.0, 0, 0)
            r.add_atom(a)
            ch.add_residue(r)
        model.add_chain(ch)
    st.add_model(model)
    st.setup_entities()
    return st


def test_trim_preserves_numbering_and_renames_target_first():
    st = _toy_structure()
    res = trim_structure(st, keep_chains=["H", "L"], keep_ranges=[Range(11, 13)],
                         target_chain="H")
    # target 'H' becomes 'A'
    assert res.target_chain == "A"
    assert res.rename["H"] == "A"
    # numbering preserved (11, 12A, 13) and out-of-range 10/14 dropped
    out = res.structure[0]["A"]
    got = [(r.seqid.num, r.seqid.icode.strip()) for r in out]
    assert got == [(11, ""), (12, "A"), (13, "")]
