"""Offline unit tests for the strict structure-selection priority ladder.

These build synthetic candidates (metrics set directly) and exercise
select_by_ladder with no network, proving the ordering of the rungs:
gate -> partner tier -> coverage -> completeness -> (experimental, resolution).
"""
from __future__ import annotations

from bindscout.structures import Candidate, select_by_ladder

ECD_TOTAL = 100  # coverage tolerance = 4/100 = 0.04


def mk(pdb, *, tier=0, cov=1.0, compl=1.0, res=2.0, predicted=False,
       label="X-ray", chain="A"):
    c = Candidate(pdb_id=pdb, chain_id=chain, unp_start=1, unp_end=100,
                  coverage=cov, resolution=res, method=label, tax_id=9606)
    c.ecd_coverage = cov
    c.completeness = compl
    c.partner_tier = tier
    c.predicted = predicted
    c.method_label = label
    return c


def pick(cands, **kw):
    chosen, survivors, _ = select_by_ladder(cands, ECD_TOTAL, **kw)
    return chosen, survivors


# ---- gate -------------------------------------------------------------------
def test_gate_filters_low_coverage_fragment():
    frag = mk("frag", cov=0.16, res=1.2)   # high-res fragment
    full = mk("full", cov=1.0, res=2.4)
    chosen, surv = pick([frag, full], min_ecd_coverage=0.40)
    assert chosen.pdb_id == "full"
    assert "frag" not in {c.pdb_id for c in surv}  # never reaches resolution compare


def test_gate_relaxes_if_it_would_empty_pool():
    a = mk("aaaa", cov=0.25, res=2.0)
    b = mk("bbbb", cov=0.24, res=1.5)
    chosen, surv = pick([a, b], min_ecd_coverage=0.40)
    assert {c.pdb_id for c in surv} >= {"aaaa"}  # not empty; gate relaxed


def test_usability_gate_drops_pathological_resolution():
    blob = mk("blob", tier=1, cov=1.0, res=25.0)   # 25Å partner complex
    apo = mk("apo_", tier=0, cov=1.0, res=1.4)
    chosen, _ = pick([blob, apo], max_resolution=9.0)
    assert chosen.pdb_id == "apo_"  # blob filtered, apo wins despite lower tier


# ---- subordination: biological criteria dominate quality --------------------
def test_partner_tier_beats_resolution():
    ab_lowres = mk("abab", tier=2, res=3.5)
    apo_hires = mk("apo_", tier=0, res=1.0)
    chosen, _ = pick([ab_lowres, apo_hires])
    assert chosen.pdb_id == "abab"  # better partner wins, resolution does NOT override


def test_coverage_beats_resolution_within_tier():
    hicov = mk("hicv", cov=1.0, res=3.0)
    locov = mk("locv", cov=0.5, res=1.0)   # 0.5 gap >> 0.04 tol
    chosen, _ = pick([hicov, locov])
    assert chosen.pdb_id == "hicv"


def test_completeness_breaks_coverage_tie_over_resolution():
    # CD44-shaped: equal coverage, NMR fully observed beats incomplete X-ray.
    nmr = mk("nmrr", cov=1.0, compl=1.0, res=None, label="NMR")
    xray = mk("xray", cov=1.0, compl=0.90, res=1.5, label="X-ray")  # 0.10 > 0.04 tol
    chosen, _ = pick([nmr, xray])
    assert chosen.pdb_id == "nmrr"


# ---- final tiebreaker: resolution only separates true ties ------------------
def test_resolution_breaks_tie_when_coverage_within_tolerance():
    a = mk("aaaa", cov=1.00, res=2.5)
    b = mk("bbbb", cov=0.98, res=1.5)   # 0.02 < 0.04 tol -> tied coverage
    chosen, _ = pick([a, b])
    assert chosen.pdb_id == "bbbb"      # higher resolution wins the tiebreak


def test_xray_and_em_tie_on_method_resolution_decides():
    em_better = mk("emem", res=1.8, label="cryo-EM")
    xray_worse = mk("xray", res=2.4, label="X-ray")
    chosen, _ = pick([em_better, xray_worse])
    assert chosen.pdb_id == "emem"      # method TYPE doesn't matter, resolution does


def test_nmr_loses_resolution_tiebreak_to_xray_when_otherwise_tied():
    nmr = mk("nmrr", res=None, label="NMR")
    xray = mk("xray", res=3.0, label="X-ray")
    chosen, _ = pick([nmr, xray])       # equal cov+compl -> step 4
    assert chosen.pdb_id == "xray"      # NMR has no resolution -> sorts worst


def test_predicted_never_beats_experimental_survivor():
    pred = mk("pred", cov=1.0, compl=1.0, res=None, predicted=True, label="predicted")
    exp = mk("xray", cov=1.0, compl=1.0, res=3.5, label="X-ray")
    chosen, _ = pick([pred, exp])
    assert chosen.pdb_id == "xray"      # experimental wins even at worse resolution


def test_predicted_chosen_only_as_sole_survivor():
    pred = mk("pred", cov=1.0, predicted=True, label="predicted", res=None)
    chosen, _ = pick([pred])
    assert chosen.pdb_id == "pred"
