"""Shared fixtures. Pipeline runs are network-backed but cached on disk, so the
suite is fast on re-run. Each fixture's pipeline result is memoised per session.
"""
from __future__ import annotations

import os
import pathlib

import pytest

# Keep every cache/download inside the project sandbox.
_ROOT = pathlib.Path(__file__).resolve().parents[1]
os.environ.setdefault("TRIMPROT_CACHE", str(_ROOT / "cache"))

from trimprot.pipeline import run_pipeline  # noqa: E402

# `pdb`  = curated soft hint (a known good structure that should appear in the list).
# `pick` = the exact structure the priority ladder must choose (HARD regression).
FIXTURES = [
    dict(id="EGFR", acc="P00533", pdb="1yy9", pick="9z9e", topo="type I", mp="C", prefer=True),
    dict(id="CD38", acc="P28907", pdb="7duo", pick="4cmh", topo="type II", mp="N", prefer=True),
    dict(id="CD44", acc="P16070", pdb="4pz3", pick="1poz", topo="type I", mp="C", prefer=False),
    dict(id="CD209", acc="Q9NNX6", pdb="1k9i", pick="2xr6", topo="type II", mp="N", prefer=False),
    dict(id="CD70", acc="P32970", pdb="7kx0", pick="7kx0", topo="type II", mp="N", prefer=False),
]

_RESULTS: dict[str, object] = {}


def run_for(meta: dict):
    acc = meta["acc"]
    if acc not in _RESULTS:
        _RESULTS[acc] = run_pipeline(None, accession=acc, prefer_antibody=meta["prefer"])
    return _RESULTS[acc]


@pytest.fixture(params=FIXTURES, ids=[f["id"] for f in FIXTURES])
def fixture_meta(request) -> dict:
    return request.param


@pytest.fixture
def prepared(fixture_meta):
    return fixture_meta, run_for(fixture_meta)


@pytest.fixture
def egfr():
    """The EGFR pipeline result (antibody-bound), for epitope-specific tests."""
    return run_for(next(f for f in FIXTURES if f["id"] == "EGFR"))
