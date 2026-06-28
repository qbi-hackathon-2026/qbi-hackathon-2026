"""Glycosylation prediction: annotated CARBOHYD sites + N-X-[S/T] sequon scan."""
from __future__ import annotations

from .uniprot import Feature


def scan_nglyc_sequons(seq: str) -> set[int]:
    """N-X-[S/T] with X != P. Returns 1-based positions of the Asn (the N)."""
    out: set[int] = set()
    n = len(seq)
    for i in range(n - 2):
        if seq[i] != "N":
            continue
        x = seq[i + 1]
        s = seq[i + 2]
        if x != "P" and s in ("S", "T"):
            out.add(i + 1)  # 1-based position of the N
    return out


def predict_glycosylation(seq: str, carbohyd: list[Feature]) -> set[int]:
    """Union of annotated CARBOHYD sites and predicted N-glycosylation sequons.

    Returns a set of 1-based UniProt positions. Annotated sites cover both N- and
    O-linked (and confirmed sites); the sequon scan adds predicted N-linked Asn.
    """
    sites: set[int] = {f.start for f in carbohyd if f.kind == "CARBOHYD"}
    sites |= scan_nglyc_sequons(seq)
    return sites
