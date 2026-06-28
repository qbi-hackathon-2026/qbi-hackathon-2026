"""Topology logic: extracellular ranges, type I/II, membrane-proximal terminus.

Everything here reads UniProt TOPO_DOM / TRANSMEM annotations. Topology is NEVER
inferred from raw residue position; the type I vs type II determination comes
from the ordering of the annotated extracellular domain relative to the annotated
transmembrane segment, which is itself an annotation-derived fact.

All ranges are 1-based UniProt numbering.
"""
from __future__ import annotations

from dataclasses import dataclass

from .uniprot import Feature, UniProtRecord


@dataclass(frozen=True)
class Range:
    start: int
    end: int

    def __contains__(self, pos: int) -> bool:
        return self.start <= pos <= self.end

    def __iter__(self):
        return iter(range(self.start, self.end + 1))


@dataclass
class MembraneProximal:
    terminus: str            # 'N' or 'C' end of the ECD that abuts the TM
    buffer: Range            # UniProt range excluded from hotspots
    ecd: Range               # the ECD range this applies to
    transmem: Range          # the adjacent TM segment
    topo_type: str           # 'type I' or 'type II'
    agrees: bool             # adjacency-derived terminus matches topo type
    note: str = ""


def get_extracellular_ranges(features: list[Feature]) -> list[Range]:
    """UniProt TOPO_DOM ranges flagged Extracellular, sorted by start."""
    out = [Range(f.start, f.end) for f in features
           if f.kind == "TOPO_DOM" and "extracellular" in f.description.lower()]
    out.sort(key=lambda r: r.start)
    return out


def get_transmem_ranges(features: list[Feature]) -> list[Range]:
    out = [Range(f.start, f.end) for f in features if f.kind == "TRANSMEM"]
    out.sort(key=lambda r: r.start)
    return out


def infer_topology_type(ecd: Range, transmem: Range) -> str:
    """type I = ECD N-terminal (before TM); type II = ECD C-terminal (after TM).

    Derived from the annotated ordering of the two segments, not raw position
    semantics.
    """
    if ecd.end <= transmem.start:
        return "type I"
    if ecd.start >= transmem.end:
        return "type II"
    # Overlap shouldn't happen for single-pass; fall back to midpoint compare.
    return "type I" if (ecd.start + ecd.end) < (transmem.start + transmem.end) else "type II"


def get_membrane_proximal_terminus(ecd_ranges: list[Range],
                                   transmem: list[Range],
                                   buffer: int = 12) -> MembraneProximal | None:
    """Find which ECD terminus is adjacent to a TM segment and the buffer to drop.

    Strategy: among all (ECD range, TM segment) pairs, pick the pair with the
    smallest sequence gap between an ECD terminus and a TM terminus. The ECD
    terminus on the small-gap side is membrane-proximal. The result is then
    cross-checked against the type I/II topology (type I -> C terminus,
    type II -> N terminus) and `agrees` is set accordingly.
    """
    if not ecd_ranges or not transmem:
        return None

    best = None  # (gap, terminus, ecd, tm)
    for ecd in ecd_ranges:
        for tm in transmem:
            # Gap from ECD's N-terminal end (start) to a TM segment on its left.
            gap_n = ecd.start - tm.end          # TM ends just before ECD starts
            # Gap from ECD's C-terminal end to a TM segment on its right.
            gap_c = tm.start - ecd.end          # TM starts just after ECD ends
            for terminus, gap in (("N", gap_n), ("C", gap_c)):
                if gap < 0:
                    continue  # TM is on the wrong side for this terminus
                if best is None or gap < best[0]:
                    best = (gap, terminus, ecd, tm)

    if best is None:
        return None

    _gap, terminus, ecd, tm = best
    topo_type = infer_topology_type(ecd, tm)
    expected = "C" if topo_type == "type I" else "N"
    agrees = terminus == expected

    if terminus == "N":
        b_start = ecd.start
        b_end = min(ecd.end, ecd.start + buffer - 1)
    else:  # 'C'
        b_end = ecd.end
        b_start = max(ecd.start, ecd.end - buffer + 1)

    note = (
        f"membrane-proximal {terminus}-terminal end of ECD {ecd.start}-{ecd.end} "
        f"abuts TM {tm.start}-{tm.end} (gap {_gap}); {topo_type}; "
        f"{'consistent' if agrees else 'INCONSISTENT'} with topology. "
        f"Sequence-proximity proxy (no explicit membrane model)."
    )
    return MembraneProximal(
        terminus=terminus,
        buffer=Range(b_start, b_end),
        ecd=ecd,
        transmem=tm,
        topo_type=topo_type,
        agrees=agrees,
        note=note,
    )


def summarize_topology(rec: UniProtRecord, buffer: int = 12) -> dict:
    """Convenience bundle used by the pipeline/summary."""
    ecd = get_extracellular_ranges(rec.features)
    tm = get_transmem_ranges(rec.features)
    mp = get_membrane_proximal_terminus(ecd, tm, buffer=buffer)
    return {
        "ecd_ranges": [(r.start, r.end) for r in ecd],
        "transmem_ranges": [(r.start, r.end) for r in tm],
        "membrane_proximal": mp,
    }
