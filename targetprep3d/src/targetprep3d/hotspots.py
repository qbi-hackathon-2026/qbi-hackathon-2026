"""Final hotspot selection = interface candidates minus the avoid set.

Every removal is logged with the reasons that disqualified the residue (including
membrane-proximal removals), so the summary can explain exactly why a contacting
residue did not become a hotspot.
"""
from __future__ import annotations

from dataclasses import dataclass

import gemmi

from .avoid import AvoidSet
from .interface import InterfaceResidue


@dataclass
class Hotspot:
    chain: str
    num: int
    icode: str
    uniprot: int | None
    contact_count: int
    source: str


@dataclass
class Removal:
    chain: str
    num: int
    icode: str
    reasons: list[str]


def filter_hotspots(candidates: list[InterfaceResidue], avoid: AvoidSet,
                    chain: str, pdb_to_unp: dict[tuple[int, str], int]
                    ) -> tuple[list[Hotspot], list[Removal]]:
    hotspots: list[Hotspot] = []
    removals: list[Removal] = []
    for c in candidates:
        entry = avoid.get(chain, c.num, c.icode)
        if entry is not None:
            removals.append(Removal(chain, c.num, c.icode, list(entry.reasons)))
            continue
        hotspots.append(Hotspot(
            chain=chain, num=c.num, icode=c.icode,
            uniprot=pdb_to_unp.get((c.num, c.icode.strip())),
            contact_count=c.contact_count, source=c.source,
        ))
    return hotspots, removals


def _cb_positions(structure: gemmi.Structure, target_chain: str
                  ) -> dict[tuple[int, str], gemmi.Position]:
    """Per residue: Cβ position (Cα for glycine / missing Cβ)."""
    pos: dict[tuple[int, str], gemmi.Position] = {}
    model = structure[0]
    for chain in model:
        if chain.name != target_chain:
            continue
        for res in chain:
            atom = res.find_atom("CB", "*") or res.get_ca()
            if atom is not None:
                pos[(res.seqid.num, res.seqid.icode.strip())] = atom.pos
    return pos


def select_patch(hotspots: list[Hotspot], structure: gemmi.Structure,
                 target_chain: str, radius: float = 11.0, cap: int = 8
                 ) -> list[Hotspot]:
    """Collapse hotspots to one small, spatially contiguous epitope patch.

    Seed with the highest-scoring residue (contact count to the antibody; for apo
    targets this is the exposure/contact proxy). Then iteratively add the
    next-highest-scoring residue whose Cβ (Cα for glycine) lies within ``radius``
    of any current patch member, until ``cap`` residues or nothing else is close
    enough. The result is a single radial patch, not every interface residue.
    """
    if not hotspots:
        return []
    pos = _cb_positions(structure, target_chain)
    ranked = sorted(hotspots, key=lambda h: (-h.contact_count, h.num, h.icode))
    ranked = [h for h in ranked if (h.num, h.icode.strip()) in pos]
    if not ranked:
        return hotspots[:cap]

    hub = ranked[0]
    patch = [hub]
    patch_keys = {(hub.num, hub.icode.strip())}

    while len(patch) < cap:
        added = False
        for h in ranked:
            k = (h.num, h.icode.strip())
            if k in patch_keys:
                continue
            if any(pos[k].dist(pos[pk]) <= radius for pk in patch_keys):
                patch.append(h)
                patch_keys.add(k)
                added = True
                break  # restart so the next addition is again the top-scoring one
        if not added:
            break
    return patch

