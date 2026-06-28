"""Build the avoid set (PDB-auth numbered) with a reason per residue.

Union of: glycosylation sites, residues near structural glycan HETATM, disulfide
cysteines, modified residues, unobserved UniProt residues, transmembrane residues,
and the membrane-proximal buffer. Everything is mapped through SIFTS into PDB
author numbering; the glycan-proximity term is read directly off coordinates.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional

import gemmi

from .sifts import SiftsMapping
from .structio import sugar_atoms
from .topology import Range


@dataclass
class AvoidEntry:
    chain: str
    num: int
    icode: str
    uniprot: Optional[int]
    reasons: list[str] = field(default_factory=list)

    @property
    def key(self) -> tuple[str, int, str]:
        return (self.chain, self.num, self.icode)


class AvoidSet:
    def __init__(self):
        self._entries: dict[tuple[str, int, str], AvoidEntry] = {}

    def add(self, chain: str, num: int, icode: str, uniprot: Optional[int],
            reason: str) -> None:
        key = (chain, num, icode.strip())
        e = self._entries.get(key)
        if e is None:
            e = AvoidEntry(chain, num, icode.strip(), uniprot, [])
            self._entries[key] = e
        if uniprot is not None and e.uniprot is None:
            e.uniprot = uniprot
        if reason not in e.reasons:
            e.reasons.append(reason)

    def contains(self, chain: str, num: int, icode: str = "") -> bool:
        return (chain, num, icode.strip()) in self._entries

    def get(self, chain: str, num: int, icode: str = "") -> Optional[AvoidEntry]:
        return self._entries.get((chain, num, icode.strip()))

    def entries(self) -> list[AvoidEntry]:
        return sorted(self._entries.values(), key=lambda e: (e.chain, e.num, e.icode))

    def __len__(self) -> int:
        return len(self._entries)


def build_avoid_set(mapping: SiftsMapping, structure: gemmi.Structure,
                    target_chain: str, *,
                    glyco_sites: set[int],
                    disulfide_cys: list[int],
                    mod_res: list[int],
                    tm_ranges: list[Range],
                    membrane_buffer: Optional[Range],
                    ecd_ranges: list[Range],
                    glycan_cutoff: float = 5.0) -> AvoidSet:
    avoid = AvoidSet()
    u2p = mapping.unp_to_pdb(target_chain)
    p2u = mapping.pdb_to_unp(target_chain)

    def add_unp(unp: int, reason: str):
        hit = u2p.get(unp)
        if hit is not None:
            num, icode = hit
            avoid.add(target_chain, num, icode, unp, reason)

    for site in sorted(glyco_sites):
        add_unp(site, "glycosylation")
    for cys in sorted(set(disulfide_cys)):
        add_unp(cys, "disulfide")
    for mr in sorted(set(mod_res)):
        add_unp(mr, "modified-residue")
    for r in tm_ranges:
        for unp in r:
            add_unp(unp, "transmembrane")
    if membrane_buffer is not None:
        for unp in membrane_buffer:
            add_unp(unp, "membrane-proximal")

    # Unobserved UniProt residues that fall inside the ECD.
    def in_ecd(unp: int) -> bool:
        return any(unp in r for r in ecd_ranges) if ecd_ranges else True

    for rm in mapping.for_chain(target_chain):
        if not rm.observed and in_ecd(rm.unp):
            avoid.add(target_chain, rm.pdb_num, rm.icode, rm.unp, "unobserved")

    # Residues with any heavy atom within cutoff of a structural glycan.
    _add_glycan_proximal(avoid, structure, target_chain, p2u, glycan_cutoff)
    return avoid


def _add_glycan_proximal(avoid: AvoidSet, structure: gemmi.Structure,
                         target_chain: str, p2u: dict, cutoff: float) -> None:
    positions, sugar_keys = sugar_atoms(structure)
    if not positions:
        return
    model = structure[0]
    ns = gemmi.NeighborSearch(model, structure.cell, cutoff).populate()
    flagged: set[tuple[int, str]] = set()
    for pos in positions:
        for mark in ns.find_atoms(pos, "\0", radius=cutoff):
            hit = mark.to_cra(model)
            if hit.chain.name != target_chain:
                continue
            # don't flag sugar-on-sugar
            if (hit.chain.name, hit.residue.seqid.num) in sugar_keys:
                continue
            num = hit.residue.seqid.num
            icode = hit.residue.seqid.icode.strip()
            if (num, icode) in flagged:
                continue
            flagged.add((num, icode))
            avoid.add(target_chain, num, icode, p2u.get((num, icode)),
                      "glycan-proximal")
