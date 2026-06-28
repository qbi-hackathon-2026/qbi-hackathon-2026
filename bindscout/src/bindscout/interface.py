"""Interface / epitope detection on the target chain.

Primary path: target residues with any heavy atom within `cutoff` of a partner
heavy atom, ranked by contact count. Apo fallback: when there is no partner, flag
exposed secondary-structure residues (helix/strand from the deposited header,
loops skipped) so designed binders aim at a structured surface patch.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

import gemmi

from .structio import SUGAR_CCDS
from .topology import Range


@dataclass
class InterfaceResidue:
    num: int
    icode: str
    contact_count: int
    source: str            # 'interface' or 'apo_exposed'


def _in_ranges(num: int, ranges: Optional[list[Range]]) -> bool:
    if not ranges:
        return True
    return any(num in r for r in ranges)


def is_amino_acid(res: gemmi.Residue) -> bool:
    """Standard (or modified, e.g. MSE) amino-acid residue, never sugar/ion/water."""
    info = gemmi.find_tabulated_residue(res.name)
    return info is not None and info.is_amino_acid()


def _is_polymer_residue(res: gemmi.Residue) -> bool:
    return not res.is_water() and res.name not in SUGAR_CCDS and len(res) > 0


def detect_interface(structure: gemmi.Structure, target_chain: str,
                     partner_chains: list[str], cutoff: float = 5.0, *,
                     apo_fallback: bool = False,
                     ecd_auth_ranges: Optional[list[Range]] = None
                     ) -> list[InterfaceResidue]:
    if apo_fallback or not partner_chains:
        return _detect_apo_exposed(structure, target_chain, ecd_auth_ranges)
    counts = partner_contact_counts(structure, target_chain, partner_chains,
                                    cutoff, ecd_auth_ranges, amino_acid_only=True)
    out = [InterfaceResidue(num, icode, n, "interface")
           for (num, icode), n in counts.items() if n > 0]
    out.sort(key=lambda r: (-r.contact_count, r.num, r.icode))
    return out


def partner_contact_counts(structure: gemmi.Structure, target_chain: str,
                           partner_chains: list[str], cutoff: float,
                           ecd_auth_ranges: Optional[list[Range]],
                           amino_acid_only: bool = True
                           ) -> dict[tuple[int, str], int]:
    """Heavy-atom contacts from each target amino-acid residue to the partner set.

    A partner atom contributes ONLY when its chain is one of the designated
    partner/antibody chains AND (when amino_acid_only) its residue is a standard
    amino acid — so the target's own glycans, ligands, ions, and waters that share
    a partner chain id never manufacture an interface contact. amino_acid_only is
    exposed so tests can compare amino-acid-only vs raw counts.
    """
    model = structure[0]
    partner_set = set(partner_chains)
    ns = gemmi.NeighborSearch(model, structure.cell, cutoff).populate()

    counts: dict[tuple[int, str], int] = {}
    for chain in model:
        if chain.name != target_chain:
            continue
        for res in chain:
            if not is_amino_acid(res):
                continue
            if not _in_ranges(res.seqid.num, ecd_auth_ranges):
                continue
            count = 0
            for atom in res:
                if atom.is_hydrogen():
                    continue
                for mark in ns.find_atoms(atom.pos, "\0", radius=cutoff):
                    cra = mark.to_cra(model)
                    if cra.chain.name not in partner_set:
                        continue
                    if cra.atom.is_hydrogen():
                        continue
                    if amino_acid_only and not is_amino_acid(cra.residue):
                        continue
                    count += 1
            if count > 0:
                counts[(res.seqid.num, res.seqid.icode.strip())] = count
    return counts


def _ss_residue_keys(structure: gemmi.Structure, chain_name: str) -> set[tuple[int, str]]:
    """Residue (num, icode) keys that lie in a helix or strand per the header."""
    keys: set[tuple[int, str]] = set()

    def add_span(start: gemmi.AtomAddress, end: gemmi.AtomAddress):
        if start.chain_name != chain_name and end.chain_name != chain_name:
            return
        lo, hi = start.res_id.seqid.num, end.res_id.seqid.num
        for n in range(min(lo, hi), max(lo, hi) + 1):
            keys.add((n, ""))

    for h in structure.helices:
        add_span(h.start, h.end)
    for sheet in structure.sheets:
        for strand in sheet.strands:
            add_span(strand.start, strand.end)
    return keys


def _detect_apo_exposed(structure, target_chain, ecd_auth_ranges
                        ) -> list[InterfaceResidue]:
    """Exposed structured residues: in a helix/strand and solvent-exposed.

    Exposure is approximated by counting same-chain Cα neighbours within 10 Å;
    fewer neighbours -> more surface-exposed. Loops are skipped via the SS filter.
    """
    model = structure[0]
    ss_keys = _ss_residue_keys(structure, target_chain)

    cas: list[tuple[gemmi.Residue, gemmi.Atom]] = []
    for chain in model:
        if chain.name != target_chain:
            continue
        for res in chain:
            if not _is_polymer_residue(res):
                continue
            ca = res.get_ca()
            if ca is not None:
                cas.append((res, ca))

    ns = gemmi.NeighborSearch(model, structure.cell, 10.0).populate()
    scored: list[tuple[int, gemmi.Residue]] = []
    for res, ca in cas:
        key = (res.seqid.num, res.seqid.icode.strip())
        if key not in ss_keys:
            continue
        if not _in_ranges(res.seqid.num, ecd_auth_ranges):
            continue
        neighbours = 0
        for mark in ns.find_atoms(ca.pos, "\0", radius=10.0):
            cra = mark.to_cra(model)
            if cra.chain.name == target_chain and cra.atom.name == "CA":
                neighbours += 1
        scored.append((neighbours, res))

    if not scored:
        return []
    # Exposed = below-median neighbour count. Rank most-exposed first.
    counts = sorted(n for n, _ in scored)
    median = counts[len(counts) // 2]
    out = []
    for neighbours, res in scored:
        if neighbours <= median:
            out.append(InterfaceResidue(res.seqid.num, res.seqid.icode.strip(),
                                        0, "apo_exposed"))
    out.sort(key=lambda r: (r.num, r.icode))
    return out
