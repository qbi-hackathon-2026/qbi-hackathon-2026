"""Trim a structure to the target ECD, preserving author numbering + icodes.

Keeps the target protein chains' polymer residues whose author number falls in the
kept ECD ranges; drops other chains, ligands, glycans, and waters. Numbering and
insertion codes are carried over verbatim (never renumbered). Orphan fragments
(short, isolated kept stretches at domain boundaries) are warned about.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import gemmi

from .topology import Range


_CLEAN_NAMES = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"


@dataclass
class TrimResult:
    structure: gemmi.Structure
    kept_chains: list[str]            # cleaned, single-letter names
    kept_counts: dict[str, int]
    target_chain: str                 # cleaned name of the target chain
    rename: dict[str, str]            # original kept name -> cleaned name
    warnings: list[str] = field(default_factory=list)


def _polymer_keys(chain: gemmi.Chain) -> set[tuple[int, str]]:
    poly = chain.get_polymer()
    return {(r.seqid.num, r.seqid.icode) for r in poly}


def _in_ranges(num: int, ranges: list[Range]) -> bool:
    return any(num in r for r in ranges)


def _fragment_warnings(chain_name: str, nums: list[int]) -> list[str]:
    """Warn about isolated short fragments (likely orphans at trim boundaries)."""
    if not nums:
        return []
    nums = sorted(nums)
    fragments: list[list[int]] = [[nums[0]]]
    for n in nums[1:]:
        if n == fragments[-1][-1] + 1:
            fragments[-1].append(n)
        else:
            fragments.append([n])
    warns = []
    for frag in fragments:
        if len(frag) < 3 and len(fragments) > 1:
            warns.append(
                f"chain {chain_name}: orphan fragment {frag[0]}-{frag[-1]} "
                f"({len(frag)} residue(s)) isolated from the main domain")
    return warns


def trim_structure(structure: gemmi.Structure, keep_chains: list[str],
                   keep_ranges: list[Range],
                   target_chain: str | None = None) -> TrimResult:
    src_model = structure[0]
    new = gemmi.Structure()
    new.name = structure.name
    new.cell = structure.cell
    new.spacegroup_hm = structure.spacegroup_hm
    model = gemmi.Model("1")

    # Deterministic clean naming: target chain -> 'A', the rest alphabetical.
    keep_set = set(keep_chains)
    ordered = sorted(keep_set)
    if target_chain in keep_set:
        ordered = [target_chain] + [c for c in ordered if c != target_chain]
    rename = {old: _CLEAN_NAMES[i] if i < len(_CLEAN_NAMES) else f"X{i}"
              for i, old in enumerate(ordered)}

    kept_counts: dict[str, int] = {}
    warnings: list[str] = []

    for old in ordered:
        chain = None
        for c in src_model:
            if c.name == old:
                chain = c
                break
        if chain is None:
            continue
        poly_keys = _polymer_keys(chain)
        new_name = rename[old]
        nc = gemmi.Chain(new_name)
        kept_nums: list[int] = []
        for res in chain:
            key = (res.seqid.num, res.seqid.icode)
            if key not in poly_keys:
                continue  # drop ligands/glycans/waters
            if not _in_ranges(res.seqid.num, keep_ranges):
                continue
            nc.add_residue(res.clone())
            kept_nums.append(res.seqid.num)
        if len(nc) > 0:
            model.add_chain(nc)
            kept_counts[new_name] = len(nc)
            warnings.extend(_fragment_warnings(new_name, kept_nums))

    new.add_model(model)
    new.setup_entities()
    if not kept_counts:
        warnings.append("trim produced no residues; check keep_chains/keep_ranges")

    target_new = rename.get(target_chain, sorted(kept_counts)[0] if kept_counts else "A")
    return TrimResult(structure=new, kept_chains=sorted(kept_counts),
                      kept_counts=kept_counts, target_chain=target_new,
                      rename=rename, warnings=warnings)
