"""SIFTS UniProt<->PDB(auth) residue mapping via PDBe, cross-checked vs gemmi.

The PDBe `api/mappings/{accession}` endpoint gives segment-based mappings with
author residue numbers and insertion codes. Within a segment the mapping is a
constant offset (author_residue_number - unp), which is exactly the
mature-vs-precursor numbering shift the EGFR keystone test checks for.

We never assume a 1:1 identity mapping. Observed/unobserved status is decided
against the actual coordinates. When the deposited mmCIF carries
`_pdbx_sifts_xref_db` (gemmi `Residue.sifts_unp`), we cross-check and log
mismatches.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional

import gemmi

from .cache import get_json

PDBE_MAPPINGS = "https://www.ebi.ac.uk/pdbe/api/mappings/{acc}"


@dataclass(frozen=True)
class ResidueMap:
    unp: int
    pdb_num: int
    icode: str
    chain: str          # author chain id
    observed: bool = False


@dataclass
class SiftsMapping:
    pdb_id: str
    accession: str
    residues: list[ResidueMap] = field(default_factory=list)
    # raw segments: (chain, unp_start, unp_end, author_start, offset, icode_start)
    segments: list[dict] = field(default_factory=list)
    xref_mismatches: list[str] = field(default_factory=list)

    # ---- lookups -------------------------------------------------------
    def chains(self) -> list[str]:
        return sorted({r.chain for r in self.residues})

    def for_chain(self, chain: str) -> list[ResidueMap]:
        return [r for r in self.residues if r.chain == chain]

    def unp_to_pdb(self, chain: str) -> dict[int, tuple[int, str]]:
        return {r.unp: (r.pdb_num, r.icode) for r in self.for_chain(chain)}

    def pdb_to_unp(self, chain: str) -> dict[tuple[int, str], int]:
        return {(r.pdb_num, r.icode): r.unp for r in self.for_chain(chain)}

    def segment_offsets(self, chain: str) -> list[int]:
        """Per-segment author-minus-UniProt offsets actually realised in coords."""
        offs = []
        for r in self.for_chain(chain):
            if not r.icode:
                offs.append(r.pdb_num - r.unp)
        return offs

    def consistent_offset(self, chain: str) -> Optional[int]:
        offs = set(self.segment_offsets(chain))
        return next(iter(offs)) if len(offs) == 1 else None


def available_chains(accession: str) -> dict[str, set[str]]:
    """{pdb_id(lower) -> set(author chains)} that have a SIFTS mapping.

    One PDBe call covers every structure for the accession, so this is a cheap
    way to avoid choosing a structure that lacks SIFTS (e.g. very new depositions).
    Returns {} on network failure so the caller can degrade gracefully.
    """
    try:
        data = get_json(PDBE_MAPPINGS.format(acc=accession))
    except Exception:
        return {}
    pdb_map = (data.get(accession) or {}).get("PDB", {}) if data else {}
    out: dict[str, set[str]] = {}
    for pid, segs in pdb_map.items():
        out[pid.lower()] = {s.get("chain_id") for s in segs if s.get("chain_id")}
    return out


def _fetch_segments(pdb_id: str, accession: str) -> list[dict]:
    data = get_json(PDBE_MAPPINGS.format(acc=accession))
    pdb_map = (data.get(accession) or {}).get("PDB", {})
    # PDBe keys are lowercase pdb ids
    segs_raw = pdb_map.get(pdb_id.lower()) or pdb_map.get(pdb_id.upper()) or []
    out = []
    for s in segs_raw:
        try:
            unp_start = int(s["unp_start"])
            unp_end = int(s["unp_end"])
            label_start = int(s["start"]["residue_number"])
        except (KeyError, TypeError, ValueError):
            continue
        author_start = s.get("start", {}).get("author_residue_number")
        icode = s["start"].get("author_insertion_code", "") or ""
        out.append({
            "chain": s.get("chain_id"),
            "struct_asym_id": s.get("struct_asym_id"),
            "unp_start": unp_start,
            "unp_end": unp_end,
            "label_start": label_start,          # SIFTS residue_number (label seq)
            "author_start": (int(author_start) if author_start is not None else None),
            "icode_start": icode,
            # unp -> label is a constant shift within the segment
            "label_shift": label_start - unp_start,
        })
    return out


def get_sifts_mapping(pdb_id: str, accession: str,
                      structure: Optional[gemmi.Structure] = None,
                      target_chain: Optional[str] = None) -> SiftsMapping:
    """Build a per-residue UniProt<->PDB(auth) mapping.

    If `structure` is provided, observed/unobserved is decided against its
    coordinates and a gemmi `sifts_unp` cross-check is run when available.
    """
    segments = _fetch_segments(pdb_id, accession)
    label_maps = _label_to_auth(structure) if structure is not None else {}

    residues: list[ResidueMap] = []
    for seg in segments:
        chain = seg["chain"]
        shift = seg["label_shift"]                 # label = unp + shift
        lmap = label_maps.get(chain, {})

        # Author offset for the segment (auth_num - label) from an observed,
        # icode-free residue; used to place unobserved residues in author space.
        auth_off = None
        for lab in range(seg["unp_start"] + shift, seg["unp_end"] + shift + 1):
            if lab in lmap and not lmap[lab][1]:
                auth_off = lmap[lab][0] - lab
                break
        if auth_off is None and seg["author_start"] is not None:
            auth_off = seg["author_start"] - seg["label_start"]

        for unp in range(seg["unp_start"], seg["unp_end"] + 1):
            label = unp + shift
            if label in lmap:
                num, icode = lmap[label]
                observed = True
            elif auth_off is not None:
                num, icode, observed = label + auth_off, "", False
            else:
                continue  # cannot place in author numbering without guessing
            residues.append(ResidueMap(unp=unp, pdb_num=num, icode=icode,
                                       chain=chain, observed=observed))

    mapping = SiftsMapping(pdb_id=pdb_id.lower(), accession=accession,
                           residues=residues, segments=segments)
    if structure is not None:
        _crosscheck_xref(mapping, structure, target_chain)
    return mapping


def _label_to_auth(structure: gemmi.Structure
                   ) -> dict[str, dict[int, tuple[int, str]]]:
    """Per chain: SIFTS/label seq id -> (author number, insertion code)."""
    out: dict[str, dict[int, tuple[int, str]]] = {}
    model = structure[0]
    for chain in model:
        m: dict[int, tuple[int, str]] = {}
        for res in chain.get_polymer():
            if res.label_seq is not None:
                m[int(res.label_seq)] = (res.seqid.num, res.seqid.icode.strip())
        if m:
            out[chain.name] = m
    return out


def _crosscheck_xref(mapping: SiftsMapping, structure: gemmi.Structure,
                     target_chain: Optional[str]) -> None:
    """Compare PDBe-derived unp numbers against gemmi sifts_unp when present."""
    model = structure[0]
    # Build (chain, num, icode) -> unp from gemmi xref
    xref: dict[tuple[str, int, str], int] = {}
    for chain in model:
        for res in chain:
            acc_char, unp_num, _ = res.sifts_unp
            if acc_char and acc_char != "\x00" and unp_num:
                xref[(chain.name, res.seqid.num, res.seqid.icode.strip())] = int(unp_num)
    if not xref:
        return
    for r in mapping.residues:
        if target_chain and r.chain != target_chain:
            continue
        key = (r.chain, r.pdb_num, r.icode.strip())
        if key in xref and xref[key] != r.unp:
            mapping.xref_mismatches.append(
                f"chain {r.chain} auth {r.pdb_num}{r.icode}: "
                f"PDBe unp {r.unp} vs gemmi xref unp {xref[key]}"
            )
