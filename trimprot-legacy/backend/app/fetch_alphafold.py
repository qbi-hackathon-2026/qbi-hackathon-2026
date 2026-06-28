"""
fetch_alphafold
===============
AlphaFold structure fetcher with pLDDT-based avoid-set masking.

Implements the AlphaFold fallback path:

    When pdb_search finds no PDB entry with adequate ECD coverage,
    download the AlphaFold model, use it as the design target, and apply
    pLDDT masking to extend the avoid set.

Hard correctness rules
-----------------------
Rule 1  Author numbering — AlphaFold auth_seq_id == UniProt canonical residue
        number for single-domain entries.  No renumbering is ever performed.
        Insertion codes are always preserved (empty string for AF models).

Rule 3  Unobserved residues → avoid.  Residues with pLDDT < PLDDT_UNOBS are
        the AlphaFold analogue of missing electron density; they are added to
        the avoid set with reason "alphafold_unobserved_plddt<50".

pLDDT thresholds (AlphaFold's own confidence bands)
-----------------------------------------------------
  > 90   Very high  — use freely as hotspot candidates
  70–90  Confident  — use freely
  50–70  Low        → avoid("alphafold_low_plddt<70"); reason logged
  ≤ 50   Very low   → avoid("alphafold_unobserved_plddt<50"); treat as
                       unobserved (Rule 3)

Avoid contributions cover ECD residues only.
"""

from __future__ import annotations

import logging
import os
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

import gemmi
import requests

log = logging.getLogger(__name__)

_AF_API   = "https://alphafold.ebi.ac.uk/api/prediction"
_AF_FILES = "https://alphafold.ebi.ac.uk/files"

PLDDT_LOW_THRESHOLD   = 70.0
PLDDT_UNOBS_THRESHOLD = 50.0
PLDDT_HIGH_THRESHOLD  = 90.0


@dataclass
class AvoidEntry:
    """A single residue excluded from hotspot candidacy."""
    chain: str
    auth_seq_id: int
    icode: str
    reason: str
    uniprot_residue: Optional[int] = None
    plddt: Optional[float] = None

    def key(self) -> tuple[str, int, str]:
        return (self.chain, self.auth_seq_id, self.icode)


@dataclass
class ResidueConfidence:
    """Per-residue pLDDT record for one AlphaFold chain."""
    auth_seq_id: int
    icode: str
    residue_name: str
    plddt: float
    in_ecd: bool
    confidence_band: str
    avoid_reason: Optional[str]


@dataclass
class AlphaFoldResult:
    """Return value of get_alphafold_structure()."""
    accession: str
    version: int
    structure: gemmi.Structure
    model_url: str
    chain_id: str
    residues: list[ResidueConfidence]
    avoid_contributions: list[AvoidEntry]
    ecd_mean_plddt: float
    ecd_n_low_plddt: int
    ecd_n_unobserved: int
    low_plddt_ranges: list[tuple[int, int]]
    unobs_ranges: list[tuple[int, int]]
    alphafold_fallback: bool = True
    warning: Optional[str] = None

    def summary_fields(self) -> dict:
        return {
            "alphafold_fallback":         True,
            "alphafold_version":          self.version,
            "alphafold_model_url":        self.model_url,
            "alphafold_ecd_mean_plddt":   round(self.ecd_mean_plddt, 2),
            "alphafold_ecd_n_low_plddt":  self.ecd_n_low_plddt,
            "alphafold_ecd_n_unobserved": self.ecd_n_unobserved,
            "alphafold_low_plddt_ranges": [list(r) for r in self.low_plddt_ranges],
            "alphafold_unobs_ranges":     [list(r) for r in self.unobs_ranges],
            "alphafold_warning":          self.warning,
        }


def _confidence_band(plddt: float) -> str:
    if plddt > PLDDT_HIGH_THRESHOLD:
        return "very_high"
    if plddt >= PLDDT_LOW_THRESHOLD:
        return "confident"
    if plddt > PLDDT_UNOBS_THRESHOLD:
        return "low"
    return "very_low"


def _avoid_reason(plddt: float) -> Optional[str]:
    if plddt <= PLDDT_UNOBS_THRESHOLD:
        return f"alphafold_unobserved_plddt<{PLDDT_UNOBS_THRESHOLD:.0f}"
    if plddt < PLDDT_LOW_THRESHOLD:
        return f"alphafold_low_plddt<{PLDDT_LOW_THRESHOLD:.0f}"
    return None


def _merge_ranges(positions: list[int]) -> list[tuple[int, int]]:
    if not positions:
        return []
    positions = sorted(set(positions))
    ranges: list[tuple[int, int]] = []
    start = end = positions[0]
    for pos in positions[1:]:
        if pos == end + 1:
            end = pos
        else:
            ranges.append((start, end))
            start = end = pos
    ranges.append((start, end))
    return ranges


def _residue_plddt(residue: gemmi.Residue) -> float:
    atoms = list(residue)
    if not atoms:
        return 0.0
    ca = next((a for a in atoms if a.name == "CA"), None)
    if ca is not None:
        return ca.b_iso
    return sum(a.b_iso for a in atoms) / len(atoms)


def _fetch_metadata(accession: str, timeout: int = 15) -> dict:
    url = f"{_AF_API}/{accession}"
    try:
        r = requests.get(url, timeout=timeout, headers={"Accept": "application/json"})
        r.raise_for_status()
        hits = r.json()
        if not hits:
            raise RuntimeError(f"AlphaFold has no entry for accession {accession!r}.")
        return hits[0]
    except requests.HTTPError as exc:
        if exc.response is not None and exc.response.status_code == 404:
            raise RuntimeError(f"Accession {accession!r} not found in AlphaFold DB.") from exc
        raise RuntimeError(f"AlphaFold API error for {accession!r}: {exc}") from exc
    except requests.RequestException as exc:
        raise RuntimeError(f"AlphaFold API unreachable: {exc}") from exc


def _download_cif(url: str, dest: Path, timeout: int = 60) -> None:
    log.info("Downloading AlphaFold CIF: %s", url)
    with requests.get(url, stream=True, timeout=timeout) as r:
        r.raise_for_status()
        with open(dest, "wb") as f:
            for chunk in r.iter_content(chunk_size=1 << 16):
                f.write(chunk)
    log.info("Saved CIF → %s (%d bytes)", dest, dest.stat().st_size)


def _ecd_set_from_ranges(ecd_ranges: list[tuple[int, int]]) -> frozenset[int]:
    positions: set[int] = set()
    for start, end in ecd_ranges:
        positions.update(range(start, end + 1))
    return frozenset(positions)


def _extract_residue_confidences(
    structure: gemmi.Structure,
    ecd_set: frozenset[int],
    chain_id: str = "A",
) -> list[ResidueConfidence]:
    model = structure[0]
    try:
        chain = model[chain_id]
    except (KeyError, ValueError):
        raise RuntimeError(
            f"Chain '{chain_id}' not found in AlphaFold structure "
            f"(chains present: {[c.name for c in model]})"
        )

    records: list[ResidueConfidence] = []
    for residue in chain:
        auth_seq_id = residue.seqid.num
        icode = residue.seqid.icode.strip()
        plddt = _residue_plddt(residue)
        band = _confidence_band(plddt)
        in_ecd = auth_seq_id in ecd_set
        reason = _avoid_reason(plddt) if in_ecd else None
        records.append(ResidueConfidence(
            auth_seq_id=auth_seq_id,
            icode=icode,
            residue_name=residue.name,
            plddt=plddt,
            in_ecd=in_ecd,
            confidence_band=band,
            avoid_reason=reason,
        ))
    return records


def _build_avoid_contributions(
    residues: list[ResidueConfidence],
    chain_id: str,
    accession: str,
) -> list[AvoidEntry]:
    entries: list[AvoidEntry] = []
    for rec in residues:
        if not rec.in_ecd or rec.avoid_reason is None:
            continue
        entries.append(AvoidEntry(
            chain=chain_id,
            auth_seq_id=rec.auth_seq_id,
            icode=rec.icode,
            reason=rec.avoid_reason,
            uniprot_residue=rec.auth_seq_id,
            plddt=rec.plddt,
        ))
    return entries


def get_alphafold_structure(
    accession: str,
    ecd_ranges: list[tuple[int, int]],
    *,
    cache_dir: Path = Path("cache") / "alphafold",
    plddt_warn_threshold: float = PLDDT_LOW_THRESHOLD,
    plddt_unobs_threshold: float = PLDDT_UNOBS_THRESHOLD,
    version: Optional[int] = None,
    chain_id: str = "A",
    timeout: int = 60,
) -> AlphaFoldResult:
    """
    Download (or load from cache) the AlphaFold model for *accession* and
    apply pLDDT masking to produce avoid-set contributions.

    ecd_ranges must not be empty — call get_extracellular_domain() first.
    """
    if not ecd_ranges:
        raise RuntimeError(
            f"get_alphafold_structure({accession!r}): ecd_ranges must not be "
            "empty — call get_extracellular_ranges() first."
        )

    _low = plddt_warn_threshold
    _unob = plddt_unobs_threshold
    if _unob >= _low:
        raise ValueError(
            f"plddt_unobs_threshold ({_unob}) must be < plddt_warn_threshold ({_low})."
        )

    meta = _fetch_metadata(accession, timeout=min(timeout, 15))
    latest_version = meta.get("latestVersion", 4)
    chosen_version = version if version is not None else latest_version
    cif_url = (
        meta.get("cifUrl")
        if chosen_version == latest_version
        else f"{_AF_FILES}/AF-{accession}-F1-model_v{chosen_version}.cif"
    )
    log.info("AlphaFold %s: version=%d url=%s", accession, chosen_version, cif_url)

    cache_dir = Path(cache_dir)
    cache_dir.mkdir(parents=True, exist_ok=True)
    cif_path = cache_dir / f"AF-{accession}-F1-model_v{chosen_version}.cif"

    if cif_path.exists():
        log.info("AlphaFold cache hit: %s", cif_path)
    else:
        _download_cif(cif_url, cif_path, timeout=timeout)

    structure = gemmi.read_structure(str(cif_path))
    structure.setup_entities()

    ecd_set = _ecd_set_from_ranges(ecd_ranges)
    residues = _extract_residue_confidences(structure, ecd_set, chain_id)

    avoid_contributions: list[AvoidEntry] = []
    for rec in residues:
        if not rec.in_ecd:
            continue
        reason: Optional[str] = None
        if rec.plddt <= _unob:
            reason = f"alphafold_unobserved_plddt<{_unob:.0f}"
        elif rec.plddt < _low:
            reason = f"alphafold_low_plddt<{_low:.0f}"
        if reason:
            avoid_contributions.append(AvoidEntry(
                chain=chain_id,
                auth_seq_id=rec.auth_seq_id,
                icode=rec.icode,
                reason=reason,
                uniprot_residue=rec.auth_seq_id,
                plddt=rec.plddt,
            ))

    ecd_records = [r for r in residues if r.in_ecd]
    ecd_plddts  = [r.plddt for r in ecd_records]
    ecd_mean    = sum(ecd_plddts) / len(ecd_plddts) if ecd_plddts else 0.0
    n_low  = sum(1 for p in ecd_plddts if _unob < p < _low)
    n_unob = sum(1 for p in ecd_plddts if p <= _unob)

    low_positions  = [r.auth_seq_id for r in ecd_records if _unob < r.plddt < _low]
    unob_positions = [r.auth_seq_id for r in ecd_records if r.plddt <= _unob]
    low_ranges  = _merge_ranges(low_positions)
    unob_ranges = _merge_ranges(unob_positions)

    for s, e in low_ranges:
        log.warning("AlphaFold %s: low pLDDT ECD region %d–%d (< %.0f); added to avoid set.",
                    accession, s, e, _low)
    for s, e in unob_ranges:
        log.warning("AlphaFold %s: very low pLDDT ECD region %d–%d (≤ %.0f); treated as unobserved.",
                    accession, s, e, _unob)

    warning: Optional[str] = None
    if ecd_mean < _low:
        warning = (
            f"Mean ECD pLDDT is {ecd_mean:.1f}, below the {_low:.0f} confidence threshold. "
            "The AlphaFold model may not be reliable for this extracellular domain."
        )
        log.warning("AlphaFold %s: %s", accession, warning)

    return AlphaFoldResult(
        accession=accession,
        version=chosen_version,
        structure=structure,
        model_url=cif_url,
        chain_id=chain_id,
        residues=residues,
        avoid_contributions=avoid_contributions,
        ecd_mean_plddt=ecd_mean,
        ecd_n_low_plddt=n_low,
        ecd_n_unobserved=n_unob,
        low_plddt_ranges=low_ranges,
        unobs_ranges=unob_ranges,
        alphafold_fallback=True,
        warning=warning,
    )
