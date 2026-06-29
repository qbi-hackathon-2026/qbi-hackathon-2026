"""ESMFold ECD prediction via the public ESM Atlas API.

Only the extracellular domain (ECD) sequence is submitted — concatenated if
the ECD spans multiple disjoint TOPO_DOM ranges. Residue numbering in the
returned PDB starts at 1 and covers only the submitted ECD residues.

Public API limit: 400 residues. Proteins with longer ECDs raise ValueError.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import requests

from .cache import cache_dir
from .topology import Range

ESM_API = "https://api.esmatlas.com/foldSequence/v1/pdb/"
ESM_MAX_LEN = 400


@dataclass
class ESMFoldModel:
    accession: str
    ecd_ranges: list[Range]
    ecd_length: int
    pdb_bytes: bytes
    mean_plddt: float
    warnings: list[str] = field(default_factory=list)

    def summary_fields(self) -> dict:
        return {
            "esmfold_available": True,
            "esmfold_ecd_length": self.ecd_length,
            "esmfold_ecd_ranges": [(r.start, r.end) for r in self.ecd_ranges],
            "esmfold_mean_plddt": round(self.mean_plddt, 1),
        }


def _extract_ecd_sequence(full_seq: str, ecd_ranges: list[Range]) -> str:
    parts = []
    for r in ecd_ranges:
        parts.append(full_seq[r.start - 1 : r.end])
    return "".join(parts)


def _mean_plddt_from_pdb(pdb_bytes: bytes) -> float:
    values = []
    for line in pdb_bytes.decode(errors="replace").splitlines():
        if line.startswith(("ATOM", "HETATM")):
            try:
                values.append(float(line[60:66]))
            except (ValueError, IndexError):
                pass
    return sum(values) / len(values) if values else 0.0


def fetch_esmfold(accession: str, full_sequence: str,
                  ecd_ranges: list[Range]) -> ESMFoldModel:
    """Fold the ECD sequence via the public ESM Atlas API.

    Raises ValueError if no ECD ranges, ECD > ESM_MAX_LEN, or API fails.
    """
    if not ecd_ranges:
        raise ValueError("no ECD ranges — ESMFold requires an extracellular domain")

    ecd_seq = _extract_ecd_sequence(full_sequence, ecd_ranges)
    if len(ecd_seq) > ESM_MAX_LEN:
        raise ValueError(
            f"ECD is {len(ecd_seq)} residues, exceeding the "
            f"{ESM_MAX_LEN}-residue public ESMFold API limit")

    range_tag = "_".join(f"{r.start}-{r.end}" for r in ecd_ranges)
    cache_path = cache_dir() / "esmfold" / f"ESM-{accession}-ECD-{range_tag}.pdb"
    if cache_path.exists():
        pdb_bytes = cache_path.read_bytes()
    else:
        resp = requests.post(
            ESM_API, data=ecd_seq,
            headers={"Content-Type": "application/x-www-form-urlencoded"},
            timeout=120,
        )
        resp.raise_for_status()
        pdb_bytes = resp.content
        cache_path.parent.mkdir(parents=True, exist_ok=True)
        cache_path.write_bytes(pdb_bytes)

    mean_plddt = _mean_plddt_from_pdb(pdb_bytes)
    return ESMFoldModel(
        accession=accession,
        ecd_ranges=ecd_ranges,
        ecd_length=len(ecd_seq),
        pdb_bytes=pdb_bytes,
        mean_plddt=mean_plddt,
        warnings=[f"ESMFold ECD prediction for {accession} "
                  f"({len(ecd_seq)} aa; mean pLDDT {mean_plddt:.0f})"],
    )
