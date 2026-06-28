"""AlphaFold fallback: fetch a predicted model when PDBe/RCSB have no usable
experimental structure, and build the same objects the pipeline expects so the
rest of the flow (chain classification, avoid set, interface, trim, hotspots,
outputs) runs unchanged.

Design notes
------------
- AlphaFold author numbering == UniProt canonical residue number (single
  fragment, chain "A"). The SIFTS mapping is therefore the identity, built here
  directly instead of from PDBe.
- There are no partner chains, so the pipeline takes its apo / surface-exposed
  hotspot path automatically.
- pLDDT masking: per-residue confidence is stored in the CA B-factor column.
  Residues <= 50 are treated as unobserved, 50-70 as low-confidence; both are
  reported so the pipeline can fold them into the avoid set.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

import gemmi

from .cache import cache_dir, get_json
from .sifts import ResidueMap, SiftsMapping
from .structures import Candidate, StructureChoice
from .structio import LoadedStructure
from .topology import Range

AF_API = "https://alphafold.ebi.ac.uk/api/prediction/{acc}"
AF_FILES = "https://alphafold.ebi.ac.uk/files/AF-{acc}-F1-model_v{ver}.cif"

PLDDT_LOW = 70.0      # 50-70 -> low-confidence, added to avoid
PLDDT_UNOBS = 50.0    # <=50 -> treated as unobserved, added to avoid
AF_CHAIN = "A"


@dataclass
class AlphaFoldModel:
    """Everything the pipeline needs to treat an AF model like a chosen structure."""
    accession: str
    version: int
    model_url: str
    structure: gemmi.Structure          # gemmi structure, chain "A"
    choice: StructureChoice
    loaded: LoadedStructure
    mapping: SiftsMapping
    low_plddt_unp: list[int] = field(default_factory=list)    # 50-70
    unobs_plddt_unp: list[int] = field(default_factory=list)  # <=50
    ecd_mean_plddt: float = 0.0
    warnings: list[str] = field(default_factory=list)

    def summary_fields(self) -> dict:
        return {
            "alphafold_fallback": True,
            "alphafold_version": self.version,
            "alphafold_model_url": self.model_url,
            "alphafold_ecd_mean_plddt": round(self.ecd_mean_plddt, 1),
            "alphafold_n_low_plddt": len(self.low_plddt_unp),
            "alphafold_n_unobserved": len(self.unobs_plddt_unp),
        }


def _residue_plddt(res: gemmi.Residue) -> float:
    ca = next((a for a in res if a.name == "CA"), None)
    if ca is not None:
        return ca.b_iso
    atoms = list(res)
    return sum(a.b_iso for a in atoms) / len(atoms) if atoms else 0.0


def fetch_alphafold(accession: str, ecd_ranges: list[Range],
                    *, version: Optional[int] = None) -> AlphaFoldModel:
    """Download the AlphaFold model for `accession` and wrap it for the pipeline.

    Raises ValueError if AlphaFold has no entry for the accession.
    """
    meta = get_json(AF_API.format(acc=accession), allow_404=True)
    if not meta:
        raise ValueError(
            f"no experimental structure and no AlphaFold model for {accession}")
    entry = meta[0]
    ver = version or int(entry.get("latestVersion", 4))
    cif_url = entry.get("cifUrl") or AF_FILES.format(acc=accession, ver=ver)

    # Cache the CIF on disk (same cache root the rest of the engine uses).
    cif_path = cache_dir() / "alphafold" / f"AF-{accession}-F1-model_v{ver}.cif"
    if not cif_path.exists():
        import requests
        cif_path.parent.mkdir(parents=True, exist_ok=True)
        resp = requests.get(cif_url, timeout=60)
        resp.raise_for_status()
        cif_path.write_bytes(resp.content)

    structure = gemmi.read_structure(str(cif_path))
    structure.setup_entities()
    structure.spacegroup_hm = "P 1"
    # Normalise the (single) polymer chain to author id "A".
    model = structure[0]
    if model and model[0].name != AF_CHAIN:
        model[0].name = AF_CHAIN

    chain = model[AF_CHAIN]

    # Identity SIFTS mapping (auth == UniProt), all residues observed.
    ecd_set: set[int] = set()
    for r in ecd_ranges:
        ecd_set.update(range(r.start, r.end + 1))

    residues: list[ResidueMap] = []
    ecd_plddts: list[float] = []
    low_unp: list[int] = []
    unobs_unp: list[int] = []
    for res in chain:
        num = res.seqid.num
        residues.append(ResidueMap(unp=num, pdb_num=num, icode="",
                                   chain=AF_CHAIN, observed=True))
        if num in ecd_set:
            p = _residue_plddt(res)
            ecd_plddts.append(p)
            if p <= PLDDT_UNOBS:
                unobs_unp.append(num)
            elif p < PLDDT_LOW:
                low_unp.append(num)

    mapping = SiftsMapping(pdb_id=f"AF-{accession}", accession=accession,
                           residues=residues)

    # Synthetic candidate / choice so downstream code reads uniform fields.
    span = (residues[0].unp, residues[-1].unp) if residues else (0, 0)
    cand = Candidate(
        pdb_id=f"AF-{accession}", chain_id=AF_CHAIN,
        unp_start=span[0], unp_end=span[1], coverage=1.0, resolution=None,
        method="AlphaFold (predicted)", tax_id=None,
        observed_regions=[span], ecd_coverage=1.0, completeness=1.0,
        partner_tier=0, predicted=True, method_label="predicted",
    )
    choice = StructureChoice(
        candidates=[cand], chosen=cand, target_chain=AF_CHAIN,
        partner_chains=[], antibody_chains=[], apo_fallback=True,
        reasons=[f"AlphaFold model AF-{accession} v{ver} "
                 "(no experimental structure with adequate ECD coverage)"],
    )
    loaded = LoadedStructure(
        analysis=structure, display=structure, pdb_id=f"AF-{accession}",
        assembly="protomer", assembly_applied=False, warnings=[],
    )

    mean_plddt = sum(ecd_plddts) / len(ecd_plddts) if ecd_plddts else 0.0
    warnings = [f"AlphaFold model used for {accession} "
                f"(v{ver}; mean ECD pLDDT {mean_plddt:.0f})"]
    if mean_plddt and mean_plddt < PLDDT_LOW:
        warnings.append(
            f"mean ECD pLDDT {mean_plddt:.0f} is below {PLDDT_LOW:.0f}; the "
            "AlphaFold model may be unreliable for this extracellular domain")

    return AlphaFoldModel(
        accession=accession, version=ver, model_url=cif_url,
        structure=structure, choice=choice, loaded=loaded, mapping=mapping,
        low_plddt_unp=sorted(low_unp), unobs_plddt_unp=sorted(unobs_unp),
        ecd_mean_plddt=mean_plddt, warnings=warnings,
    )
