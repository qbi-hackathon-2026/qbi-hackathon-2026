"""Structure search + ranking.

PDBe graph-api `best_structures/{accession}` is the ranked spine (already quality
ordered and UniProt-mapped). Each candidate is enriched via the RCSB Data API
(GraphQL) with polymer-entity descriptions, chains, source organism, and ligands,
then scored deterministically. The chosen structure favours: human, ECD coverage,
an antibody/partner-bound complex, good resolution. If nothing is partner-bound,
the best apo structure is chosen and apo_fallback is set.
"""
from __future__ import annotations

import re
from dataclasses import dataclass, field
from typing import Optional

from .cache import get_json, post_json
from .sifts import available_chains
from .topology import Range

PDBE_BEST = "https://www.ebi.ac.uk/pdbe/graph-api/uniprot/best_structures/{acc}"
PDBE_PFAM = "https://www.ebi.ac.uk/pdbe/api/mappings/pfam/{pdb}"
RCSB_GRAPHQL = "https://data.rcsb.org/graphql"


def pfam_domains(pdb_id: str, chain: str) -> list[dict]:
    """Pfam domains for a chain via PDBe, in PDB author numbering (best-effort).

    Returns [{name, start, end}] sorted by start. Author residue numbers match the
    trimmed structure, so the frontend can label them directly.
    """
    try:
        data = get_json(PDBE_PFAM.format(pdb=pdb_id.lower()), allow_404=True)
    except Exception:
        return []
    fam = (data.get(pdb_id.lower()) or {}).get("Pfam", {}) if data else {}
    out: list[dict] = []
    for fid, info in fam.items():
        name = info.get("name") or info.get("identifier") or fid
        for seg in info.get("mappings", []) or []:
            if seg.get("chain_id") != chain:
                continue
            s = (seg.get("start") or {}).get("author_residue_number")
            e = (seg.get("end") or {}).get("author_residue_number")
            if s is None or e is None:
                continue
            out.append({"name": name, "start": int(s), "end": int(e)})
    out.sort(key=lambda d: d["start"])
    return out

_ANTIBODY_RE = re.compile(
    r"\b(fab|fv|scfv|nanobody|vhh|antibody|immunoglobulin|"
    r"heavy chain|light chain|igg|sybody)\b", re.I)

_RCSB_QUERY = """
query($ids:[String!]!){
  entries(entry_ids:$ids){
    rcsb_id
    rcsb_entry_info{ structure_determination_methodology resolution_combined }
    exptl{ method }
    polymer_entities{
      rcsb_polymer_entity{ pdbx_description }
      entity_poly{ rcsb_sample_sequence_length }
      rcsb_polymer_entity_container_identifiers{ auth_asym_ids }
      rcsb_entity_source_organism{ ncbi_taxonomy_id }
    }
    nonpolymer_entities{
      nonpolymer_comp{ chem_comp{ id name } }
    }
  }
}
"""

# Minimum partner length to count as a real binding partner (drops poly-Ala /
# "undefined peptide" crystallisation artifacts).
MIN_PARTNER_LENGTH = 20
_JUNK_PARTNER_RE = re.compile(r"undefined|unknown|uncharacteri[sz]ed|"
                              r"poly-?(ala|gly|unk)|modeled as", re.I)


@dataclass
class Candidate:
    pdb_id: str
    chain_id: str           # author chain id of the target
    unp_start: int
    unp_end: int
    coverage: float
    resolution: Optional[float]
    method: str
    tax_id: Optional[int]
    observed_regions: list[tuple[int, int]] = field(default_factory=list)
    partner_chains: list[str] = field(default_factory=list)
    antibody_chains: list[str] = field(default_factory=list)
    antibody_partner: bool = False
    ligands: list[str] = field(default_factory=list)
    enriched: bool = False
    # ladder metrics (computed in search_structures)
    ecd_coverage: float = 0.0              # fraction of the ECD spanned (mapped)
    completeness: float = 1.0              # observed / mapped residues
    partner_tier: int = 0                  # 2 antibody, 1 receptor/ligand, 0 apo
    methodology: Optional[str] = None      # RCSB structure_determination_methodology
    predicted: bool = False
    method_label: str = ""                 # "X-ray" / "cryo-EM" / "NMR" / "predicted"
    reasons: list[str] = field(default_factory=list)

    @property
    def method_display(self) -> str:
        if self.predicted:
            return "predicted"
        res = f" {self.resolution}Å" if self.resolution is not None else ""
        return f"{self.method_label or self.method}{res}"


@dataclass
class StructureChoice:
    candidates: list[Candidate]
    chosen: Candidate
    target_chain: str
    partner_chains: list[str]
    antibody_chains: list[str]
    apo_fallback: bool
    reasons: list[str]


def _overlap(a: Range, s: int, e: int) -> int:
    return max(0, min(a.end, e) - max(a.start, s) + 1)


def fetch_best_structures(accession: str) -> list[Candidate]:
    data = get_json(PDBE_BEST.format(acc=accession))
    items = data.get(accession, []) if data else []
    out = []
    for it in items:
        observed = [(int(o["unp_start"]), int(o["unp_end"]))
                    for o in (it.get("observed_regions") or [])
                    if o.get("unp_start") is not None and o.get("unp_end") is not None]
        out.append(Candidate(
            pdb_id=str(it["pdb_id"]).lower(),
            chain_id=str(it["chain_id"]),
            unp_start=int(it.get("unp_start", 0) or 0),
            unp_end=int(it.get("unp_end", 0) or 0),
            coverage=float(it.get("coverage", 0) or 0),
            resolution=(float(it["resolution"]) if it.get("resolution") else None),
            method=str(it.get("experimental_method", "")),
            tax_id=(int(it["tax_id"]) if it.get("tax_id") else None),
            observed_regions=observed,
        ))
    return out


def _apply_entry(cand: Candidate, entry: dict) -> None:
    partner_chains: list[str] = []
    antibody_chains: list[str] = []
    for pe in entry.get("polymer_entities") or []:
        ids = ((pe.get("rcsb_polymer_entity_container_identifiers") or {})
               .get("auth_asym_ids") or [])
        desc = ((pe.get("rcsb_polymer_entity") or {}).get("pdbx_description") or "")
        length = ((pe.get("entity_poly") or {}).get("rcsb_sample_sequence_length") or 0)
        if cand.chain_id in ids:
            continue  # the target's own entity
        # Skip crystallisation-artifact / undefined short peptides.
        if length and length < MIN_PARTNER_LENGTH:
            continue
        if _JUNK_PARTNER_RE.search(desc):
            continue
        partner_chains.extend(ids)
        if _ANTIBODY_RE.search(desc):
            antibody_chains.extend(ids)
    ligands = []
    for ne in entry.get("nonpolymer_entities") or []:
        comp = ((ne.get("nonpolymer_comp") or {}).get("chem_comp") or {})
        if comp.get("id"):
            ligands.append(comp["id"])

    cand.partner_chains = sorted(set(partner_chains))
    cand.antibody_chains = sorted(set(antibody_chains))
    cand.antibody_partner = bool(antibody_chains)
    cand.ligands = sorted(set(ligands))

    # method + experimental/computational flag for tiering
    info = entry.get("rcsb_entry_info") or {}
    cand.methodology = info.get("structure_determination_methodology")
    exptl = entry.get("exptl") or []
    if exptl and exptl[0].get("method") and not cand.method:
        cand.method = exptl[0]["method"]
    if cand.resolution is None:
        rc = info.get("resolution_combined")
        if rc:
            cand.resolution = float(rc[0])
    cand.enriched = True


# Method label/provenance. The real axis is EXPERIMENTAL vs PREDICTED — X-ray
# (crystallography) and cryo-EM are both top experimental methods, not separate
# categories. Returns (predicted_flag, label).
def method_label(method: str, methodology: Optional[str]) -> tuple[bool, str]:
    m = (method or "").upper()
    meth = (methodology or "").lower()
    predicted = ("comput" in meth) or any(
        k in m for k in ("ALPHAFOLD", "ESMFOLD", "ROSETTAFOLD", "PREDICTED",
                         "COMPUTATIONAL", "DEEP LEARNING", "CSM"))
    if predicted:
        return True, "predicted"
    if "X-RAY" in m or "DIFFRACTION" in m or "CRYSTAL" in m:
        return False, "X-ray"
    if "ELECTRON" in m or "CRYO" in m or m == "EM":
        return False, "cryo-EM"
    if "NMR" in m:
        return False, "NMR"
    return False, (method or "experimental")


def ecd_coverage(cand: Candidate, ecd_ranges: list[Range]) -> float:
    """Fraction of the extracellular domain spanned by the candidate's mapped range."""
    total = sum(r.end - r.start + 1 for r in ecd_ranges) or 1
    covered = sum(_overlap(r, cand.unp_start, cand.unp_end) for r in ecd_ranges)
    return covered / total


def completeness(cand: Candidate) -> float:
    """Observed residues / mapped residues (1.0 when no gap info -> assume complete)."""
    mapped = cand.unp_end - cand.unp_start + 1
    if mapped <= 0:
        return 0.0
    if not cand.observed_regions:
        return 1.0
    obs = sum(e - s + 1 for s, e in cand.observed_regions)
    return min(1.0, obs / mapped)


def partner_tier(cand: Candidate) -> int:
    """2 = antibody-bound, 1 = receptor/ligand-bound, 0 = apo."""
    if cand.antibody_partner:
        return 2
    if cand.partner_chains:
        return 1
    return 0


def _res_sort_key(cand: Candidate) -> float:
    """Resolution for the final tiebreak; NMR/predicted/unknown sort worst."""
    if cand.resolution is not None and not cand.predicted and cand.method_label != "NMR":
        return cand.resolution
    return float("inf")


def enrich_rcsb_batch(candidates: list[Candidate]) -> None:
    """Enrich all candidates with one batched RCSB `entries` query (best-effort)."""
    ids = sorted({c.pdb_id.upper() for c in candidates})
    if not ids:
        return
    by_id: dict[str, dict] = {}
    # RCSB caps batch size; chunk to be safe.
    for i in range(0, len(ids), 50):
        chunk = ids[i:i + 50]
        try:
            d = post_json(RCSB_GRAPHQL, {"query": _RCSB_QUERY,
                                         "variables": {"ids": chunk}})
            for e in (d.get("data") or {}).get("entries") or []:
                by_id[str(e.get("rcsb_id", "")).upper()] = e
        except Exception:  # network/enrichment is best-effort
            continue
    for c in candidates:
        entry = by_id.get(c.pdb_id.upper())
        if entry:
            _apply_entry(c, entry)


# Default ladder thresholds.
MIN_ECD_COVERAGE = 0.40          # Gate 1: must cover at least this fraction of ECD
MAX_USABLE_RESOLUTION = 9.0      # Gate 2: drop pathological resolution (e.g. 25Å)
COVERAGE_TOL_RESIDUES = 4        # "tied" coverage: within this many ECD residues
COMPLETENESS_TOL = 0.04          # "tied" completeness: within this fraction


def _ladder_sort_key(c: Candidate, ecd_total: int):
    """Full-ladder ordering for display/ranking of the candidate list."""
    return (-c.partner_tier, -round(c.ecd_coverage, 6), -round(c.completeness, 6),
            1 if c.predicted else 0, _res_sort_key(c), c.pdb_id, c.chain_id)


def search_structures(accession: str, ecd_ranges: list[Range], *,
                      prefer_antibody: bool = False,
                      min_ecd_coverage: float = MIN_ECD_COVERAGE,
                      max_resolution: float = MAX_USABLE_RESOLUTION
                      ) -> StructureChoice:
    """Pick the best PDB by a STRICT PRIORITY LADDER (not a weighted sum).

    Order: (1) ECD-coverage gate, (2) usability gate [SIFTS + resolution sanity],
    (3) best partner tier [antibody > receptor/ligand > apo], (4) most ECD
    coverage, (5) most complete (fewest missing), (6) FINAL tiebreaker only:
    experimental over predicted, then better resolution. Method/resolution can
    never outrank a biological criterion — they only separate structures already
    tied (within tolerance) on coverage and completeness.
    """
    candidates = fetch_best_structures(accession)
    if not candidates:
        raise ValueError(f"no PDB structures found for {accession}")

    enrich_rcsb_batch(candidates)
    ecd_total = sum(r.end - r.start + 1 for r in ecd_ranges) or 1
    for c in candidates:
        c.ecd_coverage = ecd_coverage(c, ecd_ranges)
        c.completeness = completeness(c)
        c.partner_tier = partner_tier(c)
        c.predicted, c.method_label = method_label(c.method, c.methodology)

    candidates.sort(key=lambda c: _ladder_sort_key(c, ecd_total))

    sifts_idx = available_chains(accession)

    def has_sifts(c: Candidate) -> bool:
        if not sifts_idx:           # network failure: don't over-filter
            return True
        return c.chain_id in sifts_idx.get(c.pdb_id, set())

    chosen, _survivors, reasons = select_by_ladder(
        candidates, ecd_total, has_sifts=has_sifts,
        min_ecd_coverage=min_ecd_coverage, max_resolution=max_resolution,
        sifts_known=bool(sifts_idx))

    return StructureChoice(
        candidates=candidates,
        chosen=chosen,
        target_chain=chosen.chain_id,
        partner_chains=chosen.partner_chains,
        antibody_chains=chosen.antibody_chains,
        apo_fallback=chosen.partner_tier == 0,
        reasons=reasons,
    )


def select_by_ladder(candidates: list[Candidate], ecd_total: int, *,
                     has_sifts=lambda c: True,
                     min_ecd_coverage: float = MIN_ECD_COVERAGE,
                     max_resolution: float = MAX_USABLE_RESOLUTION,
                     sifts_known: bool = True):
    """Pure priority-ladder selection over candidates with metrics already set.

    Each candidate must have ecd_coverage, completeness, partner_tier, predicted,
    and method_label populated. Returns (chosen, survivors, reasons). No network.
    """
    reasons: list[str] = []

    # Gate 2 — usability: SIFTS numbering + sane resolution (NMR/predicted: None).
    usable = [c for c in candidates if has_sifts(c)
              and (c.resolution is None or c.resolution <= max_resolution)]
    pool = usable or [c for c in candidates if has_sifts(c)] or candidates

    # Gate 1 — ECD coverage (relax if it would empty the pool, e.g. CD44's HABD).
    gated = [c for c in pool if c.ecd_coverage >= min_ecd_coverage]
    if gated:
        survivors = gated
        reasons.append(f"ECD gate >= {min_ecd_coverage}: {len(gated)}/{len(pool)} pass")
    else:
        survivors = pool
        reasons.append(f"ECD gate >= {min_ecd_coverage} emptied the pool -> relaxed")

    # Step 2 — partner tier preference (antibody > receptor/ligand > apo).
    best_tier = max(c.partner_tier for c in survivors)
    survivors = [c for c in survivors if c.partner_tier == best_tier]
    tier_name = {2: "antibody-bound", 1: "receptor/ligand-bound", 0: "apo"}[best_tier]
    reasons.append(f"partner tier = {tier_name} ({len(survivors)} candidates)")

    # Step 3 — within the tier: most ECD coverage, then most complete (with a
    # tolerance, so trivial differences defer to the resolution tiebreaker).
    cov_tol = COVERAGE_TOL_RESIDUES / max(1, ecd_total)
    best_cov = max(c.ecd_coverage for c in survivors)
    survivors = [c for c in survivors if c.ecd_coverage >= best_cov - cov_tol]
    best_compl = max(c.completeness for c in survivors)
    survivors = [c for c in survivors if c.completeness >= best_compl - COMPLETENESS_TOL]
    reasons.append(f"coverage~{best_cov:.3f} (+/-{COVERAGE_TOL_RESIDUES}res), "
                   f"completeness~{best_compl:.3f} (+/-{COMPLETENESS_TOL}): "
                   f"{len(survivors)} tied")

    # Step 4 — FINAL tiebreaker only: experimental before predicted, then best
    # resolution, then deterministic id.
    survivors.sort(key=lambda c: (1 if c.predicted else 0, _res_sort_key(c),
                                  c.pdb_id, c.chain_id))
    chosen = survivors[0]

    head = f"chose {chosen.pdb_id} chain {chosen.chain_id} [{tier_name}, "
    head += "predicted model" if chosen.predicted else f"{chosen.method_label} {chosen.resolution}"
    head += f", ecd_cov={chosen.ecd_coverage:.2f}, completeness={chosen.completeness:.2f}]"
    reasons.insert(0, head)
    if chosen.predicted:
        reasons.append("WARNING: no experimental structure; predicted model used")
    if sifts_known and not has_sifts(chosen):
        reasons.append("WARNING: chosen structure lacks a SIFTS mapping")
    return chosen, survivors, reasons
