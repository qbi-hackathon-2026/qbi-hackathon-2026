"""Assemble the structured design-prep summary described in the TrimProt spec."""
from __future__ import annotations

from typing import Optional


def compact_ranges(positions: list[int]) -> str:
    """Collapse a sorted list of residue positions into a compact range string,
    e.g. [56, 73, 128, 129, 130] -> "56, 73, 128-130". Used to keep the default
    (non-verbose) summary readable instead of dumping every residue individually.
    """
    positions = sorted(set(positions))
    if not positions:
        return "none"
    ranges = []
    start = prev = positions[0]
    for pos in positions[1:]:
        if pos == prev + 1:
            prev = pos
            continue
        ranges.append(f"{start}-{prev}" if start != prev else f"{start}")
        start = prev = pos
    ranges.append(f"{start}-{prev}" if start != prev else f"{start}")
    return ", ".join(ranges)


def _condense(entries: list[dict], key: str = "unp_position") -> dict:
    positions = [e[key] for e in entries]
    return {"count": len(entries), "positions_summary": compact_ranges(positions)}


def build_summary(
    *,
    accession: str,
    uniprot_id: str,
    isoform_note: str,
    domain_start: int,
    domain_end: int,
    needs_trim: bool,
    trim_reason: str,
    chosen_pdb: Optional[dict],
    ranked_candidates: list[dict],
    trim_result: dict,
    avoid: dict,
    hotspots_unp: list[int],
    hotspot_source: str,
    partner_chains: list[str],
    other_interface_contacts_unp: list[int] | None = None,
    numbering_mismatches: dict | None = None,
    verbose: bool = False,
    alphafold_fields: Optional[dict] = None,
) -> dict:
    other_interface_contacts_unp = other_interface_contacts_unp or []
    numbering_mismatches = numbering_mismatches or {}

    if chosen_pdb is not None:
        structure_selection = {
            "chosen_pdb_id": chosen_pdb["pdb_id"],
            "resolution": chosen_pdb["resolution"],
            "ecd_coverage": chosen_pdb["ecd_coverage"],
            "has_known_partner": chosen_pdb["has_partner"],
            "rank_score": chosen_pdb["score"],
            "reason": (
                f"Selected {chosen_pdb['pdb_id']} from {len(ranked_candidates)} candidate "
                f"structures: best combination of extracellular-domain coverage "
                f"({chosen_pdb['ecd_coverage']*100:.0f}%), resolution ({chosen_pdb['resolution']} A)"
                + (", and presence of a bound partner chain for interface-based hotspot detection."
                   if chosen_pdb["has_partner"] else ".")
            ),
            "top_alternatives": ranked_candidates[1:4],
        }
    else:
        n = len(ranked_candidates)
        structure_selection = {
            "chosen_pdb_id": None,
            "resolution": None,
            "ecd_coverage": None,
            "has_known_partner": False,
            "rank_score": None,
            "reason": (
                f"No crystal structure with adequate ECD coverage found "
                f"({n} candidate{'s' if n != 1 else ''} examined); "
                "falling back to AlphaFold prediction."
            ),
            "top_alternatives": ranked_candidates[:4],
        }

    if verbose:
        avoid_section: dict = {
            "glycosylation_sites": avoid["glycosylation"],
            "disulfide_cysteines": avoid["disulfide_cysteines"],
            "other_ptms": avoid["other_ptms"],
            "missing_unresolved": avoid["missing_unresolved"],
        }
        if "low_plddt" in avoid:
            avoid_section["low_plddt"] = avoid["low_plddt"]
        avoid_section["counts"] = {k: len(v) for k, v in avoid.items()}
    else:
        avoid_section = {
            "glycosylation_sites": _condense(avoid["glycosylation"]),
            "disulfide_cysteines": _condense(avoid["disulfide_cysteines"]),
            "other_ptms": _condense(avoid["other_ptms"]),
            "missing_unresolved": _condense(avoid["missing_unresolved"]),
            "counts": {k: len(v) for k, v in avoid.items()},
            "note": "Per-residue detail omitted; request /api/run?verbose=true for the full list.",
        }
        if "low_plddt" in avoid:
            avoid_section["low_plddt"] = _condense(avoid["low_plddt"])

    out = {
        "target": {
            "uniprot_accession": accession,
            "uniprot_id": uniprot_id,
            "isoform": isoform_note,
        },
        "extracellular_domain": {
            "start": domain_start,
            "end": domain_end,
            "trimming_decision": {
                "trimmed": needs_trim,
                "reason": trim_reason,
            },
        },
        "structure_selection": structure_selection,
        "trimming": {
            "auth_chain": trim_result["auth_chain"],
            "residues_in_original_chain": trim_result["n_residues_in_original_chain"],
            "residues_kept": trim_result["n_residues_kept"],
            "residues_trimmed_away": trim_result["n_residues_trimmed_away"],
            "missing_unresolved_label_positions": trim_result["missing_label_seq_in_range"],
        },
        "avoid_residues": avoid_section,
        "hotspots": {
            "source": hotspot_source,
            "partner_chains": partner_chains,
            "candidate_residues_unp": hotspots_unp,
            "candidate_residues_summary": compact_ranges(hotspots_unp),
            "count": len(hotspots_unp),
            "other_interface_contacts_unp": other_interface_contacts_unp if verbose else compact_ranges(other_interface_contacts_unp),
            "other_interface_contacts_count": len(other_interface_contacts_unp),
        },
        "numbering_warnings": [
            {"unp_position": pos, **info} for pos, info in sorted(numbering_mismatches.items())
        ],
    }

    if alphafold_fields:
        out.update(alphafold_fields)

    return out
