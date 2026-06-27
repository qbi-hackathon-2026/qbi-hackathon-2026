"""Assemble the structured design-prep summary described in the TrimProt spec."""


def build_summary(
    *,
    accession: str,
    uniprot_id: str,
    isoform_note: str,
    domain_start: int,
    domain_end: int,
    needs_trim: bool,
    trim_reason: str,
    chosen_pdb: dict,
    ranked_candidates: list[dict],
    trim_result: dict,
    avoid: dict,
    hotspots_unp: list[int],
    hotspot_source: str,
    partner_chains: list[str],
) -> dict:
    return {
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
        "structure_selection": {
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
        },
        "trimming": {
            "auth_chain": trim_result["auth_chain"],
            "residues_in_original_chain": trim_result["n_residues_in_original_chain"],
            "residues_kept": trim_result["n_residues_kept"],
            "residues_trimmed_away": trim_result["n_residues_trimmed_away"],
            "missing_unresolved_label_positions": trim_result["missing_label_seq_in_range"],
        },
        "avoid_residues": {
            "glycosylation_sites": avoid["glycosylation"],
            "disulfide_cysteines": avoid["disulfide_cysteines"],
            "other_ptms": avoid["other_ptms"],
            "missing_unresolved": avoid["missing_unresolved"],
            "counts": {k: len(v) for k, v in avoid.items()},
        },
        "hotspots": {
            "source": hotspot_source,
            "partner_chains": partner_chains,
            "candidate_residues_unp": hotspots_unp,
            "count": len(hotspots_unp),
        },
    }
