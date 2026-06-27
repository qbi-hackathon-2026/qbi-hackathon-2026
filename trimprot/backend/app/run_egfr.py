"""End-to-end TrimProt pipeline.

Orchestrates UniProt feature extraction → PDB structure search/ranking →
structure trimming → avoid-residue/hotspot annotation → summary output.

When no crystal structure with adequate ECD coverage is found the pipeline
falls back to the AlphaFold model (fetch_alphafold.py), skips SIFTS, and
applies pLDDT-based masking in place of missing-residue avoidance.
"""
import json
import os
import shutil
from pathlib import Path

import uniprot
import pdb_search
import sifts
import trim
import annotate
import summary as summary_mod
import fetch_alphafold as af_mod
import gemmi

ACCESSION = "P00533"
OUT_DIR = os.path.join(os.path.dirname(__file__), "..", "output")

# Fall back to AlphaFold when the best-ranked crystal structure covers less
# than this fraction of the extracellular domain.
MIN_ECD_COVERAGE = 0.5


def run(accession: str = ACCESSION, verbose: bool = False) -> dict:
    os.makedirs(OUT_DIR, exist_ok=True)

    # ── Stage 1: UniProt feature extraction ──────────────────────────────────
    entry = uniprot.fetch_entry(accession)
    features = uniprot.extract_features(entry)
    needs_trim, trim_reason = uniprot.needs_trimming(features)
    domain = uniprot.get_extracellular_domain(features)
    if needs_trim and domain is None:
        raise ValueError("Target has a transmembrane region but no annotated extracellular domain")
    domain_start, domain_end = (domain["start"], domain["end"]) if domain else (1, features["length"])

    isoform_note = (
        f"Canonical sequence {entry['uniProtkbId']} ({features['length']} aa); "
        "no isoform selection needed (single reviewed canonical isoform used)."
    )

    # ── Stage 2: PDB structure search and ranking ─────────────────────────────
    entity_ids = pdb_search.find_entity_ids(accession)
    entry_ids = sorted({eid.split("_")[0] for eid in entity_ids})
    details = pdb_search.fetch_entry_details(entry_ids)
    domain_length = domain_end - domain_start + 1
    ranked = pdb_search.rank_structures(details, accession, domain_length)

    use_alphafold = (not ranked) or (ranked[0]["ecd_coverage"] < MIN_ECD_COVERAGE)

    # ── Stage 3a: AlphaFold fallback path ────────────────────────────────────
    if use_alphafold:
        af_cache = Path(OUT_DIR) / "alphafold_cache"
        af_result = af_mod.get_alphafold_structure(
            accession,
            ecd_ranges=[(domain_start, domain_end)],
            cache_dir=af_cache,
        )
        chain_id = af_result.chain_id  # always "A"

        full_pdb_path = os.path.join(OUT_DIR, f"AF-{accession}_full.pdb")
        af_result.structure.write_pdb(full_pdb_path)

        # Copy the cached CIF into output so /api/files/ can serve it
        cif_src = af_cache / f"AF-{accession}-F1-model_v{af_result.version}.cif"
        cif_dest = os.path.join(OUT_DIR, f"AF-{accession}_full.cif")
        shutil.copy2(str(cif_src), cif_dest)

        trim_result = trim.trim_alphafold_to_domain(
            af_result.structure, chain_id, domain_start, domain_end
        )
        trimmed_path = os.path.join(OUT_DIR, f"AF-{accession}_ECD_trimmed.pdb")
        trim.write_structure(trim_result["structure"], trimmed_path)

        # AF auth_seq_id == UniProt position — identity mapping, no SIFTS needed
        u2a = {i: i for i in range(domain_start, domain_end + 1)}

        # pLDDT ≤ 50 residues are the AF equivalent of missing electron density
        plddt_unobs = {
            e.auth_seq_id for e in af_result.avoid_contributions
            if "unobserved" in e.reason
        }

        avoid = annotate.avoid_residues(features, u2a, domain_start, domain_end, plddt_unobs)

        # Low-pLDDT (50 < pLDDT < 70) as a distinct avoid category
        avoid["low_plddt"] = [
            {
                "unp_position": e.auth_seq_id,
                "auth_seq_num": e.auth_seq_id,
                "plddt": round(e.plddt, 1),
                "reason": e.reason,
            }
            for e in af_result.avoid_contributions
            if "low_plddt" in e.reason
        ]

        avoid_auth_nums = {
            e["auth_seq_num"]
            for cat in ("glycosylation", "disulfide_cysteines", "other_ptms", "low_plddt")
            for e in avoid[cat]
            if e["auth_seq_num"] is not None
        }
        avoid_auth_nums |= plddt_unobs

        auth_to_unp = {v: k for k, v in u2a.items()}

        # No partner chains in an AF monomer — always use surface-exposure fallback
        hotspots_auth = annotate.surface_exposed_hotspots(
            trim_result["structure"], chain_id, avoid_auth_nums
        )
        hotspots_auth -= avoid_auth_nums
        hotspot_source = (
            f"inferred surface-exposed residues (AlphaFold model v{af_result.version}; "
            f"pLDDT masking applied, mean ECD pLDDT={af_result.ecd_mean_plddt:.1f})"
        )
        hotspots_unp = sorted(auth_to_unp[a] for a in hotspots_auth if a in auth_to_unp)

        isoform_note = (
            f"Canonical sequence {entry['uniProtkbId']} ({features['length']} aa); "
            "AlphaFold model used (no crystal structure with adequate ECD coverage found)."
        )

        result_summary = summary_mod.build_summary(
            accession=accession,
            uniprot_id=entry["uniProtkbId"],
            isoform_note=isoform_note,
            domain_start=domain_start,
            domain_end=domain_end,
            needs_trim=needs_trim,
            trim_reason=trim_reason,
            chosen_pdb=None,
            ranked_candidates=ranked,
            trim_result=trim_result,
            avoid=avoid,
            hotspots_unp=hotspots_unp,
            hotspot_source=hotspot_source,
            partner_chains=[],
            alphafold_fields=af_result.summary_fields(),
        )

        avoid_glyco_auth     = sorted({e["auth_seq_num"] for e in avoid["glycosylation"]      if e["auth_seq_num"] is not None})
        avoid_disulfide_auth = sorted({e["auth_seq_num"] for e in avoid["disulfide_cysteines"] if e["auth_seq_num"] is not None})
        avoid_ptm_auth       = sorted({e["auth_seq_num"] for e in avoid["other_ptms"]          if e["auth_seq_num"] is not None})
        avoid_low_plddt_auth = sorted({e["auth_seq_num"] for e in avoid["low_plddt"]           if e["auth_seq_num"] is not None})

        result_summary["viewer"] = {
            "pdb_id": f"AF-{accession}",
            "auth_chain": chain_id,
            "partner_chains": [],
            "original_file": f"AF-{accession}_full.pdb",
            "original_cif_file": f"AF-{accession}_full.cif",
            "trimmed_file": f"AF-{accession}_ECD_trimmed.pdb",
            "hotspot_auth_residues": sorted(hotspots_auth),
            "avoid_glycosylation_auth_residues": avoid_glyco_auth,
            "avoid_disulfide_auth_residues": avoid_disulfide_auth,
            "avoid_other_ptm_auth_residues": avoid_ptm_auth,
            "avoid_low_plddt_auth_residues": avoid_low_plddt_auth,
            "missing_auth_residues_note": "pLDDT ≤ 50 residues treated as unobserved",
        }

        summary_path = os.path.join(OUT_DIR, f"AF-{accession}_summary.json")
        with open(summary_path, "w") as f:
            json.dump(result_summary, f, indent=2)

        print(f"Wrote trimmed structure (AlphaFold): {trimmed_path}")
        print(f"Wrote summary: {summary_path}")
        return result_summary

    # ── Stage 3b: Crystal structure path ─────────────────────────────────────
    chosen = ranked[0]
    pdb_id = chosen["pdb_id"]

    cif_path = os.path.join(OUT_DIR, f"{pdb_id}_full.cif")
    trim.download_cif(pdb_id, cif_path)

    trim_result = trim.trim_to_domain(cif_path, pdb_id, accession, domain_start, domain_end)
    trimmed_path = os.path.join(OUT_DIR, f"{pdb_id}_ECD_trimmed.pdb")
    trim.write_structure(trim_result["structure"], trimmed_path)

    full_structure = gemmi.read_structure(cif_path)
    full_structure.setup_entities()

    full_pdb_path = os.path.join(OUT_DIR, f"{pdb_id}_full.pdb")
    full_structure.write_pdb(full_pdb_path)

    segments = sifts.fetch_unp_segments(pdb_id, accession)
    struct_asym = trim_result["struct_asym"]
    auth_chain = trim_result["auth_chain"]
    scheme = sifts.parse_poly_seq_scheme(cif_path, struct_asym)

    u2a = annotate.unp_to_auth_map(scheme, segments, struct_asym)
    offset = segments[0]["label_seq_start"] - segments[0]["unp_start"]
    missing_unp = {ls - offset for ls in trim_result["missing_label_seq_in_range"]}

    numbering_mismatches = annotate.find_numbering_mismatches(features["sequence"], u2a, full_structure, auth_chain)
    if numbering_mismatches:
        # Don't trust UniProt-annotated positions where the structure's actual residue
        # doesn't match the canonical UniProt sequence (engineered mutation, construct
        # variant, etc.) - drop them from the mapping entirely rather than silently
        # mislabeling a residue.
        u2a = {pos: auth for pos, auth in u2a.items() if pos not in numbering_mismatches}

    avoid = annotate.avoid_residues(features, u2a, domain_start, domain_end, missing_unp)
    avoid_auth_nums = {
        e["auth_seq_num"] for cat in ("glycosylation", "disulfide_cysteines", "other_ptms") for e in avoid[cat]
        if e["auth_seq_num"] is not None
    }

    partner_chains = annotate.find_partner_chains(full_structure, auth_chain)
    auth_to_unp = {v: k for k, v in u2a.items()}

    if partner_chains:
        all_contacts_auth = annotate.interface_hotspots(full_structure, auth_chain, partner_chains)
        all_contacts_auth -= avoid_auth_nums
        hydrophobic_auth, other_contacts_auth = annotate.prefer_hydrophobic(all_contacts_auth, full_structure, auth_chain)
        # Prefer the hydrophobic subset of the real interface as hotspots; only fall
        # back to the full contact set if none of the contacts are hydrophobic.
        hotspots_auth = hydrophobic_auth if hydrophobic_auth else all_contacts_auth
        other_interface_contacts_unp = sorted(auth_to_unp[a] for a in other_contacts_auth if a in auth_to_unp)
        hotspot_source = (
            f"known partner interface (chains {partner_chains}), preferring hydrophobic contact residues"
            if hydrophobic_auth else
            f"known partner interface (chains {partner_chains}); no hydrophobic contacts found, showing all contacts"
        )
    else:
        hotspots_auth = annotate.surface_exposed_hotspots(full_structure, auth_chain, avoid_auth_nums)
        hotspot_source = "inferred surface-exposed residues (no bound partner in chosen structure)"
        other_interface_contacts_unp = []

    hotspots_auth -= avoid_auth_nums
    hotspots_unp = sorted(auth_to_unp[a] for a in hotspots_auth if a in auth_to_unp)

    result_summary = summary_mod.build_summary(
        accession=accession,
        uniprot_id=entry["uniProtkbId"],
        isoform_note=isoform_note,
        domain_start=domain_start,
        domain_end=domain_end,
        needs_trim=needs_trim,
        trim_reason=trim_reason,
        chosen_pdb=chosen,
        ranked_candidates=ranked,
        trim_result=trim_result,
        avoid=avoid,
        hotspots_unp=hotspots_unp,
        hotspot_source=hotspot_source,
        partner_chains=partner_chains,
        other_interface_contacts_unp=other_interface_contacts_unp,
        numbering_mismatches=numbering_mismatches,
        verbose=verbose,
    )

    # Surface protein display fields (name, gene symbols, organism, length) on the
    # summary so the frontend title block renders on every path (search-select,
    # reload, direct) without depending on the transient search-result object.
    result_summary["target"].update(uniprot._candidate_from_entry(entry))

    avoid_glyco_auth = sorted({e["auth_seq_num"] for e in avoid["glycosylation"] if e["auth_seq_num"] is not None})
    avoid_disulfide_auth = sorted({e["auth_seq_num"] for e in avoid["disulfide_cysteines"] if e["auth_seq_num"] is not None})
    avoid_ptm_auth = sorted({e["auth_seq_num"] for e in avoid["other_ptms"] if e["auth_seq_num"] is not None})
    other_contacts_auth = sorted(u2a[p] for p in other_interface_contacts_unp if p in u2a)
    result_summary["viewer"] = {
        "pdb_id": pdb_id,
        "auth_chain": auth_chain,
        "partner_chains": partner_chains,
        "original_file": f"{pdb_id}_full.pdb",
        "original_cif_file": f"{pdb_id}_full.cif",
        "trimmed_file": f"{pdb_id}_ECD_trimmed.pdb",
        "hotspot_auth_residues": sorted(hotspots_auth),
        "other_interface_contact_auth_residues": other_contacts_auth,
        "avoid_glycosylation_auth_residues": avoid_glyco_auth,
        "avoid_disulfide_auth_residues": avoid_disulfide_auth,
        "avoid_other_ptm_auth_residues": avoid_ptm_auth,
        "missing_auth_residues_note": "missing residues have no atoms and cannot be highlighted",
    }

    summary_path = os.path.join(OUT_DIR, f"{pdb_id}_summary.json")
    with open(summary_path, "w") as f:
        json.dump(result_summary, f, indent=2)

    print(f"Wrote trimmed structure: {trimmed_path}")
    print(f"Wrote summary: {summary_path}")
    return result_summary


if __name__ == "__main__":
    result = run()
    print(json.dumps(result, indent=2)[:2000])
