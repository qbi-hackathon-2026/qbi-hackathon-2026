"""Download a PDB structure and trim it to the extracellular-domain residue range."""
import requests
import gemmi

from sifts import fetch_unp_segments, unp_range_to_label_seq, parse_poly_seq_scheme


def download_cif(pdb_id: str, dest_path: str) -> str:
    r = requests.get(f"https://files.rcsb.org/download/{pdb_id.upper()}.cif", timeout=30)
    r.raise_for_status()
    with open(dest_path, "wb") as f:
        f.write(r.content)
    return dest_path


def trim_to_domain(
    cif_path: str,
    pdb_id: str,
    accession: str,
    unp_start: int,
    unp_end: int,
) -> dict:
    """Returns trimmed gemmi.Structure plus residue-level bookkeeping for annotation/summary."""
    segments = fetch_unp_segments(pdb_id, accession)
    if not segments:
        raise ValueError(f"No SIFTS mapping found for {pdb_id} / {accession}")

    struct_asym = segments[0]["struct_asym_id"]
    auth_chain = segments[0]["chain_id"]

    label_ranges = unp_range_to_label_seq(segments, unp_start, unp_end, struct_asym)
    if not label_ranges:
        raise ValueError(f"UniProt range {unp_start}-{unp_end} not covered by {pdb_id} chain {struct_asym}")

    scheme = parse_poly_seq_scheme(cif_path, struct_asym)

    keep_auth_nums = set()
    for lo, hi in label_ranges:
        for label_seq in range(lo, hi + 1):
            info = scheme.get(label_seq)
            if info is not None and info["present"]:
                keep_auth_nums.add(info["auth_seq_num"])

    full = gemmi.read_structure(cif_path)
    full.setup_entities()

    trimmed = gemmi.Structure()
    trimmed.name = f"{pdb_id}_ECD_trim"
    model = gemmi.Model("1")
    src_model = full[0]
    src_chain = src_model[auth_chain]

    new_chain = gemmi.Chain(auth_chain)
    n_total, n_kept = 0, 0
    for res in src_chain:
        if res.het_flag == "H" or res.name == "HOH":
            continue
        n_total += 1
        if res.seqid.num in keep_auth_nums:
            new_chain.add_residue(res)
            n_kept += 1
    model.add_chain(new_chain)
    trimmed.add_model(model)
    trimmed.setup_entities()

    return {
        "structure": trimmed,
        "pdb_id": pdb_id,
        "auth_chain": auth_chain,
        "struct_asym": struct_asym,
        "unp_range": (unp_start, unp_end),
        "n_residues_in_original_chain": n_total,
        "n_residues_kept": n_kept,
        "n_residues_trimmed_away": n_total - n_kept,
        "missing_label_seq_in_range": [
            ls for lo, hi in label_ranges for ls in range(lo, hi + 1)
            if scheme.get(ls) and not scheme[ls]["present"]
        ],
    }


def trim_alphafold_to_domain(
    structure: gemmi.Structure,
    chain_id: str,
    domain_start: int,
    domain_end: int,
) -> dict:
    """Trim an AlphaFold structure to ECD bounds.

    AlphaFold auth_seq_id == UniProt position, so no SIFTS mapping is needed —
    we filter directly by residue number within [domain_start, domain_end].
    """
    trimmed = gemmi.Structure()
    trimmed.name = "AF_ECD_trim"
    model = gemmi.Model("1")
    new_chain = gemmi.Chain(chain_id)

    src_chain = structure[0][chain_id]
    n_total = n_kept = 0
    for res in src_chain:
        if res.het_flag == "H" or res.name == "HOH":
            continue
        n_total += 1
        if domain_start <= res.seqid.num <= domain_end:
            new_chain.add_residue(res)
            n_kept += 1

    model.add_chain(new_chain)
    trimmed.add_model(model)
    trimmed.setup_entities()

    return {
        "structure": trimmed,
        "auth_chain": chain_id,
        "struct_asym": chain_id,
        "unp_range": (domain_start, domain_end),
        "n_residues_in_original_chain": n_total,
        "n_residues_kept": n_kept,
        "n_residues_trimmed_away": n_total - n_kept,
        "missing_label_seq_in_range": [],  # pLDDT-unobserved handled separately
    }


def write_structure(structure: gemmi.Structure, out_path: str):
    if out_path.endswith(".pdb"):
        structure.write_pdb(out_path)
    else:
        structure.make_mmcif_document().write_file(out_path)


if __name__ == "__main__":
    cif_path = download_cif("1YY9", "1YY9_full.cif")
    result = trim_to_domain(cif_path, "1YY9", "P00533", 25, 645)
    print({k: v for k, v in result.items() if k != "structure"})
    write_structure(result["structure"], "1YY9_trimmed.pdb")
    print("wrote 1YY9_trimmed.pdb")
