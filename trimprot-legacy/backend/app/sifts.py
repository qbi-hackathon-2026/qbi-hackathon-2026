"""UniProt <-> PDB residue mapping, built on PDBe SIFTS + mmCIF poly_seq_scheme.

Mapping is done via the entity's sequential `label_seq_id` (stable, gap-free
over the full construct) rather than author residue numbers, which can carry
author-defined gaps/insertion codes. Author numbers are recovered per-residue
from `_pdbx_poly_seq_scheme` for writing out the trimmed structure.
"""
import requests
import gemmi

SIFTS_API = "https://www.ebi.ac.uk/pdbe/api/mappings/uniprot_segments"


def fetch_unp_segments(pdb_id: str, accession: str) -> list[dict]:
    """SIFTS segments mapping UniProt residue ranges to label_seq_id ranges, per chain."""
    r = requests.get(f"{SIFTS_API}/{pdb_id.lower()}", timeout=30)
    r.raise_for_status()
    data = r.json()
    entry = data.get(pdb_id.lower(), {}).get("UniProt", {}).get(accession)
    if not entry:
        return []
    segments = []
    for m in entry["mappings"]:
        segments.append({
            "chain_id": m["chain_id"],  # author chain id
            "struct_asym_id": m["struct_asym_id"],  # label asym id
            "unp_start": m["unp_start"],
            "unp_end": m["unp_end"],
            "label_seq_start": m["start"]["residue_number"],
            "label_seq_end": m["end"]["residue_number"],
        })
    return segments


def unp_range_to_label_seq(segments: list[dict], unp_start: int, unp_end: int, struct_asym_id: str | None = None) -> list[tuple[int, int]]:
    """Translate a UniProt residue range into label_seq_id ranges (possibly several segments)."""
    out = []
    for seg in segments:
        if struct_asym_id and seg["struct_asym_id"] != struct_asym_id:
            continue
        lo = max(unp_start, seg["unp_start"])
        hi = min(unp_end, seg["unp_end"])
        if lo > hi:
            continue
        offset = seg["label_seq_start"] - seg["unp_start"]
        out.append((lo + offset, hi + offset))
    return out


def parse_poly_seq_scheme(cif_path: str, struct_asym_id: str) -> dict[int, dict]:
    """label_seq_id -> {auth_seq_num, mon_id, present} for one chain (label asym id)."""
    doc = gemmi.cif.read(cif_path)
    block = doc.sole_block()
    table = block.find(
        "_pdbx_poly_seq_scheme.",
        ["asym_id", "seq_id", "mon_id", "auth_seq_num", "pdb_strand_id"],
    )
    out = {}
    for row in table:
        asym_id, seq_id, mon_id, auth_seq_num, _strand_id = row
        if asym_id != struct_asym_id:
            continue
        seq_id = int(seq_id)
        present = auth_seq_num != "?"
        out[seq_id] = {
            "auth_seq_num": int(auth_seq_num) if present else None,
            "mon_id": mon_id,
            "present": present,
        }
    return out


def unp_position_for_label_seq(segments: list[dict], struct_asym_id: str, label_seq: int) -> int | None:
    for seg in segments:
        if seg["struct_asym_id"] != struct_asym_id:
            continue
        if seg["label_seq_start"] <= label_seq <= seg["label_seq_end"]:
            offset = seg["label_seq_start"] - seg["unp_start"]
            return label_seq - offset
    return None


if __name__ == "__main__":
    import requests as _r
    cif_path = "1YY9.cif"
    open(cif_path, "wb").write(_r.get("https://files.rcsb.org/download/1YY9.cif", timeout=30).content)

    segs = fetch_unp_segments("1YY9", "P00533")
    print("segments:", segs)

    label_ranges = unp_range_to_label_seq(segs, 25, 645)
    print("ECD (unp 25-645) -> label_seq ranges:", label_ranges)

    struct_asym = segs[0]["struct_asym_id"]
    scheme = parse_poly_seq_scheme(cif_path, struct_asym)
    missing = [s for s, info in scheme.items() if not info["present"]]
    print(f"chain {struct_asym}: {len(scheme)} residues total, {len(missing)} missing/unresolved:", missing)
