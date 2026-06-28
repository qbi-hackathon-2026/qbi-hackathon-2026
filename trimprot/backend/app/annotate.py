"""Avoid-residue and hotspot annotation for a trimmed structure."""
import gemmi

CONTACT_CUTOFF = 4.5  # Angstrom, heavy-atom distance for interface contacts

AA3TO1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
    "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
    "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
    "TYR": "Y", "VAL": "V",
}

HYDROPHOBIC_RESNAMES = {"ALA", "VAL", "LEU", "ILE", "MET", "PHE", "TRP", "PRO", "TYR"}


def unp_to_auth_map(scheme: dict, segments: list[dict], struct_asym: str) -> dict[int, int]:
    """UniProt residue position -> auth_seq_num, for residues actually present in the structure."""
    seg = next(s for s in segments if s["struct_asym_id"] == struct_asym)
    offset = seg["label_seq_start"] - seg["unp_start"]
    out = {}
    for label_seq, info in scheme.items():
        if info["present"]:
            out[label_seq - offset] = info["auth_seq_num"]
    return out


def find_numbering_mismatches(
    unp_sequence: str,
    unp_to_auth: dict[int, int],
    structure: gemmi.Structure,
    auth_chain: str,
) -> dict[int, dict]:
    """Spot-check the UniProt<->structure mapping by comparing expected residue type
    (from the UniProt sequence) against the actual residue at the mapped auth position.

    Real systemic numbering bugs show up as a run of consecutive mismatches (a
    frameshift); an isolated single mismatch with correct neighbors on either side is
    a genuine sequence difference between canonical UniProt and the crystallized
    construct (engineered mutation, polymorphism, etc.), not a mapping error - but
    either way it's not safe to treat as the UniProt-annotated residue, so callers
    should exclude these positions from avoid-residue/hotspot reporting.
    """
    chain = structure[0][auth_chain]
    auth_to_resname = {res.seqid.num: res.name for res in chain.get_polymer()}

    mismatches = {}
    for unp_pos, auth_num in unp_to_auth.items():
        expected = unp_sequence[unp_pos - 1] if 1 <= unp_pos <= len(unp_sequence) else None
        actual_resname = auth_to_resname.get(auth_num)
        actual = AA3TO1.get(actual_resname)
        if expected and actual and expected != actual:
            mismatches[unp_pos] = {
                "auth_seq_num": auth_num,
                "expected_residue": expected,
                "structure_residue": actual_resname,
            }
    return mismatches


def avoid_residues(
    features: dict,
    unp_to_auth: dict[int, int],
    unp_start: int,
    unp_end: int,
    missing_unp_positions: set[int],
) -> dict[str, list]:
    """Residues that should not be picked as design hotspots, keyed by reason."""
    avoid = {
        "glycosylation": [],
        "disulfide_cysteines": [],
        "other_ptms": [],
        "missing_unresolved": [],
    }

    for site in features["glycosylation_sites"]:
        pos = site["position"]
        if unp_start <= pos <= unp_end:
            avoid["glycosylation"].append({
                "unp_position": pos,
                "auth_seq_num": unp_to_auth.get(pos),
                "description": site["description"],
            })

    for a, b in features["disulfide_bonds"]:
        for pos in (a, b):
            if unp_start <= pos <= unp_end:
                avoid["disulfide_cysteines"].append({
                    "unp_position": pos,
                    "auth_seq_num": unp_to_auth.get(pos),
                })

    for ptm in features["other_ptms"]:
        pos = ptm["start"]  # modified-residue/lipidation/cross-link sites are single-residue
        if unp_start <= pos <= unp_end:
            avoid["other_ptms"].append({
                "unp_position": pos,
                "auth_seq_num": unp_to_auth.get(pos),
                "ptm_type": ptm["type"],
                "description": ptm["description"],
            })

    for pos in sorted(missing_unp_positions):
        if unp_start <= pos <= unp_end:
            avoid["missing_unresolved"].append({"unp_position": pos})

    return avoid


def find_partner_chains(structure: gemmi.Structure, target_auth_chain: str) -> list[str]:
    model = structure[0]
    return [
        ch.name for ch in model
        if ch.name != target_auth_chain and len(ch.get_polymer()) > 0
    ]


def interface_hotspots(
    full_structure: gemmi.Structure,
    target_auth_chain: str,
    partner_chains: list[str],
    cutoff: float = CONTACT_CUTOFF,
) -> set[int]:
    """Auth_seq_num set of target-chain residues with a heavy-atom contact to any partner chain."""
    model = full_structure[0]
    target_chain = model[target_auth_chain]

    ns = gemmi.NeighborSearch(model, full_structure.cell, cutoff + 0.1).populate()
    contacts = set()
    for res in target_chain.get_polymer():
        for atom in res:
            if atom.element.is_hydrogen:
                continue
            marks = ns.find_atoms(atom.pos, "\0", radius=cutoff)
            for mark in marks:
                cra = mark.to_cra(model)
                if cra.chain.name in partner_chains and cra.residue.het_flag != "H" and cra.residue.name != "HOH":
                    contacts.add(res.seqid.num)
                    break
    return contacts


def prefer_hydrophobic(
    contact_auth_nums: set[int],
    structure: gemmi.Structure,
    target_auth_chain: str,
) -> tuple[set[int], set[int]]:
    """Split an interface-contact set into a hydrophobic-preferred subset and the rest.

    Hydrophobic residues at an existing antibody/partner interface are generally
    preferred binder-design hotspots over charged/polar contact residues. If none of
    the contacts are hydrophobic, the preferred set is empty and callers should fall
    back to the full contact set rather than reporting no hotspots at all.
    """
    chain = structure[0][target_auth_chain]
    resname_by_auth = {res.seqid.num: res.name for res in chain.get_polymer()}

    hydrophobic = {
        n for n in contact_auth_nums
        if resname_by_auth.get(n) in HYDROPHOBIC_RESNAMES
    }
    rest = contact_auth_nums - hydrophobic
    return hydrophobic, rest


def confidence_by_auth(structure: gemmi.Structure, target_auth_chain: str) -> dict[int, float]:
    """Per-residue B-factor (CA atom, or all-atom mean if no CA) for the target chain.

    This is a relative flexibility/disorder proxy for crystal structures, not a
    calibrated confidence score like AlphaFold's pLDDT - lower B-factor means more
    rigidly/consistently resolved across the crystal, not "more correct" in an
    absolute sense. Used by the frontend's confidence-coloring viz mode.
    """
    chain = structure[0][target_auth_chain]
    out = {}
    for res in chain.get_polymer():
        ca = res.find_atom("CA", "\0")
        if ca is not None:
            out[res.seqid.num] = ca.b_iso
        else:
            atoms = list(res)
            if atoms:
                out[res.seqid.num] = sum(a.b_iso for a in atoms) / len(atoms)
    return out


def surface_exposed_hotspots(
    structure: gemmi.Structure,
    target_auth_chain: str,
    avoid_auth_nums: set[int],
    contact_radius: float = 10.0,
    max_contacts: int = 18,
) -> set[int]:
    """Fallback when no bound partner exists: approximate surface exposure via CB
    contact-number (gemmi has no built-in SASA). Fewer CB neighbors within
    `contact_radius` => more solvent-exposed; below `max_contacts` counts as exposed.
    Excludes avoid-listed residues.
    """
    model = structure[0]
    chain = model[target_auth_chain]

    cb_positions = []
    res_for_pos = []
    for res in chain.get_polymer():
        atom = res.find_atom("CB", "\0") or res.find_atom("CA", "\0")
        if atom is not None:
            cb_positions.append(atom.pos)
            res_for_pos.append(res.seqid.num)

    ns = gemmi.NeighborSearch(model, structure.cell, contact_radius + 0.1).populate()
    exposed = set()
    for pos, num in zip(cb_positions, res_for_pos):
        if num in avoid_auth_nums:
            continue
        marks = ns.find_atoms(pos, "\0", radius=contact_radius)
        neighbor_residues = {mark.to_cra(model).residue.seqid.num for mark in marks}
        neighbor_residues.discard(num)
        if len(neighbor_residues) <= max_contacts:
            exposed.add(num)
    return exposed


if __name__ == "__main__":
    import uniprot
    from sifts import fetch_unp_segments, parse_poly_seq_scheme

    cif_path = "1YY9_full.cif"
    full = gemmi.read_structure(cif_path)
    full.setup_entities()

    entry = uniprot.fetch_entry("P00533")
    features = uniprot.extract_features(entry)

    segments = fetch_unp_segments("1YY9", "P00533")
    struct_asym = segments[0]["struct_asym_id"]
    auth_chain = segments[0]["chain_id"]
    scheme = parse_poly_seq_scheme(cif_path, struct_asym)

    u2a = unp_to_auth_map(scheme, segments, struct_asym)
    missing_unp = {
        ls - (segments[0]["label_seq_start"] - segments[0]["unp_start"])
        for ls, info in scheme.items() if not info["present"]
    }

    avoid = avoid_residues(features, u2a, 25, 645, missing_unp)
    print("avoid counts:", {k: len(v) for k, v in avoid.items()})

    partners = find_partner_chains(full, auth_chain)
    print("partner chains:", partners)

    hotspots_auth = interface_hotspots(full, auth_chain, partners)
    auth_to_unp = {v: k for k, v in u2a.items()}
    hotspots_unp = sorted(auth_to_unp.get(a) for a in hotspots_auth if a in auth_to_unp)
    print(f"{len(hotspots_unp)} interface hotspot residues (unp numbering):", hotspots_unp)
