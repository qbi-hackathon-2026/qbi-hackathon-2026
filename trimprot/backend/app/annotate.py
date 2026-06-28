"""Avoid-residue and hotspot annotation for a trimmed structure."""
import math

import gemmi
import numpy as np

CONTACT_CUTOFF = 4.5  # Angstrom, heavy-atom distance for interface contacts

# Bondi van-der-Waals radii (A) for the elements we see in protein heavy atoms.
VDW_RADII = {"C": 1.70, "N": 1.55, "O": 1.52, "S": 1.80, "P": 1.80}
PROBE_RADIUS = 1.4  # water probe

# Tien et al. (2013) theoretical per-residue maximum SASA (A^2), for normalizing
# absolute SASA into a 0..1 relative exposure that's comparable across residue types.
MAX_ASA = {
    "ALA": 129.0, "ARG": 274.0, "ASN": 195.0, "ASP": 193.0, "CYS": 167.0,
    "GLU": 223.0, "GLN": 225.0, "GLY": 104.0, "HIS": 224.0, "ILE": 197.0,
    "LEU": 201.0, "LYS": 236.0, "MET": 224.0, "PHE": 240.0, "PRO": 159.0,
    "SER": 155.0, "THR": 172.0, "TRP": 285.0, "TYR": 263.0, "VAL": 174.0,
}

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


def _fibonacci_sphere(n: int) -> np.ndarray:
    """`n` roughly-uniform points on the unit sphere (golden-spiral construction).
    Used as the test directions for Shrake-Rupley SASA sampling.
    """
    i = np.arange(n) + 0.5
    phi = np.arccos(1.0 - 2.0 * i / n)        # polar angle, uniform in cos
    theta = math.pi * (1.0 + 5.0 ** 0.5) * i  # golden-angle azimuth
    return np.column_stack([
        np.cos(theta) * np.sin(phi),
        np.sin(theta) * np.sin(phi),
        np.cos(phi),
    ])


_SPHERE_POINTS = _fibonacci_sphere(96)


def residue_rsasa(
    structure: gemmi.Structure,
    target_auth_chain: str,
    n_points: int = 96,
) -> dict[int, float]:
    """Per-residue relative solvent-accessible surface area (rSASA, ~0..1) for the
    target chain, via a Shrake-Rupley rolling-probe calculation.

    This is real SASA (each heavy atom's probe-inflated sphere is sampled, and a
    point counts as accessible if no other target-chain atom occludes it), not the
    CB-contact-number proxy `surface_exposed_hotspots` uses. SASA is computed for
    the isolated target chain (partner/other chains ignored) because binder design
    cares about the monomer surface that's available to a new binder. Absolute
    atom SASA is summed per residue and divided by the residue's theoretical max
    (Tien et al.) so 0 = fully buried, ~1 = fully exposed.
    """
    model = structure[0]
    chain = model[target_auth_chain]

    coords, radii, res_of_atom, residues = [], [], [], []
    for res_idx, res in enumerate(chain.get_polymer()):
        residues.append(res)
        for atom in res:
            if atom.element.is_hydrogen:
                continue
            coords.append([atom.pos.x, atom.pos.y, atom.pos.z])
            radii.append(VDW_RADII.get(atom.element.name, 1.70) + PROBE_RADIUS)
            res_of_atom.append(res_idx)

    if not coords:
        return {}

    coords = np.asarray(coords)
    radii = np.asarray(radii)
    res_of_atom = np.asarray(res_of_atom)
    sphere = _SPHERE_POINTS if n_points == 96 else _fibonacci_sphere(n_points)
    max_r = float(radii.max())

    ns = gemmi.NeighborSearch(model, structure.cell, 2.0 * max_r).populate()

    atom_sasa = np.zeros(len(coords))
    for i in range(len(coords)):
        test_pts = coords[i] + radii[i] * sphere  # (P, 3) points on atom i's probe sphere
        accessible = np.ones(len(sphere), dtype=bool)
        marks = ns.find_atoms(gemmi.Position(*coords[i]), "\0", radius=radii[i] + max_r)
        for mark in marks:
            cra = mark.to_cra(model)
            if cra.chain.name != target_auth_chain or cra.atom.element.is_hydrogen:
                continue
            jp = np.array([cra.atom.pos.x, cra.atom.pos.y, cra.atom.pos.z])
            rj = VDW_RADII.get(cra.atom.element.name, 1.70) + PROBE_RADIUS
            # Self never occludes: its own points sit at distance == radii[i] == rj,
            # and the strict `<` below leaves those accessible.
            d2 = ((test_pts - jp) ** 2).sum(axis=1)
            accessible &= d2 >= rj * rj
        atom_sasa[i] = accessible.mean() * 4.0 * math.pi * radii[i] ** 2

    res_sasa = np.zeros(len(residues))
    np.add.at(res_sasa, res_of_atom, atom_sasa)

    out = {}
    for res_idx, res in enumerate(residues):
        max_asa = MAX_ASA.get(res.name)
        if max_asa:
            out[res.seqid.num] = float(res_sasa[res_idx] / max_asa)
    return out


def candidate_surface_patches(
    structure: gemmi.Structure,
    target_auth_chain: str,
    rsasa: dict[int, float],
    avoid_auth_nums: set[int],
    exposed_threshold: float = 0.25,
    patch_radius: float = 10.0,
    min_size: int = 4,
) -> list[list[int]]:
    """Generate candidate surface patches by the sliding-neighborhood method.

    Each solvent-exposed, non-avoid residue seeds one patch = the exposed residues
    whose CB lies within `patch_radius` of the seed's CB. Unlike connected-component
    clustering (which merges the whole contiguous surface sheet into one giant
    blob), this yields compact, localized, fixed-scale patches the size of a real
    epitope; the overlapping candidates are deduplicated and later thinned by
    non-maximum suppression in `surface_patch_hotspots`. `exposed_threshold` follows
    the common 25% relative-accessibility cut for "surface-exposed".
    """
    chain = structure[0][target_auth_chain]
    nums, points = [], []
    for res in chain.get_polymer():
        num = res.seqid.num
        if num in avoid_auth_nums or rsasa.get(num, 0.0) < exposed_threshold:
            continue
        atom = res.find_atom("CB", "\0") or res.find_atom("CA", "\0")
        if atom is not None:
            nums.append(num)
            points.append([atom.pos.x, atom.pos.y, atom.pos.z])

    if not nums:
        return []

    pts = np.asarray(points)
    radius2 = patch_radius * patch_radius

    seen: set[tuple[int, ...]] = set()
    patches: list[list[int]] = []
    for i in range(len(nums)):
        d2 = ((pts - pts[i]) ** 2).sum(axis=1)
        members = tuple(sorted(nums[j] for j in np.nonzero(d2 <= radius2)[0]))
        if len(members) >= min_size and members not in seen:
            seen.add(members)
            patches.append(list(members))
    return patches


def score_patches(
    patches: list[list[int]],
    structure: gemmi.Structure,
    target_auth_chain: str,
    rsasa: dict[int, float],
) -> list[dict]:
    """Rank surface patches for binder-design developability.

    Each patch is scored 0..1 from three normalized terms a de-novo binder cares
    about: hydrophobic fraction (hydrophobic patches make better anchor sites),
    size (a larger contiguous epitope gives a binder more to grip), and mean
    exposure. Returns patch dicts sorted best-first.
    """
    chain = structure[0][target_auth_chain]
    resname_by_auth = {res.seqid.num: res.name for res in chain.get_polymer()}

    scored = []
    for patch in patches:
        size = len(patch)
        hydrophobic = sum(1 for n in patch if resname_by_auth.get(n) in HYDROPHOBIC_RESNAMES)
        hydro_frac = hydrophobic / size
        mean_rsasa = sum(rsasa.get(n, 0.0) for n in patch) / size

        size_score = min(size / 12.0, 1.0)   # ~12 exposed residues ≈ a full epitope
        exposure_score = min(mean_rsasa, 1.0)
        score = 0.40 * hydro_frac + 0.35 * size_score + 0.25 * exposure_score

        scored.append({
            "residues_auth": patch,
            "size": size,
            "hydrophobic_fraction": round(hydro_frac, 2),
            "mean_rsasa": round(mean_rsasa, 2),
            "score": round(score, 3),
        })

    scored.sort(key=lambda p: p["score"], reverse=True)
    return scored


def _select_non_overlapping(ranked: list[dict], max_patches: int, overlap_frac: float = 0.5) -> list[dict]:
    """Greedy non-maximum suppression over score-sorted patches: accept the best
    patch, then skip any later patch that shares more than `overlap_frac` of its
    residues with one already accepted, so the returned set covers distinct regions
    of the surface instead of `max_patches` near-duplicate views of one hotspot.
    """
    selected: list[dict] = []
    for patch in ranked:
        members = set(patch["residues_auth"])
        if any(len(members & set(s["residues_auth"])) / len(members) > overlap_frac for s in selected):
            continue
        selected.append(patch)
        if len(selected) >= max_patches:
            break
    return selected


def surface_patch_hotspots(
    structure: gemmi.Structure,
    target_auth_chain: str,
    avoid_auth_nums: set[int],
    max_patches: int = 3,
) -> tuple[set[int], list[dict]]:
    """SASA-based replacement for `surface_exposed_hotspots` when no partner exists.

    Computes real per-residue SASA, generates sliding-neighborhood surface patches,
    ranks them by developability, applies non-maximum suppression to keep the top
    distinct patches, and returns (union of their residues, ranked patch dicts). The
    residue set keeps the existing hotspot-as-residue-set interface working; the
    patch list is the richer output for the summary/UI.
    """
    rsasa = residue_rsasa(structure, target_auth_chain)
    candidates = candidate_surface_patches(structure, target_auth_chain, rsasa, avoid_auth_nums)
    ranked = score_patches(candidates, structure, target_auth_chain, rsasa)
    top = _select_non_overlapping(ranked, max_patches)
    residues = {n for p in top for n in p["residues_auth"]}
    return residues, top


def surface_exposed_hotspots(
    structure: gemmi.Structure,
    target_auth_chain: str,
    avoid_auth_nums: set[int],
    contact_radius: float = 10.0,
    max_contacts: int = 18,
) -> set[int]:
    """Legacy fallback: approximate surface exposure via CB contact-number. Superseded
    by `surface_patch_hotspots` (real SASA + patch clustering); kept for reference.
    Fewer CB neighbors within `contact_radius` => more solvent-exposed; below
    `max_contacts` counts as exposed. Excludes avoid-listed residues.
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
