"""gemmi structure loading: assembly expansion, chain classification, glycans.

Loads the deposited mmCIF and (by default) expands the biological assembly with
`transform_to_assembly`, which preserves author numbering and keeps the original
copy's chain id stable (copies get a numeric suffix). Chain classification is done
by sequence identity against the target chain so we never need fragile entity
metadata: identical-sequence chains are homo-oligomer copies; differing chains are
true partners (antibody/receptor).
"""
from __future__ import annotations

from dataclasses import dataclass

import gemmi

from .cache import download_cif

# PDB chemical-component ids for common glycan / sugar residues seen as HETATM.
SUGAR_CCDS = {
    "NAG", "NDG", "BMA", "MAN", "BGC", "GLC", "GAL", "GLA", "FUC", "FUL",
    "A2G", "NGA", "SIA", "NAN", "NEU", "GLB", "XYP", "XYS", "RAM", "RIB",
    "MMA", "BNX", "BNG", "GCS", "IDS", "SGN", "GCU", "ADA", "KDN",
}


@dataclass
class LoadedStructure:
    analysis: gemmi.Structure   # asymmetric unit (model 1): SIFTS/interface/avoid
    display: gemmi.Structure    # assembly-expanded (or AU): trim + viewer
    pdb_id: str
    assembly: str               # 'bioassembly' or 'protomer'
    assembly_applied: bool
    warnings: list[str]


def _first_model_only(st: gemmi.Structure) -> None:
    while len(st) > 1:
        del st[1]


def load_structure(pdb_id: str, assembly: str = "bioassembly") -> LoadedStructure:
    """Read a structure.

    Analysis always uses the deposited asymmetric unit so chain ids line up with
    SIFTS. For ``bioassembly`` the display/trim copy is expanded with
    ``transform_to_assembly`` (auth numbering preserved; the first copy keeps its
    original chain id, additional copies get a numeric suffix).
    """
    path = download_cif(pdb_id)
    analysis = gemmi.read_structure(str(path))
    analysis.setup_entities()
    _first_model_only(analysis)

    warnings: list[str] = []
    applied = False
    display = analysis
    if assembly == "bioassembly":
        names = [a.name for a in analysis.assemblies]
        if names:
            display = gemmi.read_structure(str(path))
            display.setup_entities()
            _first_model_only(display)
            try:
                display.transform_to_assembly(
                    names[0], gemmi.HowToNameCopiedChain.AddNumber)
                applied = True
            except Exception as exc:  # pragma: no cover - defensive
                warnings.append(f"assembly expansion failed ({exc}); using AU for display")
                display = analysis
        else:
            warnings.append("no assembly records; using deposited asymmetric unit")

    return LoadedStructure(analysis=analysis, display=display, pdb_id=pdb_id.lower(),
                           assembly=assembly, assembly_applied=applied,
                           warnings=warnings)


def polymer_sequence(chain: gemmi.Chain) -> str:
    poly = chain.get_polymer()
    if poly.length() == 0:
        return ""
    return poly.make_one_letter_sequence()


def _identity(a: str, b: str, k: int = 4) -> float:
    """Alignment-free sequence similarity via shared k-mers (Jaccard-like).

    Robust to register shifts and ragged termini: two copies of the same protein
    that differ only by a missing N-terminal residue (e.g. BCMA chains F vs K in
    4ZFO) still score ~1.0, whereas a positional comparison would collapse to ~0.
    """
    if not a or not b or len(a) < k or len(b) < k:
        # too short for k-mers: fall back to exact / substring containment
        return 1.0 if a and b and (a in b or b in a) else 0.0
    ka = {a[i:i + k] for i in range(len(a) - k + 1)}
    kb = {b[i:i + k] for i in range(len(b) - k + 1)}
    if not ka or not kb:
        return 0.0
    return len(ka & kb) / min(len(ka), len(kb))


@dataclass
class ChainClasses:
    target_chain: str
    same_protein: list[str]   # homo-oligomer copies (incl. target)
    partners: list[str]       # different-protein polymer chains


def classify_chains(structure: gemmi.Structure, target_chain: str,
                    identity_threshold: float = 0.9,
                    target_seq: str | None = None) -> ChainClasses:
    """Split polymer chains into same-protein copies vs partners by sequence.

    By default the target's sequence is looked up by chain name within
    ``structure`` itself (matching the deposited asymmetric unit, where the
    selected chain id is always present). Pass ``target_seq`` explicitly when
    classifying a *different* structure than the one ``target_chain`` was chosen
    from - e.g. an assembly-expanded copy, where biological-assembly generation
    can legitimately rename or drop the originally-selected chain id entirely
    (it may not exist under that name, or any name sharing its prefix). Without
    this, the name-based lookup silently fails, target_seq stays empty, identity
    against it is always 0.0, and every chain - including real homo-oligomer
    copies of the target - gets misclassified as a "partner" (or no chain is
    classified as the target's copy at all).
    """
    model = structure[0]
    seqs: dict[str, str] = {}
    for ch in model:
        s = polymer_sequence(ch)
        if s:
            seqs[ch.name] = s

    if target_seq:
        tgt_seq = target_seq
        same = [n for n, s in seqs.items() if _identity(tgt_seq, s) >= identity_threshold]
        partners = [n for n in seqs if n not in same]
        tgt_name = sorted(same)[0] if same else target_chain
        return ChainClasses(target_chain=tgt_name,
                            same_protein=sorted(same),
                            partners=sorted(partners))

    # The target chain name may carry an assembly suffix; match by prefix too.
    tgt_name = target_chain
    if tgt_name not in seqs:
        cands = [n for n in seqs if n == target_chain or n.startswith(target_chain)]
        if cands:
            tgt_name = sorted(cands)[0]
    tgt_seq = seqs.get(tgt_name, "")

    same, partners = [], []
    for name, s in seqs.items():
        if name == tgt_name:
            same.append(name)
            continue
        if tgt_seq and _identity(tgt_seq, s) >= identity_threshold:
            same.append(name)
        else:
            partners.append(name)
    return ChainClasses(target_chain=tgt_name,
                        same_protein=sorted(same),
                        partners=sorted(partners))


def sugar_atoms(structure: gemmi.Structure):
    """Heavy atoms of glycan/sugar HETATM residues.

    Returns (positions, sugar_keys) where positions is a list of gemmi.Position
    and sugar_keys is the set of (chain_name, residue_number) the sugars occupy
    (used to avoid flagging sugar-on-sugar contacts).
    """
    positions: list[gemmi.Position] = []
    sugar_keys: set[tuple[str, int]] = set()
    model = structure[0]
    for chain in model:
        for res in chain:
            if res.name in SUGAR_CCDS:
                sugar_keys.add((chain.name, res.seqid.num))
                for atom in res:
                    if not atom.is_hydrogen():
                        positions.append(atom.pos)
    return positions, sugar_keys
