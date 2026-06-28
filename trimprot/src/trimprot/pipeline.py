"""End-to-end orchestration. Pure/deterministic; no LLM calls.

Wires the standalone pipeline functions into one PipelineResult that the CLI, the
MCP server, and the test suite all consume through the same code path.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional

import gemmi

from .avoid import AvoidSet, build_avoid_set
from .glyco import predict_glycosylation
from .hotspots import Hotspot, Removal, filter_hotspots, select_patch
from .interface import InterfaceResidue, detect_interface
from .sifts import SiftsMapping, get_sifts_mapping
from .structio import ChainClasses, LoadedStructure, classify_chains, load_structure
from .structures import StructureChoice, pfam_domains, search_structures
from .topology import (
    MembraneProximal,
    Range,
    get_extracellular_ranges,
    get_membrane_proximal_terminus,
    get_transmem_ranges,
)
from .trim import TrimResult, trim_structure
from .uniprot import UniProtRecord, resolve_uniprot

# Above this many homo-oligomer copies the assembly is almost certainly crystal
# packing rather than biology; trim to the protomer instead.
MAX_OLIGOMER_CHAINS = 12


@dataclass
class PipelineResult:
    record: UniProtRecord
    ecd_ranges: list[Range]
    transmem_ranges: list[Range]
    membrane_proximal: Optional[MembraneProximal]
    glyco_sites: set[int]
    choice: StructureChoice
    loaded: LoadedStructure
    chains: ChainClasses
    display_chains: ChainClasses
    mapping: SiftsMapping
    ecd_auth_ranges: list[Range]
    membrane_buffer_auth: Optional[Range]
    avoid: AvoidSet
    interface: list[InterfaceResidue]
    interface_chains: list[str]
    hotspots: list[Hotspot]
    patch: list[Hotspot]
    domains: list[dict]
    removals: list[Removal]
    trim: TrimResult
    assembly: str
    membrane_buffer: int
    warnings: list[str] = field(default_factory=list)


def _ecd_to_auth(ecd_ranges: list[Range], mapping: SiftsMapping,
                 chain: str) -> list[Range]:
    """Map each UniProt ECD range to a PDB-auth range using observed/mapped residues."""
    u2p = mapping.unp_to_pdb(chain)
    out = []
    for r in ecd_ranges:
        auth_nums = [u2p[u][0] for u in r if u in u2p]
        if auth_nums:
            out.append(Range(min(auth_nums), max(auth_nums)))
    return out


def _range_to_auth(r: Optional[Range], mapping: SiftsMapping,
                   chain: str) -> Optional[Range]:
    if r is None:
        return None
    u2p = mapping.unp_to_pdb(chain)
    auth_nums = [u2p[u][0] for u in r if u in u2p]
    if not auth_nums:
        return None
    return Range(min(auth_nums), max(auth_nums))


def run_pipeline(name: Optional[str] = None, *, accession: Optional[str] = None,
                 organism_id: int = 9606, prefer_antibody: bool = False,
                 assembly: str = "bioassembly", membrane_buffer: int = 12,
                 interface_cutoff: float = 5.0, patch_radius: float = 11.0,
                 patch_size: int = 8, no_patch: bool = False,
                 min_ecd_coverage: float = 0.40) -> PipelineResult:
    warnings: list[str] = []

    # 1. UniProt + features
    rec = resolve_uniprot(name, accession=accession, organism_id=organism_id)

    # 2. Topology
    ecd_ranges = get_extracellular_ranges(rec.features)
    tm_ranges = get_transmem_ranges(rec.features)
    mp = get_membrane_proximal_terminus(ecd_ranges, tm_ranges, buffer=membrane_buffer)
    if mp is not None and not mp.agrees:
        warnings.append(f"membrane-proximal terminus {mp.terminus} disagrees with "
                        f"{mp.topo_type} topology expectation")

    # 3. Glycosylation
    glyco_sites = predict_glycosylation(rec.sequence, rec.of("CARBOHYD"))

    # 4. Structure search (strict priority ladder)
    choice = search_structures(rec.accession, ecd_ranges,
                               prefer_antibody=prefer_antibody,
                               min_ecd_coverage=min_ecd_coverage)

    if choice.chosen.predicted:
        warnings.append("no experimental structure; predicted model used")

    # 5. Load (assembly-aware). Analysis uses the AU; display/trim use the assembly.
    loaded = load_structure(choice.chosen.pdb_id, assembly=assembly)
    warnings.extend(loaded.warnings)
    structure = loaded.analysis

    # 6. Classify chains (target vs homo-copies vs partners) from AU coordinates
    chains = classify_chains(structure, choice.target_chain)
    target_chain = chains.target_chain
    au_chain_names = {c.name for c in structure[0]}

    # Interface chains: prefer the designated antibody chains (so we measure the
    # epitope, not crystal contacts or assembly-neighbour copies of the target).
    # Fall back to genuine non-antibody partners, then to the apo path. Same-
    # protein (homo-oligomer) copies are never counted.
    antibody = sorted(set(choice.antibody_chains) & set(chains.partners)
                      & au_chain_names)
    if antibody:
        interface_chains = antibody
        apo_fallback = False
    elif chains.partners and not choice.apo_fallback:
        interface_chains = chains.partners
        apo_fallback = False
    else:
        interface_chains = []
        apo_fallback = True
    partner_chains = interface_chains

    # 7. SIFTS mapping (auth numbering, observed flags, xref cross-check)
    mapping = get_sifts_mapping(choice.chosen.pdb_id, rec.accession,
                                structure=structure, target_chain=target_chain)
    if mapping.xref_mismatches:
        warnings.append(f"{len(mapping.xref_mismatches)} SIFTS/xref mismatches "
                        f"(see summary)")

    # 8. Map ECD + buffer to auth numbering
    ecd_auth = _ecd_to_auth(ecd_ranges, mapping, target_chain)
    membrane_buffer_auth = _range_to_auth(mp.buffer if mp else None,
                                          mapping, target_chain)

    # 9. Avoid set
    disulfide_cys = [f.start for f in rec.of("DISULFID")]
    mod_res = [f.start for f in rec.of("MOD_RES")]
    avoid = build_avoid_set(
        mapping, structure, target_chain,
        glyco_sites=glyco_sites,
        disulfide_cys=disulfide_cys,
        mod_res=mod_res,
        tm_ranges=tm_ranges,
        membrane_buffer=(mp.buffer if mp else None),
        ecd_ranges=ecd_ranges,
    )

    # 10. Interface / epitope (contacts to antibody chains, or apo fallback)
    interface = detect_interface(structure, target_chain, interface_chains,
                                 cutoff=interface_cutoff,
                                 apo_fallback=apo_fallback,
                                 ecd_auth_ranges=ecd_auth)

    # 11. Filter hotspots, then collapse to a single contiguous epitope patch.
    hotspots, removals = filter_hotspots(
        interface, avoid, target_chain, mapping.pdb_to_unp(target_chain))
    if not hotspots:
        why = ("no interface contacts and no exposed structured residues found"
               if not interface else
               f"all {len(interface)} candidates fell in the avoid set")
        warnings.append(f"no hotspots emitted: {why}")
    # --no-patch (transparency mode) emits the full ranked list as the "patch".
    if no_patch:
        patch = list(hotspots)
    else:
        patch = select_patch(hotspots, structure, target_chain,
                             radius=patch_radius, cap=patch_size)

    # Pfam domains (author numbering) for labelling in the viewer
    domains = pfam_domains(choice.chosen.pdb_id, target_chain)

    # 12. Keep chains for trim, classified against the DISPLAY (assembly) coords.
    if loaded.assembly_applied:
        display_chains = classify_chains(loaded.display, choice.target_chain)
    else:
        display_chains = chains
    if assembly == "bioassembly":
        keep_chains = display_chains.same_protein
        # Guard against pathological crystal assemblies (dozens of copies): keep
        # the target protomer only, with a warning, so outputs stay design-usable.
        if len(keep_chains) > MAX_OLIGOMER_CHAINS:
            warnings.append(
                f"assembly has {len(keep_chains)} target copies (> "
                f"{MAX_OLIGOMER_CHAINS}); trimming to the target protomer only")
            keep_chains = [display_chains.target_chain]
    else:
        keep_chains = [target_chain]

    # 13. Trim the display structure (preserves auth numbering + icodes)
    trim = trim_structure(loaded.display, keep_chains, ecd_auth,
                          target_chain=display_chains.target_chain)
    warnings.extend(trim.warnings)

    # Re-label hotspots (and the patch) to the trimmed target chain id so the
    # chain-prefixed strings match the residue numbering in trimmed.pdb.
    for h in hotspots:
        h.chain = trim.target_chain
    for h in patch:
        h.chain = trim.target_chain

    return PipelineResult(
        record=rec, ecd_ranges=ecd_ranges, transmem_ranges=tm_ranges,
        membrane_proximal=mp, glyco_sites=glyco_sites, choice=choice,
        loaded=loaded, chains=chains, display_chains=display_chains,
        mapping=mapping, ecd_auth_ranges=ecd_auth,
        membrane_buffer_auth=membrane_buffer_auth, avoid=avoid,
        interface=interface, interface_chains=interface_chains,
        hotspots=hotspots, patch=patch, domains=domains, removals=removals, trim=trim,
        assembly=assembly, membrane_buffer=membrane_buffer, warnings=warnings,
    )
