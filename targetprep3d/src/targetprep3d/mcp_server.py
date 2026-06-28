"""FastMCP server exposing every pipeline stage as a tool.

Same deterministic code paths as the CLI — an agent can orchestrate the
individual stages or run the whole pipeline. Run with:

    python -m targetprep3d.mcp_server
"""
from __future__ import annotations

from dataclasses import asdict
from pathlib import Path
from typing import Optional

from fastmcp import FastMCP

from .avoid import build_avoid_set
from .glyco import predict_glycosylation
from .hotspots import filter_hotspots
from .interface import detect_interface
from .outputs import build_summary, emit_outputs, format_hotspot_string
from .pipeline import run_pipeline
from .sifts import get_sifts_mapping
from .structio import classify_chains, load_structure
from .structures import search_structures
from .topology import (
    Range,
    get_extracellular_ranges,
    get_membrane_proximal_terminus,
    get_transmem_ranges,
)
from .uniprot import resolve_uniprot

mcp = FastMCP("targetprep3d")

DEFAULT_OUTDIR = Path(__file__).resolve().parents[2] / "outputs"


@mcp.tool
def resolve_target(name: Optional[str] = None, accession: Optional[str] = None,
                   organism_id: int = 9606) -> dict:
    """Resolve a gene/protein name (or accession) to a UniProt record + features."""
    rec = resolve_uniprot(name, accession=accession, organism_id=organism_id)
    return {
        "accession": rec.accession, "gene": rec.gene,
        "protein_name": rec.protein_name, "organism_id": rec.organism_id,
        "sequence_length": len(rec.sequence),
        "features": {k: [(f.start, f.end, f.description) for f in rec.of(k)]
                     for k in ("TRANSMEM", "TOPO_DOM", "CARBOHYD", "DISULFID",
                               "MOD_RES", "SIGNAL")},
    }


@mcp.tool
def extracellular_ranges(accession: str) -> dict:
    """UniProt extracellular (TOPO_DOM) ranges + transmembrane ranges."""
    rec = resolve_uniprot(accession=accession)
    ecd = get_extracellular_ranges(rec.features)
    tm = get_transmem_ranges(rec.features)
    return {"ecd_ranges": [(r.start, r.end) for r in ecd],
            "transmem_ranges": [(r.start, r.end) for r in tm]}


@mcp.tool
def membrane_proximal(accession: str, membrane_buffer: int = 12) -> dict:
    """Membrane-proximal ECD terminus + excluded buffer (UniProt numbering)."""
    rec = resolve_uniprot(accession=accession)
    ecd = get_extracellular_ranges(rec.features)
    tm = get_transmem_ranges(rec.features)
    mp = get_membrane_proximal_terminus(ecd, tm, buffer=membrane_buffer)
    if mp is None:
        return {"membrane_proximal": None}
    return {"terminus": mp.terminus, "buffer": [mp.buffer.start, mp.buffer.end],
            "topo_type": mp.topo_type, "agrees": mp.agrees, "note": mp.note}


@mcp.tool
def glycosylation(accession: str) -> dict:
    """Predicted glycosylation sites: CARBOHYD union N-X-[S/T] sequon scan."""
    rec = resolve_uniprot(accession=accession)
    sites = predict_glycosylation(rec.sequence, rec.of("CARBOHYD"))
    return {"sites_uniprot": sorted(sites)}


@mcp.tool
def structures(accession: str, prefer_antibody: bool = False) -> dict:
    """Ranked candidate structures with the chosen target/partner chains."""
    rec = resolve_uniprot(accession=accession)
    ecd = get_extracellular_ranges(rec.features)
    choice = search_structures(accession, ecd, prefer_antibody=prefer_antibody)
    return {
        "chosen_pdb": choice.chosen.pdb_id, "target_chain": choice.target_chain,
        "partner_chains": choice.partner_chains,
        "antibody_chains": choice.antibody_chains,
        "apo_fallback": choice.apo_fallback, "reasons": choice.reasons,
        "candidates": [{"pdb_id": c.pdb_id, "chain": c.chain_id,
                        "partner_tier": c.partner_tier,
                        "ecd_coverage": round(c.ecd_coverage, 3),
                        "completeness": round(c.completeness, 3),
                        "method": c.method_label, "resolution": c.resolution,
                        "antibody": c.antibody_partner} for c in choice.candidates[:15]],
    }


@mcp.tool
def sifts_mapping(pdb_id: str, accession: str, assembly: str = "bioassembly") -> dict:
    """Per-residue UniProt<->PDB(auth) mapping with observed flags + offset."""
    loaded = load_structure(pdb_id, assembly=assembly)
    rec = resolve_uniprot(accession=accession)
    ecd = get_extracellular_ranges(rec.features)
    choice = search_structures(accession, ecd)
    chains = classify_chains(loaded.structure, choice.target_chain)
    m = get_sifts_mapping(pdb_id, accession, structure=loaded.structure,
                          target_chain=chains.target_chain)
    chain = chains.target_chain
    return {
        "target_chain": chain,
        "consistent_offset": m.consistent_offset(chain),
        "n_residues": len(m.for_chain(chain)),
        "n_observed": sum(1 for r in m.for_chain(chain) if r.observed),
        "xref_mismatches": m.xref_mismatches,
        "sample": [{"unp": r.unp, "pdb_auth": r.pdb_num, "icode": r.icode,
                    "observed": r.observed} for r in m.for_chain(chain)[:10]],
    }


@mcp.tool
def prepare_target(name: Optional[str] = None, accession: Optional[str] = None,
                   prefer_antibody: bool = False, assembly: str = "bioassembly",
                   membrane_buffer: int = 12, patch_radius: float = 11.0,
                   patch_size: int = 8, no_patch: bool = False,
                   outdir: Optional[str] = None) -> dict:
    """Run the full pipeline and emit the design-ready artifact set."""
    result = run_pipeline(name, accession=accession, prefer_antibody=prefer_antibody,
                          assembly=assembly, membrane_buffer=membrane_buffer,
                          patch_radius=patch_radius, patch_size=patch_size,
                          no_patch=no_patch)
    out = emit_outputs(result, Path(outdir) if outdir else DEFAULT_OUTDIR)
    summary = build_summary(result)
    return {"output_dir": str(out),
            "hotspot_string_full": format_hotspot_string(result.hotspots),
            "patch_string": format_hotspot_string(result.patch),
            "summary": summary}


def main() -> None:
    mcp.run()


if __name__ == "__main__":
    main()
