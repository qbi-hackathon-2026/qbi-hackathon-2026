"""Emit the design-ready artifact set for a PipelineResult.

Writes: original.cif, trimmed.pdb, summary.json, hotspots.csv,
avoid_residues.csv, bindcraft_config.json, viewer.html.
"""
from __future__ import annotations

import csv
import json
from pathlib import Path

from .hotspots import Hotspot
from .pipeline import PipelineResult
from .structures import MIN_ECD_COVERAGE
from .viewer import build_viewer_html


def _auth_label(num: int, icode: str) -> str:
    return f"{num}{icode}" if icode else f"{num}"


def format_hotspot_string(hotspots: list[Hotspot]) -> str:
    """Chain-prefixed, range-compressed auth numbering: 'A180,A182,A190-A193'."""
    by_chain: dict[str, list[Hotspot]] = {}
    for h in hotspots:
        by_chain.setdefault(h.chain, []).append(h)

    tokens: list[str] = []
    for chain in sorted(by_chain):
        plain = sorted({h.num for h in by_chain[chain] if not h.icode})
        iced = sorted((h.num, h.icode) for h in by_chain[chain] if h.icode)
        # compress consecutive plain numbers
        i = 0
        while i < len(plain):
            j = i
            while j + 1 < len(plain) and plain[j + 1] == plain[j] + 1:
                j += 1
            if j > i:
                tokens.append(f"{chain}{plain[i]}-{chain}{plain[j]}")
            else:
                tokens.append(f"{chain}{plain[i]}")
            i = j + 1
        for num, icode in iced:
            tokens.append(f"{chain}{num}{icode}")
    return ",".join(tokens)


def build_summary(result: PipelineResult) -> dict:
    rec = result.record
    mp = result.membrane_proximal
    choice = result.choice
    return {
        "target": rec.gene or rec.accession,
        "accession": rec.accession,
        "protein_name": rec.protein_name,
        "function": rec.function,
        "organism_id": rec.organism_id,
        "viewer": {
            "target_chain": result.trim.target_chain,
            "kept_chains": result.trim.kept_chains,
            "hotspot_auth_residues": sorted({h.num for h in result.hotspots}),
            "patch_auth_residues": sorted({h.num for h in result.patch}),
            "domains": result.domains,
        },
        "topology": {
            "type": mp.topo_type if mp else None,
            "ecd_ranges_uniprot": [(r.start, r.end) for r in result.ecd_ranges],
            "transmem_ranges_uniprot": [(r.start, r.end) for r in result.transmem_ranges],
            "membrane_proximal_terminus": mp.terminus if mp else None,
            "membrane_proximal_agrees_with_topology": mp.agrees if mp else None,
            "membrane_buffer_residues": result.membrane_buffer,
            "membrane_buffer_uniprot": ([mp.buffer.start, mp.buffer.end] if mp else None),
            "membrane_buffer_auth": ([result.membrane_buffer_auth.start,
                                      result.membrane_buffer_auth.end]
                                     if result.membrane_buffer_auth else None),
            "assumption": (mp.note if mp else
                           "no TRANSMEM/TOPO_DOM annotation; membrane buffer skipped"),
        },
        "structure": {
            "chosen_pdb": choice.chosen.pdb_id,
            "target_chain": result.chains.target_chain,
            "partner_chains": result.chains.partners,
            "antibody_chains": result.choice.antibody_chains,
            "interface_chains": result.interface_chains,
            "homo_oligomer_chains": result.chains.same_protein,
            "assembly": result.assembly,
            "assembly_applied": result.loaded.assembly_applied,
            "apo_fallback": result.choice.apo_fallback or not result.chains.partners,
            "resolution": choice.chosen.resolution,
            "method": choice.chosen.method,
            "method_label": choice.chosen.method_display,
            "predicted": choice.chosen.predicted,
            "ecd_coverage": round(choice.chosen.ecd_coverage, 3),
            "completeness": round(choice.chosen.completeness, 3),
            "partner_tier": choice.chosen.partner_tier,
            "min_ecd_coverage": MIN_ECD_COVERAGE,
            "selection_reasons": choice.reasons,
            "candidates_top10": [
                {"pdb_id": c.pdb_id, "chain": c.chain_id,
                 "partner_tier": c.partner_tier,
                 "ecd_coverage": round(c.ecd_coverage, 3),
                 "completeness": round(c.completeness, 3),
                 "method": c.method_label, "resolution": c.resolution,
                 "predicted": c.predicted, "antibody": c.antibody_partner}
                for c in choice.candidates[:10]
            ],
        },
        "sifts": {
            "target_offset": result.mapping.consistent_offset(result.chains.target_chain),
            "n_mapped_residues": len(result.mapping.for_chain(result.chains.target_chain)),
            "xref_mismatches": result.mapping.xref_mismatches,
        },
        "ecd_auth_ranges": [(r.start, r.end) for r in result.ecd_auth_ranges],
        "glycosylation_sites_uniprot": sorted(result.glyco_sites),
        "counts": {
            "interface_candidates": len(result.interface),
            "hotspots": len(result.hotspots),
            "patch": len(result.patch),
            "avoid": len(result.avoid),
            "removed_by_avoid": len(result.removals),
        },
        "hotspot_string_full": format_hotspot_string(result.hotspots),
        "patch_string": format_hotspot_string(result.patch),
        "patch_residues": [
            {"residue": _auth_label(h.num, h.icode), "uniprot": h.uniprot,
             "contact_count": h.contact_count} for h in result.patch
        ],
        "removed_hotspots": [
            {"residue": _auth_label(r.num, r.icode), "reasons": r.reasons}
            for r in result.removals
        ],
        "trim": {"kept_chains": result.trim.kept_chains,
                 "kept_counts": result.trim.kept_counts},
        "warnings": result.warnings,
        "alphafold": result.alphafold.summary_fields() if result.alphafold else None,
    }


def build_bindcraft_config(result: PipelineResult, outdir: Path) -> dict:
    target = result.record.gene or result.record.accession
    # Only the focused epitope patch goes to BindCraft (full list lives in
    # hotspots.csv); designing against one contiguous patch is what BindCraft wants.
    return {
        "design_path": str(outdir),
        "binder_name": f"{target}_binder",
        "starting_pdb": "trimmed.pdb",
        "chains": ",".join(result.trim.kept_chains),
        "target_hotspot_residues": format_hotspot_string(result.patch),
        "lengths": [65, 150],
        "number_of_final_designs": 100,
    }


def emit_outputs(result: PipelineResult, base_dir: Path) -> Path:
    target = result.record.gene or result.record.accession
    outdir = Path(base_dir) / target
    outdir.mkdir(parents=True, exist_ok=True)

    # original assembly (mmCIF)
    original_path = outdir / "original.cif"
    result.loaded.display.make_mmcif_document().write_file(str(original_path))

    # trimmed target (PDB)
    trimmed_path = outdir / "trimmed.pdb"
    result.trim.structure.write_pdb(str(trimmed_path))

    # summary.json
    summary = build_summary(result)
    (outdir / "summary.json").write_text(json.dumps(summary, indent=2))

    # hotspots.csv — FULL ranked list, with an in_patch flag for the BindCraft patch
    patch_keys = {(h.num, h.icode.strip()) for h in result.patch}
    with (outdir / "hotspots.csv").open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["chain", "pdb_auth_residue", "uniprot_residue",
                    "contact_count", "source", "in_patch"])
        for h in result.hotspots:
            in_patch = (h.num, h.icode.strip()) in patch_keys
            w.writerow([h.chain, _auth_label(h.num, h.icode),
                        h.uniprot if h.uniprot is not None else "",
                        h.contact_count, h.source, in_patch])

    # avoid_residues.csv
    with (outdir / "avoid_residues.csv").open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["chain", "pdb_auth_residue", "uniprot_residue", "reasons"])
        out_chain = result.trim.target_chain
        for e in result.avoid.entries():
            w.writerow([out_chain, _auth_label(e.num, e.icode),
                        e.uniprot if e.uniprot is not None else "",
                        ";".join(e.reasons)])

    # bindcraft_config.json
    config = build_bindcraft_config(result, outdir)
    (outdir / "bindcraft_config.json").write_text(json.dumps(config, indent=2))

    # viewer.html (self-contained): trimmed target with hotspot residues highlighted
    html = build_viewer_html(
        target=target,
        pdb_id=result.choice.chosen.pdb_id,
        pdb_text=trimmed_path.read_text(),
        hotspot_residues=[h.num for h in result.hotspots],
    )
    (outdir / "viewer.html").write_text(html)

    return outdir
