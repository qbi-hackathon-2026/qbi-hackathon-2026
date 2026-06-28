"""Command-line interface.

    python -m trimprot.cli --protein EGFR [--accession P00533]
        [--prefer-antibody] [--assembly bioassembly] [--membrane-buffer 12]
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

from .outputs import emit_outputs, format_hotspot_string
from .pipeline import PipelineResult, run_pipeline

DEFAULT_OUTDIR = Path(__file__).resolve().parents[2] / "outputs"


def _summary_row(result: PipelineResult, outdir: Path) -> dict:
    mp = result.membrane_proximal
    return {
        "target": result.record.gene or result.record.accession,
        "accession": result.record.accession,
        "pdb": result.choice.chosen.pdb_id,
        "chain": result.chains.target_chain,
        "method": result.choice.chosen.method_display,
        "hotspots": len(result.hotspots),
        "avoid": len(result.avoid),
        "apo_fallback": result.choice.apo_fallback or not result.chains.partners,
        "mp_terminus": mp.terminus if mp else "-",
        "outdir": str(outdir),
    }


def _print_table(rows: list[dict]) -> None:
    cols = ["target", "accession", "pdb", "chain", "method", "hotspots", "avoid",
            "apo_fallback", "mp_terminus", "outdir"]
    widths = {c: max(len(c), *(len(str(r[c])) for r in rows)) for c in cols}
    line = "  ".join(c.ljust(widths[c]) for c in cols)
    print("\n" + line)
    print("  ".join("-" * widths[c] for c in cols))
    for r in rows:
        print("  ".join(str(r[c]).ljust(widths[c]) for c in cols))


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(prog="trimprot",
                                description="Deterministic protein target prep for binder design.")
    p.add_argument("--protein", help="gene or protein name (e.g. EGFR)")
    p.add_argument("--accession", help="explicit UniProt accession override (e.g. P00533)")
    p.add_argument("--organism-id", type=int, default=9606)
    p.add_argument("--prefer-antibody", action="store_true",
                   help="boost antibody-bound complexes in structure ranking")
    p.add_argument("--assembly", choices=["bioassembly", "protomer"], default="bioassembly")
    p.add_argument("--membrane-buffer", type=int, default=12,
                   help="ECD residues nearest the TM excluded from hotspots")
    p.add_argument("--interface-cutoff", type=float, default=5.0)
    p.add_argument("--min-ecd-coverage", type=float, default=0.40,
                   help="minimum ECD coverage a structure must have to be eligible")
    p.add_argument("--patch-radius", type=float, default=11.0,
                   help="Cβ contiguity radius (Å) for the BindCraft epitope patch")
    p.add_argument("--patch-size", type=int, default=8,
                   help="max residues in the BindCraft epitope patch")
    p.add_argument("--no-patch", action="store_true",
                   help="emit the full ranked hotspot list to BindCraft (no patch reduction)")
    p.add_argument("--outdir", default=str(DEFAULT_OUTDIR))
    return p


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    if not args.protein and not args.accession:
        print("error: provide --protein or --accession", file=sys.stderr)
        return 2

    result = run_pipeline(
        args.protein, accession=args.accession, organism_id=args.organism_id,
        prefer_antibody=args.prefer_antibody, assembly=args.assembly,
        membrane_buffer=args.membrane_buffer, interface_cutoff=args.interface_cutoff,
        patch_radius=args.patch_radius, patch_size=args.patch_size,
        no_patch=args.no_patch, min_ecd_coverage=args.min_ecd_coverage,
    )
    outdir = emit_outputs(result, Path(args.outdir))

    for w in result.warnings:
        print(f"[warn] {w}", file=sys.stderr)
    _print_table([_summary_row(result, outdir)])
    print(f"\nfull hotspots ({len(result.hotspots)}): "
          f"{format_hotspot_string(result.hotspots) or '(empty)'}")
    print(f"epitope patch -> bindcraft ({len(result.patch)}): "
          f"{format_hotspot_string(result.patch) or '(empty)'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
