#!/usr/bin/env python3
"""Read-only auditor for trimprot structure choices.

For every target under ./outputs/, read summary.json and independently re-query
the RCSB Data API (REST) to verify the chosen PDB: method/resolution, per-entity
identity, antibody detection (multi-signal), the designated partner chains, and
ligands. Nothing in the pipeline is touched or imported — this only reads
summary.json and talks to RCSB.

Exit status is nonzero if any target is "MISLABELED" (claims a partner-bound
complex whose partner is not actually an antibody) so it is CI-catchable.

    cd trimprot && ./.venv/bin/python validate_structures.py   # or: python validate_structures.py
"""
from __future__ import annotations

import hashlib
import json
import re
import sys
from pathlib import Path

import requests

ROOT = Path(__file__).resolve().parent
OUTPUTS = ROOT / "outputs"
CACHE = ROOT / "cache" / "rcsb"
REST = "https://data.rcsb.org/rest/v1/core"

_SESSION = requests.Session()
_SESSION.headers["User-Agent"] = "trimprot-validate/1.0"

# --- antibody signal vocabularies -------------------------------------------
# Immunoglobulin-fold domain annotations. NOTE: an Ig fold ALONE is ambiguous —
# many cell-surface receptors are Ig-superfamily — so it only yields AMBIGUOUS
# unless an explicit antibody name/UniProt signal also fires.
PFAM_IG = {"PF07686", "PF07654", "PF07679", "PF00047", "PF16190", "PF05790"}
INTERPRO_IG = {"IPR007110", "IPR013783", "IPR003599", "IPR013106",
               "IPR036179", "IPR003598", "IPR013098"}
# Explicit antibody name keywords (these identify an actual immunoglobulin chain).
AB_NAME_RE = re.compile(
    r"\b(fab|fv|scfv|nanobody|vhh|sybody|igg|igm|iga|ige)\b"
    r"|heavy chain|light chain|antibody|immunoglobulin", re.I)
# UniProt protein-name / id signal for immunoglobulin entries.
AB_UNP_RE = re.compile(r"immunoglobulin|\bIg[ GKLHMAE]", re.I)

GLYCANS = {"NAG", "NDG", "BMA", "MAN", "BGC", "GLC", "GAL", "GLA", "FUC", "FUL",
           "A2G", "NGA", "SIA", "NAN", "FUM", "XYP", "GCS", "MMA"}
IONS = {"NA", "CL", "K", "MG", "CA", "ZN", "MN", "FE", "CD", "NI", "CO", "CU",
        "SO4", "PO4", "ACT", "EDO", "GOL", "PEG", "PGE", "PG4", "SCN", "IOD",
        "BR", "FLC", "CIT", "MES", "TRS", "EPE", "DMS", "FMT", "NO3"}


# --- cached HTTP ------------------------------------------------------------
def get(url: str):
    """Cached GET returning parsed JSON, or None on 404 / any error."""
    CACHE.mkdir(parents=True, exist_ok=True)
    key = hashlib.sha256(url.encode()).hexdigest()[:20]
    path = CACHE / f"{key}.json"
    if path.exists():
        txt = path.read_text()
        return None if txt == "null" else json.loads(txt)
    try:
        resp = _SESSION.get(url, timeout=45)
    except requests.RequestException as exc:
        print(f"    [warn] request failed: {url} ({exc})")
        return None
    if resp.status_code == 404:
        path.write_text("null")
        return None
    if resp.status_code != 200:
        print(f"    [warn] {resp.status_code} for {url}")
        return None
    try:
        data = resp.json()
    except ValueError:
        return None
    path.write_text(json.dumps(data))
    return data


# --- RCSB accessors ---------------------------------------------------------
def entry_info(pdb: str) -> dict | None:
    return get(f"{REST}/entry/{pdb.lower()}")


def polymer_entity(pdb: str, eid: str) -> dict | None:
    return get(f"{REST}/polymer_entity/{pdb.lower()}/{eid}")


def nonpolymer_entity(pdb: str, eid: str) -> dict | None:
    return get(f"{REST}/nonpolymer_entity/{pdb.lower()}/{eid}")


def uniprot_core(pdb: str, eid: str) -> list | None:
    return get(f"{REST}/uniprot/{pdb.lower()}/{eid}")


# --- antibody classification ------------------------------------------------
def classify_entity(pdb: str, eid: str, pe: dict) -> dict:
    """Return identity + antibody verdict for one polymer entity."""
    rpe = pe.get("rcsb_polymer_entity", {}) or {}
    cid = pe.get("rcsb_polymer_entity_container_identifiers", {}) or {}
    name = rpe.get("pdbx_description", "") or ""
    auth = cid.get("auth_asym_ids", []) or []
    unp_ids = cid.get("uniprot_ids") or []
    refs = cid.get("reference_sequence_identifiers") or []
    if not unp_ids:
        unp_ids = [r.get("database_accession") for r in refs
                   if r.get("database_name") == "UniProt"]
    orgs = [o.get("scientific_name") or o.get("ncbi_scientific_name")
            for o in (pe.get("rcsb_entity_source_organism") or [])]

    signals: list[str] = []

    # (a) immunoglobulin-FOLD annotations. Only structural domain databases count
    # (GO terms like "B cell mediated immunity" carry "immunoglobulin" in their
    # lineage but say nothing about the fold, so they must be ignored).
    ig_anno = False
    structural = {"Pfam", "InterPro", "SCOP", "SCOP2", "SCOP2B", "CATH", "ECOD"}
    for a in pe.get("rcsb_polymer_entity_annotation") or []:
        atype = a.get("type", "") or ""
        if atype not in structural:
            continue
        aid = a.get("annotation_id", "") or ""
        aname = a.get("name", "") or ""
        lineage = " ".join((x.get("name") or "") for x in (a.get("annotation_lineage") or []))
        if (aid in PFAM_IG or aid in INTERPRO_IG
                or "immunoglobulin" in aname.lower()
                or "immunoglobulin" in lineage.lower()):
            ig_anno = True
            signals.append(f"annotation:{atype}:{aid or aname}")

    # (b) explicit antibody name keyword
    name_ab = bool(AB_NAME_RE.search(name))
    if name_ab:
        signals.append(f"name:'{name}'")

    # (c) UniProt maps to an immunoglobulin entry
    unp_ab = False
    unp_names = []
    for acc in unp_ids:
        if not acc:
            continue
        core = uniprot_core(pdb, eid)
        if isinstance(core, list):
            for o in core:
                nm = ((o.get("rcsb_uniprot_protein") or {}).get("name") or {}).get("value", "")
                uid = ((o.get("rcsb_uniprot_container_identifiers") or {})
                       .get("uniprot_id", ""))
                unp_names.append(nm or uid)
                if AB_UNP_RE.search(nm or "") or AB_UNP_RE.search(uid or ""):
                    unp_ab = True
                    signals.append(f"uniprot:{uid}:'{nm}'")
        break  # one lookup per entity is enough

    # verdict: an explicit antibody identity (name or UniProt) is decisive; an Ig
    # fold by itself is only AMBIGUOUS (could be an Ig-superfamily receptor).
    if name_ab or unp_ab:
        verdict = "ANTIBODY"
    elif ig_anno:
        verdict = "AMBIGUOUS"
    else:
        verdict = "NON-ANTIBODY"

    return {
        "entity_id": eid, "name": name, "auth_chains": auth,
        "uniprot": [a for a in unp_ids if a], "organisms": [o for o in orgs if o],
        "verdict": verdict, "signals": signals,
    }


def collect_ligands(pdb: str, entry: dict) -> list[dict]:
    out = []
    for eid in (entry.get("rcsb_entry_container_identifiers", {})
                .get("nonpolymer_entity_ids") or []):
        ne = nonpolymer_entity(pdb, eid)
        if not ne:
            continue
        comp = (ne.get("pdbx_entity_nonpoly") or {})
        cid = comp.get("comp_id", "?")
        nm = comp.get("name", "") or (ne.get("rcsb_nonpolymer_entity") or {}).get("pdbx_description", "")
        kind = ("glycosylation" if cid in GLYCANS
                else "buffer/ion" if cid in IONS else "functional?")
        out.append({"comp_id": cid, "name": nm, "kind": kind})
    return out


# --- per-target audit -------------------------------------------------------
def audit_target(target_dir: Path) -> dict:
    summary = json.loads((target_dir / "summary.json").read_text())
    st = summary.get("structure", {})
    target = summary.get("target", target_dir.name)
    target_unp = (summary.get("accession") or "").upper()
    pdb = (st.get("chosen_pdb") or "").lower()
    partner_chains = st.get("partner_chains", []) or []
    apo_fallback = bool(st.get("apo_fallback"))
    sum_res = st.get("resolution")
    sum_method = (st.get("method") or "")

    print("=" * 88)
    print(f"TARGET {target}  ({summary.get('accession','?')})   chosen PDB: {pdb.upper()}"
          f"   target_chain: {st.get('target_chain')}")
    print("-" * 88)

    result = {"target": target, "pdb": pdb.upper(), "partner_verdict": "?",
              "resolution_ok": None, "n_partner_chains": len(partner_chains),
              "n_antibody_partners": 0, "mislabeled": False}

    entry = entry_info(pdb)
    if entry is None:
        print("  [error] entry not available from RCSB (404 or fetch failed) — "
              "reporting only what summary.json claims.")
        print(f"  summary says: method={sum_method!r} resolution={sum_res}")
        result["partner_verdict"] = "UNVERIFIABLE (entry 404)"
        return result

    # 1. entry-level
    title = (entry.get("struct") or {}).get("title", "")
    methods = [m.get("method") for m in (entry.get("exptl") or [])]
    res_combined = (entry.get("rcsb_entry_info") or {}).get("resolution_combined")
    rcsb_res = res_combined[0] if res_combined else None
    print(f"  title: {title}")
    print(f"  method (RCSB): {methods}   resolution (RCSB): {rcsb_res}")
    print(f"  summary claims: method={sum_method!r}  resolution={sum_res}")

    res_ok = True
    if sum_res is not None and rcsb_res is not None and abs(sum_res - rcsb_res) > 0.1:
        print(f"  [FLAG] resolution mismatch: summary {sum_res} vs RCSB {rcsb_res}")
        res_ok = False
    if methods and sum_method and not any(
            sum_method.lower().replace("-", " ").split()[0] in (m or "").lower()
            for m in methods):
        print(f"  [FLAG] method mismatch: summary {sum_method!r} vs RCSB {methods}")
        res_ok = False
    result["resolution_ok"] = res_ok

    # 2. per-entity identity + antibody verdict
    eids = (entry.get("rcsb_entry_container_identifiers") or {}).get("polymer_entity_ids") or []
    chain_to_entity: dict[str, dict] = {}
    print("  entities:")
    for eid in eids:
        pe = polymer_entity(pdb, eid)
        if pe is None:
            print(f"    entity {eid}: [unavailable]")
            continue
        info = classify_entity(pdb, eid, pe)
        for ch in info["auth_chains"]:
            chain_to_entity[ch] = info
        org = ", ".join(info["organisms"]) or "?"
        unp = ",".join(info["uniprot"]) or "-"
        print(f"    entity {eid} [{info['verdict']}]  chains={info['auth_chains']}  "
              f"unp={unp}  org={org}")
        print(f"        name: {info['name']}")
        if info["signals"]:
            print(f"        signals: {'; '.join(info['signals'])}")

    # 3. partner verdict
    print("  designated partner_chains:", partner_chains or "(none)")
    partner_infos = []
    for ch in partner_chains:
        info = chain_to_entity.get(ch)
        if info is None:
            print(f"    chain {ch}: [no entity found in RCSB record]")
            partner_infos.append(None)
            continue
        partner_infos.append(info)
        print(f"    chain {ch}: entity {info['entity_id']} -> {info['verdict']}  "
              f"({info['name']})")

    present = [i for i in partner_infos if i]
    n_ab = sum(1 for i in present if i["verdict"] == "ANTIBODY")
    result["n_antibody_partners"] = n_ab

    # A "partner" that is actually the SAME protein as the target is a genuine
    # mislabel (a homo-oligomer copy counted as a partner). A different protein
    # is a legitimate partner: antibody -> antibody-bound, otherwise a
    # receptor/ligand-bound complex. Both are valid epitope sources.
    selfish = [i for i in present if target_unp
               and target_unp in {u.upper() for u in i["uniprot"]}]
    if apo_fallback:
        sneaky = [i for i in present if i["verdict"] == "ANTIBODY"]
        verdict = ("apo BUT an antibody chain is present (inconsistent)" if sneaky
                   else "apo (no partner expected)")
    elif not present:
        verdict = "MISLABELED: apo_fallback=False but no partner entity resolved"
        result["mislabeled"] = True
    elif selfish:
        names = sorted({i["name"] for i in selfish})
        verdict = f"MISLABELED: 'partner' is another copy of the target -> {', '.join(names)}"
        result["mislabeled"] = True
    elif all(i["verdict"] == "ANTIBODY" for i in present):
        verdict = "antibody-bound (consistent)"
    else:
        names = sorted({i["name"] for i in present})
        verdict = f"receptor_ligand-bound (consistent) -> {', '.join(names)}"
    print(f"  >>> PARTNER VERDICT: {verdict}")
    result["partner_verdict"] = verdict

    # 4. ligands
    ligs = collect_ligands(pdb, entry)
    if ligs:
        print("  ligands/heteroatoms:")
        for lg in ligs:
            print(f"    {lg['comp_id']:4} [{lg['kind']}]  {lg['name']}")
    else:
        print("  ligands/heteroatoms: none reported")
    return result


def main() -> int:
    if not OUTPUTS.is_dir():
        print(f"no outputs directory at {OUTPUTS}")
        return 2
    target_dirs = sorted(d for d in OUTPUTS.iterdir()
                         if d.is_dir() and (d / "summary.json").exists())
    if not target_dirs:
        print("no target summaries found under outputs/")
        return 2

    rows = []
    for d in target_dirs:
        try:
            rows.append(audit_target(d))
        except Exception as exc:  # never let one target crash the whole run
            print(f"  [error] auditing {d.name} failed: {exc}")
            rows.append({"target": d.name, "pdb": "?", "partner_verdict": f"ERROR: {exc}",
                         "resolution_ok": None, "n_partner_chains": 0,
                         "n_antibody_partners": 0, "mislabeled": False})

    # final summary table
    print("\n" + "=" * 88)
    print("SUMMARY")
    print("-" * 88)
    hdr = f"{'target':8} {'pdb':6} {'res_ok':6} {'#part':5} {'#ab':4}  partner_verdict"
    print(hdr)
    print("-" * 88)
    for r in rows:
        ok = "-" if r["resolution_ok"] is None else ("yes" if r["resolution_ok"] else "NO")
        print(f"{r['target']:8} {r['pdb']:6} {ok:6} {r['n_partner_chains']:<5} "
              f"{r['n_antibody_partners']:<4}  {r['partner_verdict']}")

    mislabeled = [r["target"] for r in rows if r["mislabeled"]]
    print("-" * 88)
    if mislabeled:
        print(f"MISLABELED targets (exit 1): {', '.join(mislabeled)}")
        return 1
    print("all partner verdicts consistent")
    return 0


if __name__ == "__main__":
    sys.exit(main())
