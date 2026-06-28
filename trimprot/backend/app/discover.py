"""Disease -> candidate target-protein discovery.

Given a free-text disease (and optional cell type), resolve it against the Open
Targets Platform, pull disease-associated targets ranked by evidence, and return
the ones that have a UniProt accession (so they feed straight into the TrimProt
pipeline). An optional LLM re-ranking step can reorder/annotate that fixed
candidate list when an ANTHROPIC_API_KEY is available.

Design notes:
- Open Targets is the source of truth. The deterministic ranking is fully usable
  on its own; the LLM only reorders and annotates the list it's given and can
  never introduce a gene/accession we didn't already verify.
- "Upstream/driver" is approximated by up-weighting the genetic_association and
  somatic_mutation datatype scores, which are causal evidence, over downstream
  markers (expression, known drug, etc.).
"""
import json
import os

import requests

OPENTARGETS_GQL = "https://api.platform.opentargets.org/api/v4/graphql"

# Datatype-score weights for the deterministic "driver-ness" re-rank. Causal
# evidence (germline genetics, somatic driver mutations) is weighted above
# associative/downstream evidence.
DATATYPE_WEIGHTS = {
    "genetic_association": 1.0,
    "somatic_mutation": 1.0,
    "known_drug": 0.5,
    "affected_pathway": 0.5,
    "literature": 0.3,
    "rna_expression": 0.2,
    "animal_model": 0.2,
}

ANTHROPIC_MESSAGES_URL = (
    os.environ.get("ANTHROPIC_BASE_URL", "https://api.anthropic.com").rstrip("/")
    + "/v1/messages"
)
LLM_MODEL = "claude-opus-4-8"


def resolve_disease(text: str) -> dict | None:
    """Resolve free-text disease name to an Open Targets disease (EFO/MONDO id)."""
    query = """
    query($q: String!) {
      search(queryString: $q, entityNames: ["disease"], page: {index: 0, size: 1}) {
        hits { id name entity }
      }
    }"""
    r = requests.post(
        OPENTARGETS_GQL,
        json={"query": query, "variables": {"q": text}},
        timeout=20,
    )
    r.raise_for_status()
    hits = r.json()["data"]["search"]["hits"]
    if not hits:
        return None
    return {"efo_id": hits[0]["id"], "name": hits[0]["name"]}


def disease_targets(efo_id: str, size: int = 25) -> list[dict]:
    """Fetch disease-associated targets with per-datatype evidence scores.

    Returns one row per target that has a reviewed (Swiss-Prot) UniProt accession,
    carrying the Open Targets overall score, the datatype breakdown, and a
    derived `driver_score` (datatype scores weighted by DATATYPE_WEIGHTS).
    """
    query = """
    query($efoId: String!, $size: Int!) {
      disease(efoId: $efoId) {
        name
        associatedTargets(page: {index: 0, size: $size}) {
          rows {
            score
            datatypeScores { id score }
            target {
              approvedSymbol
              approvedName
              proteinIds { id source }
            }
          }
        }
      }
    }"""
    r = requests.post(
        OPENTARGETS_GQL,
        json={"query": query, "variables": {"efoId": efo_id, "size": size}},
        timeout=25,
    )
    r.raise_for_status()
    disease = r.json()["data"]["disease"]
    if not disease:
        return []

    out = []
    for row in disease["associatedTargets"]["rows"]:
        accession = _pick_uniprot(row["target"]["proteinIds"])
        if not accession:
            continue
        datatypes = {d["id"]: d["score"] for d in row["datatypeScores"]}
        driver = sum(DATATYPE_WEIGHTS.get(k, 0.1) * v for k, v in datatypes.items())
        out.append({
            "symbol": row["target"]["approvedSymbol"],
            "name": row["target"]["approvedName"],
            "accession": accession,
            "association_score": round(row["score"], 3),
            "driver_score": round(driver, 3),
            "evidence": {k: round(v, 3) for k, v in sorted(datatypes.items(), key=lambda x: -x[1])},
        })

    out.sort(key=lambda t: t["driver_score"], reverse=True)
    return out


def _pick_uniprot(protein_ids: list[dict]) -> str | None:
    """Prefer a reviewed Swiss-Prot accession; fall back to TrEMBL if that's all."""
    swissprot = [p["id"] for p in protein_ids if p["source"] == "uniprot_swissprot"]
    if swissprot:
        return swissprot[0]
    trembl = [p["id"] for p in protein_ids if p["source"] == "uniprot_trembl"]
    return trembl[0] if trembl else None


def llm_rerank(disease_name: str, cell_type: str | None, targets: list[dict]) -> dict:
    """Optionally reorder/annotate the candidate list with an LLM.

    The model is constrained to the accessions we pass in: it returns an ordering
    of those accessions plus a one-line rationale each (factoring in cell-type
    relevance, which Open Targets can't). Any accession the model returns that
    isn't in our list is ignored, and any we passed that it omits is appended in
    the original order, so the set of candidates can never change - only the order
    and the rationales. Returns {"used": bool, "targets": [...], "note": str}.
    """
    api_key = os.environ.get("ANTHROPIC_API_KEY")
    if not api_key:
        return {"used": False, "targets": targets,
                "note": "Deterministic ranking (set ANTHROPIC_API_KEY to enable LLM re-ranking)."}

    catalog = [
        {"accession": t["accession"], "symbol": t["symbol"],
         "association_score": t["association_score"], "driver_score": t["driver_score"],
         "evidence": t["evidence"]}
        for t in targets
    ]
    cell = cell_type or "(unspecified)"
    prompt = (
        f"You are helping a protein engineer pick de novo binder design targets for the "
        f"disease '{disease_name}' in cell type '{cell}'.\n\n"
        f"Here is a fixed list of candidate targets from Open Targets (with evidence "
        f"scores). Reorder them best-first for THIS disease and cell type, preferring "
        f"upstream drivers and targets relevant to the cell type. You may ONLY use the "
        f"accessions given; do not invent any.\n\n"
        f"{json.dumps(catalog, indent=2)}\n\n"
        f"Respond with ONLY a JSON array of objects "
        f'[{{"accession": "...", "rationale": "<=15 words"}}], best first.'
    )

    try:
        r = requests.post(
            ANTHROPIC_MESSAGES_URL,
            headers={
                "x-api-key": api_key,
                "anthropic-version": "2023-06-01",
                "content-type": "application/json",
            },
            json={
                "model": LLM_MODEL,
                "max_tokens": 1500,
                "messages": [{"role": "user", "content": prompt}],
            },
            timeout=60,
        )
        r.raise_for_status()
        text = r.json()["content"][0]["text"]
        ranking = json.loads(_extract_json_array(text))
    except Exception as e:
        return {"used": False, "targets": targets,
                "note": f"LLM re-ranking unavailable ({type(e).__name__}); showing deterministic ranking."}

    by_acc = {t["accession"]: t for t in targets}
    ordered, seen = [], set()
    for item in ranking:
        acc = item.get("accession")
        if acc in by_acc and acc not in seen:
            t = dict(by_acc[acc])
            t["llm_rationale"] = item.get("rationale", "")
            ordered.append(t)
            seen.add(acc)
    # Append any candidates the model dropped, so the set is preserved.
    for t in targets:
        if t["accession"] not in seen:
            ordered.append(t)
    return {"used": True, "targets": ordered,
            "note": f"Re-ranked by {LLM_MODEL} for cell type '{cell}' (candidate set unchanged)."}


def _extract_json_array(text: str) -> str:
    """Pull the first [...] block out of an LLM response (handles code fences/prose)."""
    start = text.find("[")
    end = text.rfind("]")
    if start == -1 or end == -1 or end < start:
        raise ValueError("no JSON array in LLM response")
    return text[start:end + 1]


def discover(disease_text: str, cell_type: str | None = None, size: int = 25,
             use_llm: bool = True) -> dict:
    """Top-level: disease text -> ranked candidate targets ready for the pipeline."""
    disease = resolve_disease(disease_text)
    if disease is None:
        return {"disease": None, "query": disease_text, "targets": [],
                "note": f"No disease matched '{disease_text}'."}

    targets = disease_targets(disease["efo_id"], size=size)
    rerank = llm_rerank(disease["name"], cell_type, targets) if use_llm else \
        {"used": False, "targets": targets, "note": "Deterministic ranking (LLM disabled)."}

    return {
        "disease": disease,
        "query": disease_text,
        "cell_type": cell_type,
        "llm_used": rerank["used"],
        "note": rerank["note"],
        "targets": rerank["targets"],
    }


if __name__ == "__main__":
    import sys
    q = sys.argv[1] if len(sys.argv) > 1 else "melanoma"
    result = discover(q, cell_type=sys.argv[2] if len(sys.argv) > 2 else None, size=10)
    print("disease:", result["disease"], "| llm_used:", result["llm_used"])
    print("note:", result["note"])
    for t in result["targets"]:
        print(f"  {t['symbol']:>8} {t['accession']:>8} assoc={t['association_score']} "
              f"driver={t['driver_score']} {t.get('llm_rationale', '')}")
