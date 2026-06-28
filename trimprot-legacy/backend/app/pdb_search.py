"""Search and rank PDB structures for a UniProt target, batched via RCSB APIs."""
import requests

RCSB_SEARCH = "https://search.rcsb.org/rcsbsearch/v2/query"
RCSB_GRAPHQL = "https://data.rcsb.org/graphql"

PARTNER_KEYWORDS = ("antibody", "fab", "nanobody", "scfv", "vhh", "fragment, antigen-binding")

ENTRY_QUERY = """
query ($ids: [String!]!) {
  entries(entry_ids: $ids) {
    rcsb_id
    rcsb_entry_info {
      resolution_combined
      polymer_entity_count_protein
    }
    polymer_entities {
      rcsb_id
      rcsb_polymer_entity_container_identifiers {
        auth_asym_ids
        reference_sequence_identifiers {
          database_accession
          database_name
        }
      }
      rcsb_polymer_entity { pdbx_description }
      entity_poly { pdbx_seq_one_letter_code_can }
    }
  }
}
"""


def find_entity_ids(accession: str) -> list[str]:
    """All polymer-entity ids (PDB ENTRY_ENTITY) referencing this UniProt accession."""
    query = {
        "query": {
            "type": "terminal",
            "service": "text",
            "parameters": {
                "attribute": "rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_accession",
                "operator": "exact_match",
                "value": accession,
            },
        },
        "return_type": "polymer_entity",
        "request_options": {"results_content_type": ["experimental"], "return_all_hits": True},
    }
    r = requests.post(RCSB_SEARCH, json=query, timeout=30)
    r.raise_for_status()
    if r.status_code == 204 or not r.content:
        return []
    hits = r.json().get("result_set", [])
    return [h["identifier"] for h in hits]


def fetch_entry_details(entry_ids: list[str], batch_size: int = 50) -> list[dict]:
    """Batched GraphQL fetch of resolution, entity descriptions and sequences."""
    out = []
    for i in range(0, len(entry_ids), batch_size):
        batch = entry_ids[i:i + batch_size]
        r = requests.post(RCSB_GRAPHQL, json={"query": ENTRY_QUERY, "variables": {"ids": batch}}, timeout=30)
        r.raise_for_status()
        entries = r.json().get("data", {}).get("entries", []) or []
        out.extend(e for e in entries if e is not None)
    return out


def _has_partner(entry: dict, target_accession: str) -> bool:
    for pe in entry.get("polymer_entities", []):
        desc = (pe.get("rcsb_polymer_entity") or {}).get("pdbx_description", "") or ""
        if any(k in desc.lower() for k in PARTNER_KEYWORDS):
            return True
        refs = (pe.get("rcsb_polymer_entity_container_identifiers") or {}).get(
            "reference_sequence_identifiers"
        ) or []
        accs = {r["database_accession"] for r in refs if r.get("database_name") == "UniProt"}
        if accs and target_accession not in accs and desc:
            # A non-target polymer entity present alongside the target = a bound partner.
            return True
    return False


def _target_coverage(entry: dict, target_accession: str, domain_length: int) -> float:
    best = 0
    for pe in entry.get("polymer_entities", []):
        refs = (pe.get("rcsb_polymer_entity_container_identifiers") or {}).get(
            "reference_sequence_identifiers"
        ) or []
        if any(r.get("database_accession") == target_accession for r in refs):
            seq = (pe.get("entity_poly") or {}).get("pdbx_seq_one_letter_code_can", "") or ""
            best = max(best, len(seq))
    return min(best / domain_length, 1.0) if domain_length else 0.0


def rank_structures(entries: list[dict], target_accession: str, domain_length: int) -> list[dict]:
    """Score: ECD coverage (most important) + has-partner bonus - poor resolution penalty."""
    scored = []
    for entry in entries:
        info = entry.get("rcsb_entry_info") or {}
        res_list = info.get("resolution_combined") or []
        resolution = res_list[0] if res_list else None
        coverage = _target_coverage(entry, target_accession, domain_length)
        has_partner = _has_partner(entry, target_accession)

        score = coverage * 10
        score += 2.0 if has_partner else 0.0
        if resolution:
            score += max(0.0, (4.0 - resolution))  # better (lower) resolution -> higher score

        scored.append({
            "pdb_id": entry["rcsb_id"],
            "resolution": resolution,
            "ecd_coverage": round(coverage, 3),
            "has_partner": has_partner,
            "score": round(score, 3),
        })

    scored.sort(key=lambda x: x["score"], reverse=True)
    return scored


if __name__ == "__main__":
    accession = "P00533"
    domain_length = 645 - 25 + 1  # EGFR ECD per UniProt topological domain
    entity_ids = find_entity_ids(accession)
    entry_ids = sorted({eid.split("_")[0] for eid in entity_ids})
    print(f"{len(entry_ids)} candidate PDB entries for {accession}")

    details = fetch_entry_details(entry_ids)
    ranked = rank_structures(details, accession, domain_length)
    for r in ranked[:10]:
        print(r)
