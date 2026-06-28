"""UniProt lookup and feature extraction for a target accession, plus search."""
import re
import requests

UNIPROT_REST = "https://rest.uniprot.org/uniprotkb"
ACCESSION_RE = re.compile(r"^[OPQ][0-9][A-Z0-9]{3}[0-9]$|^[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}$", re.IGNORECASE)


def search_proteins(query: str, limit: int = 12) -> list[dict]:
    """Resolve a free-text protein/gene name (or accession) to UniProt candidates.

    Human-tagged ("_HUMAN") entries are surfaced first per spec; other organisms
    (e.g. "_MOUSE") are still returned for visibility but ranked lower.
    """
    query = query.strip()
    if ACCESSION_RE.match(query):
        try:
            entry = fetch_entry(query.upper())
            return [_candidate_from_entry(entry)]
        except requests.HTTPError:
            pass  # fall through to text search

    r = requests.get(
        f"{UNIPROT_REST}/search",
        params={
            "query": f'(gene:{query} OR protein_name:"{query}") AND reviewed:true',
            "fields": "accession,id,protein_name,gene_names,organism_name,organism_id,length",
            "format": "json",
            "size": limit * 3,
        },
        timeout=20,
    )
    r.raise_for_status()
    results = r.json().get("results", [])

    candidates = [_candidate_from_entry(item) for item in results]

    q_upper = query.upper()
    candidates.sort(key=lambda c: (not c["is_human"], not c["id"].upper().startswith(q_upper + "_")))
    return candidates[:limit]


def _candidate_from_entry(entry: dict) -> dict:
    """Shape a single UniProt entry (search hit or full fetch) into a search candidate."""
    desc = entry.get("proteinDescription", {})
    protein_name = desc.get("recommendedName", {}).get("fullName", {}).get("value")
    if not protein_name:
        submitted = desc.get("submissionNames") or desc.get("submittedNames") or []
        if submitted:
            protein_name = submitted[0].get("fullName", {}).get("value")
    if not protein_name:
        protein_name = entry["uniProtkbId"]

    gene_names: list[str] = []
    for gene in entry.get("genes", []):
        name = gene.get("geneName", {}).get("value")
        if name and name not in gene_names:
            gene_names.append(name)
        for syn in gene.get("synonyms", []):
            value = syn.get("value")
            if value and value not in gene_names:
                gene_names.append(value)

    organism = entry.get("organism", {})
    return {
        "accession": entry["primaryAccession"],
        "id": entry["uniProtkbId"],
        "protein_name": protein_name,
        "gene_names": gene_names,
        "organism": organism.get("scientificName", ""),
        "is_human": organism.get("taxonId") == 9606,
        "length": entry.get("sequence", {}).get("length", 0),
    }

FEATURE_TYPES_OF_INTEREST = {
    "Signal peptide",
    "Transmembrane",
    "Topological domain",
    "Glycosylation",
    "Disulfide bond",
    "Modified residue",  # phosphorylation, acetylation, methylation, hydroxylation, etc.
    "Lipidation",
    "Cross-link",
}


def fetch_entry(accession: str) -> dict:
    r = requests.get(f"{UNIPROT_REST}/{accession}.json", timeout=20)
    r.raise_for_status()
    return r.json()


def _loc(feature: dict) -> tuple[int, int]:
    start = feature["location"]["start"]["value"]
    end = feature["location"]["end"]["value"]
    return start, end


def extract_features(entry: dict) -> dict:
    """Pull out the feature categories needed to decide on / perform trimming."""
    sequence = entry["sequence"]["value"]
    out = {
        "accession": entry["primaryAccession"],
        "id": entry["uniProtkbId"],
        "length": entry["sequence"]["length"],
        "sequence": sequence,
        "signal_peptide": None,
        "transmembrane_regions": [],
        "topological_domains": [],
        "glycosylation_sites": [],
        "disulfide_bonds": [],
        "other_ptms": [],
    }

    for f in entry.get("features", []):
        ftype = f["type"]
        if ftype not in FEATURE_TYPES_OF_INTEREST:
            continue
        start, end = _loc(f)
        if ftype == "Signal peptide":
            out["signal_peptide"] = (start, end)
        elif ftype == "Transmembrane":
            out["transmembrane_regions"].append((start, end))
        elif ftype == "Topological domain":
            out["topological_domains"].append(
                {"start": start, "end": end, "description": f.get("description", "")}
            )
        elif ftype == "Glycosylation":
            out["glycosylation_sites"].append(
                {"position": start, "description": f.get("description", "")}
            )
        elif ftype == "Disulfide bond":
            out["disulfide_bonds"].append((start, end))
        elif ftype in ("Modified residue", "Lipidation", "Cross-link"):
            out["other_ptms"].append({
                "type": ftype,
                "start": start,
                "end": end,
                "description": f.get("description", ""),
            })

    return out


def needs_trimming(features: dict) -> tuple[bool, str]:
    """Decide whether extracellular-domain trimming is needed for this target."""
    if features["transmembrane_regions"]:
        return True, (
            f"Target has {len(features['transmembrane_regions'])} transmembrane "
            f"region(s) ({features['transmembrane_regions']}); an extracellular-only "
            f"domain must be isolated for binder design."
        )
    return False, "No transmembrane region annotated; target may already be extracellular."


def get_extracellular_domain(features: dict) -> dict | None:
    """Return the extracellular topological domain boundaries (UniProt numbering)."""
    extracellular = [
        d for d in features["topological_domains"]
        if d["description"].strip().lower() == "extracellular"
    ]
    if not extracellular:
        return None
    # Largest extracellular span (handles multi-pass TM proteins with several entries).
    domain = max(extracellular, key=lambda d: d["end"] - d["start"])
    return domain


if __name__ == "__main__":
    entry = fetch_entry("P00533")
    feats = extract_features(entry)
    trim, reason = needs_trimming(feats)
    print("needs_trimming:", trim, "-", reason)
    print("extracellular domain:", get_extracellular_domain(feats))
    print("glyco sites:", len(feats["glycosylation_sites"]))
    print("disulfide bonds:", len(feats["disulfide_bonds"]))
