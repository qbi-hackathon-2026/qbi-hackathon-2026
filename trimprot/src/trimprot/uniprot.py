"""UniProt resolution + feature parsing.

resolve_uniprot(name) -> best reviewed human accession with sequence + the
features the rest of the pipeline needs (TRANSMEM, TOPO_DOM, CARBOHYD, DISULFID,
MOD_RES, SIGNAL). A direct accession override is also supported.

All positions are 1-based UniProt numbering.
"""
from __future__ import annotations

import re
from dataclasses import dataclass, field
from typing import Optional

from .cache import get_json

UNIPROT_SEARCH = "https://rest.uniprot.org/uniprotkb/search"
UNIPROT_ENTRY = "https://rest.uniprot.org/uniprotkb/{acc}.json"

# UniProt JSON feature "type" strings we care about, mapped to short canonical keys.
_FEATURE_TYPES = {
    "Transmembrane": "TRANSMEM",
    "Topological domain": "TOPO_DOM",
    "Glycosylation": "CARBOHYD",
    "Disulfide bond": "DISULFID",
    "Modified residue": "MOD_RES",
    "Signal": "SIGNAL",
}

# Common aliases -> gene symbol, to make name resolution robust.
SYNONYMS = {
    "DC-SIGN": "CD209",
    "DCSIGN": "CD209",
    "CLEC4L": "CD209",
    "TNFSF7": "CD70",
    "CD27L": "CD70",
    "CD27LG": "CD70",
    "ERBB1": "EGFR",
    "HER1": "EGFR",
    "ERBB": "EGFR",
}


@dataclass(frozen=True)
class Feature:
    kind: str          # canonical key, e.g. "TRANSMEM"
    start: int         # 1-based UniProt
    end: int
    description: str = ""

    @property
    def length(self) -> int:
        return self.end - self.start + 1


@dataclass
class UniProtRecord:
    accession: str
    gene: str
    protein_name: str
    organism_id: int
    sequence: str
    function: str = ""
    features: list[Feature] = field(default_factory=list)

    def of(self, kind: str) -> list[Feature]:
        return [f for f in self.features if f.kind == kind]


def _parse_features(raw_features: list[dict]) -> list[Feature]:
    out: list[Feature] = []
    for f in raw_features:
        kind = _FEATURE_TYPES.get(f.get("type", ""))
        if kind is None:
            continue
        loc = f.get("location", {})
        try:
            start = int(loc["start"]["value"])
            end = int(loc["end"]["value"])
        except (KeyError, TypeError, ValueError):
            continue
        desc = f.get("description", "") or ""
        if kind == "DISULFID":
            # Disulfide features encode the two bonded cysteines as start & end.
            out.append(Feature("DISULFID", start, start, desc))
            out.append(Feature("DISULFID", end, end, desc))
        else:
            out.append(Feature(kind, start, end, desc))
    return out


def _protein_name(desc: dict) -> str:
    rec = desc.get("recommendedName") or {}
    full = (rec.get("fullName") or {}).get("value")
    if full:
        return full
    subs = desc.get("submissionNames") or []
    if subs:
        return (subs[0].get("fullName") or {}).get("value", "")
    return ""


def _function_text(entry: dict) -> str:
    for c in entry.get("comments", []) or []:
        if c.get("commentType") == "FUNCTION":
            texts = c.get("texts") or []
            if texts:
                return texts[0].get("value", "") or ""
    return ""


def _record_from_entry(entry: dict) -> UniProtRecord:
    acc = entry["primaryAccession"]
    genes = entry.get("genes") or []
    gene = ""
    if genes:
        gene = (genes[0].get("geneName") or {}).get("value", "")
    return UniProtRecord(
        accession=acc,
        gene=gene,
        protein_name=_protein_name(entry.get("proteinDescription", {})),
        organism_id=int((entry.get("organism") or {}).get("taxonId", 0) or 0),
        sequence=(entry.get("sequence") or {}).get("value", ""),
        function=_function_text(entry),
        features=_parse_features(entry.get("features", [])),
    )


def fetch_entry(accession: str) -> UniProtRecord:
    """Fetch a full UniProt entry by accession."""
    data = get_json(UNIPROT_ENTRY.format(acc=accession))
    return _record_from_entry(data)


def _search_accession(name: str, organism_id: int) -> str:
    """Search reviewed human entries over gene AND protein-name fields."""
    term = SYNONYMS.get(name.strip().upper(), name.strip())
    # Query both gene and protein name; restrict to reviewed + organism.
    query = (
        f"(gene:{term} OR gene_exact:{term} OR protein_name:{term}) "
        f"AND organism_id:{organism_id} AND reviewed:true"
    )
    params = {
        "query": query,
        "fields": "accession,gene_names,protein_name,annotation_score,organism_id",
        "format": "json",
        "size": "25",
    }
    data = get_json(UNIPROT_SEARCH, params=params)
    results = data.get("results", [])
    if not results:
        raise ValueError(f"no reviewed UniProt entry found for {name!r} (organism {organism_id})")

    upper = term.upper()

    def rank(r: dict) -> tuple:
        genes = r.get("genes") or []
        gene_syms = {(_g.get("geneName") or {}).get("value", "").upper() for _g in genes}
        # synonyms too
        for g in genes:
            for syn in g.get("synonyms", []) or []:
                gene_syms.add(syn.get("value", "").upper())
        exact_gene = upper in gene_syms
        score = float(r.get("annotationScore", 0) or 0)
        return (exact_gene, score)

    results.sort(key=rank, reverse=True)
    return results[0]["primaryAccession"]


def resolve_uniprot(name: Optional[str] = None, *, accession: Optional[str] = None,
                    organism_id: int = 9606) -> UniProtRecord:
    """Resolve a gene/protein name (or explicit accession) to a UniProtRecord.

    - accession override short-circuits the search.
    - otherwise searches reviewed entries for `organism_id` over gene + protein
      name fields, applying a small synonym map, and returns the best match.
    """
    if accession:
        return fetch_entry(accession)
    if not name:
        raise ValueError("resolve_uniprot requires either name or accession")
    acc = _search_accession(name, organism_id)
    return fetch_entry(acc)


# ── Typeahead search (UI dropdown) ────────────────────────────────────────────
# Official UniProtKB accession pattern.
ACCESSION_RE = re.compile(
    r"^([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})$",
    re.IGNORECASE,
)


def _candidate_from_search_hit(item: dict) -> dict:
    """Shape a UniProt search/entry JSON record into a typeahead candidate.

    Returns the field names the frontend dropdown consumes: accession, id,
    protein_name, gene_names, organism, is_human, length.
    """
    desc = item.get("proteinDescription", {})
    protein_name = (desc.get("recommendedName", {}).get("fullName", {}).get("value"))
    if not protein_name:
        submitted = desc.get("submissionNames") or desc.get("submittedNames") or []
        if submitted:
            protein_name = submitted[0].get("fullName", {}).get("value")
    if not protein_name:
        protein_name = item.get("uniProtkbId", item.get("primaryAccession", ""))

    gene_names: list[str] = []
    for gene in item.get("genes", []):
        gname = (gene.get("geneName") or {}).get("value")
        if gname and gname not in gene_names:
            gene_names.append(gname)
        for syn in gene.get("synonyms", []) or []:
            sval = syn.get("value")
            if sval and sval not in gene_names:
                gene_names.append(sval)

    organism = item.get("organism", {})
    return {
        "accession": item.get("primaryAccession", ""),
        "id": item.get("uniProtkbId", ""),
        "protein_name": protein_name,
        "gene_names": gene_names,
        "organism": organism.get("scientificName", ""),
        "is_human": organism.get("taxonId") == 9606,
        "length": item.get("sequence", {}).get("length", 0),
    }


def search_proteins(query: str, limit: int = 12) -> list[dict]:
    """Resolve free-text (gene/protein name or accession) to UniProt candidates.

    Human entries are surfaced first; an exact accession short-circuits to that
    single entry. Used by the web UI's typeahead dropdown.
    """
    query = (query or "").strip()
    if not query:
        return []

    if ACCESSION_RE.match(query):
        try:
            entry = get_json(UNIPROT_ENTRY.format(acc=query.upper()), allow_404=True)
            if entry:
                return [_candidate_from_search_hit(entry)]
        except Exception:
            pass  # fall through to text search

    data = get_json(
        UNIPROT_SEARCH,
        params={
            "query": f'(gene:{query} OR protein_name:"{query}") AND reviewed:true',
            "fields": "accession,id,protein_name,gene_names,organism_name,organism_id,length",
            "format": "json",
            "size": str(limit * 3),
        },
    )
    results = data.get("results", []) if data else []
    candidates = [_candidate_from_search_hit(item) for item in results]

    q_upper = query.upper()
    candidates.sort(key=lambda c: (not c["is_human"],
                                   not c["id"].upper().startswith(q_upper + "_")))
    return candidates[:limit]
