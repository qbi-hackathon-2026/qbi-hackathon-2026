"""UniProt resolution + feature parsing.

resolve_uniprot(name) -> best reviewed human accession with sequence + the
features the rest of the pipeline needs (TRANSMEM, TOPO_DOM, CARBOHYD, DISULFID,
MOD_RES, SIGNAL). A direct accession override is also supported.

All positions are 1-based UniProt numbering.
"""
from __future__ import annotations

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
