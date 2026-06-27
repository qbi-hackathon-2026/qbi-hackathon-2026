import type { ProteinSearchResult } from "./types";

interface UniProtValue<T = string> {
  value: T;
}

interface UniProtName {
  fullName?: UniProtValue;
  shortNames?: UniProtValue[];
}

interface UniProtGene {
  geneName?: UniProtValue;
  synonyms?: UniProtValue[];
  orderedLocusNames?: UniProtValue[];
  orfNames?: UniProtValue[];
}

interface UniProtEntry {
  primaryAccession?: string;
  uniProtkbId?: string;
  entryType?: string;
  proteinDescription?: {
    recommendedName?: UniProtName;
    submissionNames?: UniProtName[];
    submittedNames?: UniProtName[];
    alternativeNames?: UniProtName[];
  };
  genes?: UniProtGene[];
  organism?: {
    scientificName?: string;
    commonName?: string;
  };
  sequence?: {
    length?: number;
  };
}

interface UniProtSearchResponse {
  results?: UniProtEntry[];
}

function pickProteinName(entry: UniProtEntry): string {
  const desc = entry.proteinDescription;
  const rec = desc?.recommendedName?.fullName?.value;
  if (rec) return rec;
  const submitted =
    desc?.submittedNames?.[0]?.fullName?.value ??
    desc?.submissionNames?.[0]?.fullName?.value;
  if (submitted) return submitted;
  const alt = desc?.alternativeNames?.[0]?.fullName?.value;
  if (alt) return alt;
  return entry.uniProtkbId ?? entry.primaryAccession ?? "Unknown protein";
}

function pickGeneNames(entry: UniProtEntry): string[] {
  const out = new Set<string>();
  for (const g of entry.genes ?? []) {
    if (g.geneName?.value) out.add(g.geneName.value);
    for (const s of g.synonyms ?? []) {
      if (s.value) out.add(s.value);
    }
  }
  return Array.from(out);
}

export function parseUniProtResponse(
  raw: UniProtSearchResponse
): ProteinSearchResult[] {
  const entries = raw.results ?? [];
  return entries.map((entry) => ({
    accession: entry.primaryAccession ?? "",
    id: entry.uniProtkbId ?? "",
    proteinName: pickProteinName(entry),
    geneNames: pickGeneNames(entry),
    organism: entry.organism?.scientificName ?? "",
    length: entry.sequence?.length ?? 0,
    reviewed: entry.entryType === "UniProtKB reviewed (Swiss-Prot)",
  }));
}

// UniProtKB accession format, per https://www.uniprot.org/help/accession_numbers
const ACCESSION_RE = /^[OPQ][0-9][A-Z0-9]{3}[0-9]$|^[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2}$/i;
// Entry name (UniProtKB ID), e.g. "INS_HUMAN"
const ENTRY_ID_RE = /^[A-Z0-9]{1,10}_[A-Z0-9]{1,5}$/i;

export function buildUniProtSearchUrl(query: string): string {
  const q = query.trim();
  const clauses = [`gene:${q}`, `protein_name:${q}`];
  if (ACCESSION_RE.test(q)) clauses.push(`accession:${q}`);
  if (ENTRY_ID_RE.test(q)) clauses.push(`id:${q}`);
  const lucene = `(${clauses.join(" OR ")}) AND reviewed:true`;
  const params = new URLSearchParams({
    query: lucene,
    fields: "accession,id,protein_name,gene_names,organism_name,length,reviewed",
    size: "10",
    format: "json",
  });
  return `https://rest.uniprot.org/uniprotkb/search?${params.toString()}`;
}
