export interface ProteinSearchResult {
  accession: string;
  id: string;
  proteinName: string;
  geneNames: string[];
  organism: string;
  length: number;
  reviewed: boolean;
}

export interface SearchResponse {
  results: ProteinSearchResult[];
  error?: string;
}
