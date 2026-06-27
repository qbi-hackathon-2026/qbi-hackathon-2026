import { NextResponse } from "next/server";
import { buildUniProtSearchUrl, parseUniProtResponse } from "@/lib/uniprot";
import type { SearchResponse } from "@/lib/types";

export const runtime = "nodejs";
export const dynamic = "force-dynamic";

export async function GET(request: Request): Promise<NextResponse<SearchResponse>> {
  const { searchParams } = new URL(request.url);
  const q = (searchParams.get("q") ?? "").trim();

  if (q.length < 2) {
    return NextResponse.json({ results: [] });
  }

  try {
    const url = buildUniProtSearchUrl(q);
    const upstream = await fetch(url, {
      headers: { Accept: "application/json" },
      signal: AbortSignal.timeout(8000),
    });

    if (!upstream.ok) {
      return NextResponse.json({
        results: [],
        error: `UniProt responded with ${upstream.status}`,
      });
    }

    const json = await upstream.json();
    const results = parseUniProtResponse(json);
    return NextResponse.json({ results });
  } catch (err) {
    const message = err instanceof Error ? err.message : "Unknown error";
    return NextResponse.json({ results: [], error: message });
  }
}
