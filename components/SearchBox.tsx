"use client";

import {
  useCallback,
  useEffect,
  useId,
  useRef,
  useState,
  type KeyboardEvent,
} from "react";
import { Loader2, Search, X } from "lucide-react";
import { cn } from "@/lib/cn";
import type { ProteinSearchResult, SearchResponse } from "@/lib/types";

type Status = "idle" | "loading" | "ready" | "error";

export function SearchBox() {
  const [query, setQuery] = useState("");
  const [results, setResults] = useState<ProteinSearchResult[]>([]);
  const [status, setStatus] = useState<Status>("idle");
  const [open, setOpen] = useState(false);
  const [activeIndex, setActiveIndex] = useState(0);
  const [selected, setSelected] = useState<ProteinSearchResult | null>(null);

  const inputRef = useRef<HTMLInputElement>(null);
  const containerRef = useRef<HTMLDivElement>(null);
  const debounceRef = useRef<ReturnType<typeof setTimeout> | null>(null);
  const abortRef = useRef<AbortController | null>(null);
  const listboxId = useId();

  const runSearch = useCallback(async (q: string) => {
    abortRef.current?.abort();
    const controller = new AbortController();
    abortRef.current = controller;

    setStatus("loading");
    try {
      const res = await fetch(`/api/search?q=${encodeURIComponent(q)}`, {
        signal: controller.signal,
      });
      const data: SearchResponse = await res.json();
      if (controller.signal.aborted) return;
      setResults(data.results);
      setActiveIndex(0);
      setStatus(data.error ? "error" : "ready");
    } catch (err) {
      if ((err as Error).name === "AbortError") return;
      setResults([]);
      setStatus("error");
    }
  }, []);

  useEffect(() => {
    if (debounceRef.current) clearTimeout(debounceRef.current);

    const trimmed = query.trim();
    // If the input shows a selection, don't re-search until the user edits.
    if (selected && query === formatSelection(selected)) {
      return;
    }

    if (trimmed.length < 2) {
      abortRef.current?.abort();
      setResults([]);
      setStatus("idle");
      return;
    }

    debounceRef.current = setTimeout(() => {
      void runSearch(trimmed);
    }, 250);

    return () => {
      if (debounceRef.current) clearTimeout(debounceRef.current);
    };
  }, [query, runSearch, selected]);

  useEffect(() => {
    function onPointerDown(e: MouseEvent) {
      if (!containerRef.current?.contains(e.target as Node)) {
        setOpen(false);
      }
    }
    window.addEventListener("mousedown", onPointerDown);
    return () => window.removeEventListener("mousedown", onPointerDown);
  }, []);

  function handleSelect(result: ProteinSearchResult) {
    setSelected(result);
    setQuery(formatSelection(result));
    setOpen(false);
    setResults([]);
    setStatus("idle");
    inputRef.current?.blur();
    // eslint-disable-next-line no-console
    console.log("Selected protein:", result);
  }

  function handleKeyDown(e: KeyboardEvent<HTMLInputElement>) {
    if (e.key === "Escape") {
      setOpen(false);
      return;
    }
    if (!open && (e.key === "ArrowDown" || e.key === "ArrowUp")) {
      setOpen(true);
      return;
    }
    if (e.key === "ArrowDown") {
      e.preventDefault();
      if (results.length === 0) return;
      setActiveIndex((i) => (i + 1) % results.length);
    } else if (e.key === "ArrowUp") {
      e.preventDefault();
      if (results.length === 0) return;
      setActiveIndex((i) => (i - 1 + results.length) % results.length);
    } else if (e.key === "Enter") {
      if (results[activeIndex]) {
        e.preventDefault();
        handleSelect(results[activeIndex]);
      }
    }
  }

  function clearInput() {
    setQuery("");
    setResults([]);
    setSelected(null);
    setStatus("idle");
    inputRef.current?.focus();
  }

  const showDropdown = open && query.trim().length >= 2;

  return (
    <div ref={containerRef} className="relative">
      <label htmlFor="protein-search" className="sr-only">
        Search proteins
      </label>
      <div
        className={cn(
          "flex items-center gap-2 rounded-lg border border-slate-200 bg-white px-3 py-2.5 shadow-sm transition",
          "focus-within:border-teal-500 focus-within:ring-2 focus-within:ring-teal-500/20"
        )}
      >
        <Search className="h-4 w-4 shrink-0 text-slate-400" aria-hidden />
        <input
          id="protein-search"
          ref={inputRef}
          type="text"
          autoComplete="off"
          spellCheck={false}
          role="combobox"
          aria-expanded={showDropdown}
          aria-controls={listboxId}
          aria-autocomplete="list"
          aria-activedescendant={
            showDropdown && results[activeIndex]
              ? `${listboxId}-opt-${activeIndex}`
              : undefined
          }
          placeholder="Search by protein name, gene symbol, or UniProt accession..."
          value={query}
          onChange={(e) => {
            setQuery(e.target.value);
            setOpen(true);
            if (selected) setSelected(null);
          }}
          onFocus={() => setOpen(true)}
          onKeyDown={handleKeyDown}
          className="w-full bg-transparent text-sm text-slate-900 placeholder:text-slate-400 focus:outline-none"
        />
        {status === "loading" && (
          <Loader2 className="h-4 w-4 shrink-0 animate-spin text-slate-400" aria-hidden />
        )}
        {query.length > 0 && status !== "loading" && (
          <button
            type="button"
            onClick={clearInput}
            className="rounded p-0.5 text-slate-400 hover:bg-slate-100 hover:text-slate-600"
            aria-label="Clear search"
          >
            <X className="h-4 w-4" />
          </button>
        )}
      </div>

      {showDropdown && (
        <div
          className="absolute left-0 right-0 top-full z-20 mt-2 overflow-hidden rounded-lg border border-slate-200 bg-white shadow-lg"
          role="listbox"
          id={listboxId}
        >
          <DropdownContents
            status={status}
            results={results}
            query={query}
            activeIndex={activeIndex}
            listboxId={listboxId}
            onHover={setActiveIndex}
            onSelect={handleSelect}
          />
        </div>
      )}
    </div>
  );
}

function DropdownContents({
  status,
  results,
  query,
  activeIndex,
  listboxId,
  onHover,
  onSelect,
}: {
  status: Status;
  results: ProteinSearchResult[];
  query: string;
  activeIndex: number;
  listboxId: string;
  onHover: (i: number) => void;
  onSelect: (r: ProteinSearchResult) => void;
}) {
  if (status === "loading" && results.length === 0) {
    return (
      <div className="flex items-center gap-2 px-4 py-3 text-sm text-slate-500">
        <Loader2 className="h-4 w-4 animate-spin" aria-hidden />
        Searching UniProt…
      </div>
    );
  }

  if (status === "error" && results.length === 0) {
    return (
      <div className="px-4 py-3 text-sm text-rose-600">
        Couldn&apos;t reach UniProt. Try again in a moment.
      </div>
    );
  }

  if (results.length === 0) {
    return (
      <div className="px-4 py-3 text-sm text-slate-500">
        No proteins found matching &ldquo;{query.trim()}&rdquo;.
      </div>
    );
  }

  return (
    <ul className="max-h-96 overflow-y-auto py-1">
      {results.map((r, i) => {
        const active = i === activeIndex;
        return (
          <li
            key={`${r.accession}-${r.id}`}
            id={`${listboxId}-opt-${i}`}
            role="option"
            aria-selected={active}
            onMouseEnter={() => onHover(i)}
            onMouseDown={(e) => {
              e.preventDefault();
              onSelect(r);
            }}
            className={cn(
              "flex cursor-pointer items-start justify-between gap-4 px-4 py-2.5 text-sm",
              active ? "bg-slate-50" : "bg-white"
            )}
          >
            <div className="min-w-0 flex-1">
              <div className="flex flex-wrap items-baseline gap-x-2">
                <span className="truncate font-medium text-slate-900">
                  {r.proteinName}
                </span>
                {r.geneNames.length > 0 && (
                  <span className="truncate text-xs text-slate-500">
                    {r.geneNames.slice(0, 3).join(", ")}
                  </span>
                )}
              </div>
            </div>
            <div className="flex shrink-0 items-center gap-2 text-xs">
              <span className="italic text-slate-500">{r.organism}</span>
              <span className="rounded bg-slate-100 px-1.5 py-0.5 font-mono text-[11px] text-slate-700">
                {r.accession}
              </span>
              <span className="text-slate-400">{r.length} aa</span>
            </div>
          </li>
        );
      })}
    </ul>
  );
}

function formatSelection(r: ProteinSearchResult): string {
  return `${r.proteinName} (${r.accession})`;
}
