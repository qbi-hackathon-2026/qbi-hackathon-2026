import { Dna } from "lucide-react";
import { ColumnCard } from "./ColumnCard";

export function OriginColumn() {
  return (
    <ColumnCard
      title="Origin"
      dotClassName="bg-slate-400"
      icon={<Dna className="h-10 w-10" strokeWidth={1.5} />}
      emptyText="Select a protein to view origin data."
    />
  );
}
