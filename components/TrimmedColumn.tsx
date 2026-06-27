import { Scissors } from "lucide-react";
import { ColumnCard } from "./ColumnCard";

export function TrimmedColumn() {
  return (
    <ColumnCard
      title="Trimmed"
      dotClassName="bg-teal-500"
      icon={<Scissors className="h-10 w-10" strokeWidth={1.5} />}
      emptyText="Select a protein to view trimmed sequence."
    />
  );
}
