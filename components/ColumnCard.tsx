import type { ReactNode } from "react";
import { cn } from "@/lib/cn";

interface ColumnCardProps {
  title: string;
  dotClassName: string;
  icon: ReactNode;
  emptyText: string;
}

export function ColumnCard({
  title,
  dotClassName,
  icon,
  emptyText,
}: ColumnCardProps) {
  return (
    <div className="flex h-80 flex-col rounded-xl border border-slate-200 bg-white p-6 shadow-sm sm:h-96">
      <div className="flex items-center gap-2">
        <span
          aria-hidden
          className={cn("h-2.5 w-2.5 rounded-full", dotClassName)}
        />
        <h2 className="text-sm font-semibold uppercase tracking-wider text-slate-700">
          {title}
        </h2>
      </div>
      <div className="mt-4 flex flex-1 flex-col items-center justify-center gap-3 text-center">
        <div className="text-slate-300">{icon}</div>
        <p className="max-w-xs text-sm text-slate-500">{emptyText}</p>
      </div>
    </div>
  );
}
