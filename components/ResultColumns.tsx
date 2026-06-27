import { OriginColumn } from "./OriginColumn";
import { TrimmedColumn } from "./TrimmedColumn";

export function ResultColumns() {
  return (
    <div className="grid grid-cols-1 gap-6 lg:grid-cols-2">
      <OriginColumn />
      <TrimmedColumn />
    </div>
  );
}
