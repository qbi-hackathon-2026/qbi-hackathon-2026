export function Header() {
  return (
    <header className="flex flex-col gap-2 sm:flex-row sm:items-end sm:justify-between sm:gap-6">
      <div className="flex items-baseline gap-4">
        <h1 className="text-3xl font-semibold tracking-tight sm:text-4xl">
          <span className="text-teal-600">Trim</span>
          <span className="text-slate-900">Prot</span>
        </h1>
        <p className="hidden text-sm text-slate-500 sm:block">
          Protein trimming for de novo design pipelines.
        </p>
      </div>
      <p className="text-sm text-slate-500 sm:hidden">
        Protein trimming for de novo design pipelines.
      </p>
    </header>
  );
}
