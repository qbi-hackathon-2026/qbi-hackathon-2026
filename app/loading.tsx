export default function Loading() {
  return (
    <main className="min-h-screen">
      <div className="mx-auto max-w-7xl px-6 py-10 lg:px-10">
        <div className="h-10 w-48 animate-pulse rounded bg-slate-200" />
        <div className="mt-4 h-4 w-80 animate-pulse rounded bg-slate-200" />
        <div className="mt-10 h-11 w-full max-w-2xl animate-pulse rounded-md bg-slate-200" />
        <div className="mt-12 grid grid-cols-1 gap-6 lg:grid-cols-2">
          <div className="h-72 animate-pulse rounded-xl bg-slate-200" />
          <div className="h-72 animate-pulse rounded-xl bg-slate-200" />
        </div>
      </div>
    </main>
  );
}
