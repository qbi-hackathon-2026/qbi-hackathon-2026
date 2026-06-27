import Link from "next/link";

export default function NotFound() {
  return (
    <main className="min-h-screen">
      <div className="mx-auto max-w-2xl px-6 py-24 text-center">
        <p className="text-sm font-medium uppercase tracking-wider text-teal-600">
          404
        </p>
        <h1 className="mt-3 text-3xl font-semibold text-slate-900">
          Page not found
        </h1>
        <p className="mt-3 text-slate-600">
          The page you were looking for doesn&apos;t exist.
        </p>
        <Link
          href="/"
          className="mt-8 inline-block rounded-md bg-slate-900 px-4 py-2 text-sm font-medium text-white hover:bg-slate-800"
        >
          Back to TrimProt
        </Link>
      </div>
    </main>
  );
}
