import { Header } from "@/components/Header";
import { SearchBox } from "@/components/SearchBox";
import { ResultColumns } from "@/components/ResultColumns";

export default function HomePage() {
  return (
    <main className="min-h-screen">
      <div className="mx-auto max-w-7xl px-6 py-10 sm:py-14 lg:px-10">
        <Header />
        <section className="mt-10 max-w-2xl">
          <SearchBox />
        </section>
        <section className="mt-12">
          <ResultColumns />
        </section>
      </div>
    </main>
  );
}
