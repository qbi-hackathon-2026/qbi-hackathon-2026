<div align="center">

# TrimProt

**Deterministic protein target-preparation for de novo binder design.**

[Run it locally](#installation) ·
[How it works](#correctness-guarantees) ·
[Contributing](CONTRIBUTING.md)

</div>

Give TrimProt a protein (a gene name, protein name, or UniProt accession) and it
prepares a design-ready target: it resolves the UniProt entry, finds and ranks
structures, trims to the relevant domain, identifies candidate **hotspot**
residues and an epitope **patch**, builds an **avoid** set, and renders it all in
an interactive 3D viewer. When no suitable experimental structure exists, it falls
back to the **AlphaFold Database** automatically.

---

## Installation

With [uv](https://docs.astral.sh/uv/) and `git` installed:

```bash
git clone https://github.com/qbi-hackathon-2026/qbi-hackathon-2026.git
cd qbi-hackathon-2026/trimprot
uv sync
```

uv handles Python for you (no system interpreter needed; requires Python ≥ 3.11).

<details>
<summary><b>Never used a terminal? Full step-by-step for macOS</b></summary>

<br>

A "terminal" is a window where you type commands instead of clicking. Open it,
paste each line, and press **Enter** (paste with **Cmd+V**). Do these once.

**1. Open the Terminal.** Press **Cmd + Space**, type `Terminal`, press **Enter**.

**2. Install Homebrew** (a tool that installs other tools):
```bash
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
```
It may ask for your Mac password (typing shows nothing — that's normal). When it
finishes, **close the Terminal and open a new one** so the new tool is available.

**3. Install `git` and `uv`:**
```bash
brew install git uv
```

**4. Download the code and set it up:**
```bash
git clone https://github.com/qbi-hackathon-2026/qbi-hackathon-2026.git
cd qbi-hackathon-2026/trimprot
uv sync
```

You're ready — see [Running TrimProt](#running-trimprot) below.

</details>

<details>
<summary><b>Never used a terminal? Full step-by-step for Windows</b></summary>

<br>

A "terminal" is a window where you type commands instead of clicking. Open it,
paste each line, and press **Enter** (paste with **Ctrl+V**). Do these once.

**1. Open PowerShell.** Click **Start**, type `PowerShell`, press **Enter**.

**2. Install `git` and `uv`:**
```powershell
winget install --id Git.Git -e
winget install --id astral-sh.uv -e
```
If Windows asks for permission, allow it. When both finish, **close PowerShell and
open a new one** so the new tools are available.

**3. Download the code and set it up:**
```powershell
git clone https://github.com/qbi-hackathon-2026/qbi-hackathon-2026.git
cd qbi-hackathon-2026\trimprot
uv sync
```

You're ready — see [Running TrimProt](#running-trimprot) below.

</details>

---

## Running TrimProt

Start the web app from the `trimprot` folder:

```bash
uv run python -m trimprot.server
```

Leave that window open. When you see `Uvicorn running on http://127.0.0.1:8000`,
open **http://127.0.0.1:8000** in your browser.

Then:

1. **Search** a target — type a gene name, protein name, or UniProt accession
   (e.g. `EGFR` or `P00533`) and pick a result.
2. The pipeline runs and the page shows two 3D panels: the **chosen RCSB
   assembly** (colored by chain) beside the **trimmed target** — hotspots in
   green, the epitope patch in orange — with chain/domain labels, per-residue
   hover, and a "focus patch" button.
3. **Download** the prepared artifacts from the page (see [Output](#output)).

To stop the app, click its terminal window and press **Ctrl + C**.

> **Note:** `http://127.0.0.1:8000` runs on your own computer only — it's not
> reachable by others on the internet. Putting it online is a separate, more
> advanced setup.

<details>
<summary><b>If something goes wrong</b></summary>

<br>

- **`command not found` / `not recognized`** right after installing — open a
  **new** terminal window (the old one didn't pick up the newly installed tool).
- **The browser can't connect** — make sure the terminal still shows the
  `Uvicorn running…` line. If you closed it, run the start command again.
- **A target shows an error** — some proteins have no usable structure; try a
  well-known one like `EGFR` to confirm the app itself is working.

</details>

---

## Output

For each target, TrimProt emits files under `trimprot/outputs/<TARGET>/`,
downloadable from the page:

- `trimmed.pdb` — the design-ready trimmed target (author numbering preserved)
- `original.cif` — the chosen assembly (or AlphaFold model)
- `summary.json` — structured result: chosen structure + reasoning, topology, ECD
  ranges, hotspots, epitope patch, avoid set, warnings
- `hotspots.csv` — ranked hotspots with an `in_patch` flag
- a BindCraft config for the epitope patch

---

## How TrimProt picks a structure

<details>
<summary><b>How the selection ladder works</b></summary>

<br>

Candidates come from **PDBe `best_structures`** (UniProt-mapped), enriched via the
**RCSB Data API** (partner chains, method, resolution) and validated against
**PDBe SIFTS** numbering. Selection is a priority ladder, not a weighted
score: it keeps usable structures (valid numbering, sane resolution) that cover
the extracellular domain, prefers them by partner type (**antibody-bound ▸
ligand-bound ▸ apo**), then by coverage and completeness. Method and resolution
only break ties. The reasoning behind each pick is written to `summary.json`,
which can be downloaded.

</details>

---

## Correctness guarantees

<details>
<summary><b>The structural details TrimProt handles correctly</b></summary>

<br>

1. Everything operates on **author numbering** (`gemmi residue.seqid.num` + insertion code). Never renumbered.
2. UniProt↔PDB(auth) mapping comes from **PDBe SIFTS** (segment-based, offsets + gaps), cross-checked against `_pdbx_sifts_xref_db`.
3. Unobserved UniProt residues → `avoid`.
4. Topology from UniProt **TOPO_DOM** only (never inferred from position). Handles type I and type II correctly.
5. Homo-oligomers download the **biological assembly** by default.
6. **Membrane-proximal exclusion:** the ECD terminus adjacent to TRANSMEM is excluded by a buffer (default 12 residues).
7. **AlphaFold fallback:** when no experimental structure exists (or none covers
   ≥ 50% of the ECD), the AlphaFold model is used, with pLDDT masking: residues
   ≤ 50 treated as unobserved and 50–70 as low-confidence, both folded into the
   avoid set.

</details>

---

## Development

```bash
cd trimprot
uv sync --extra dev
uv run pytest -m "not network"     # offline unit tests (run in CI)
uv run pytest                       # full suite (hits UniProt/PDBe/RCSB)
```

See [`CONTRIBUTING.md`](CONTRIBUTING.md) for how to contribute.

---

## License

[MIT](LICENSE) © TrimProt contributors.
