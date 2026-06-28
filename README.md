<div align="center">

# Bonsai

**Deterministic protein target-preparation for de novo binder design.**

[Run it locally](#installation) ·
[How it works](#correctness-guarantees) ·
[Contributing](CONTRIBUTING.md)

</div>

To design a new binder (an antibody or mini-protein that sticks to a disease
target), you first have to prepare that target structurally. Today a structural
biologist does this by hand, and it can take hours per target.

Bonsai automates it. Give Bonsai a protein (a gene name, protein name, or
UniProt accession) and it prepares a design-ready target: it resolves the UniProt
entry, finds and ranks structures, trims to the relevant domain, identifies
candidate **hotspot** residues and an epitope **patch**, builds an **avoid** set,
and renders it all in an interactive 3D viewer. When no suitable experimental
structure exists, it falls back to the **AlphaFold Database** automatically. A
built-in **chat assistant** lets you ask why a structure was chosen or re-run with
new settings in plain English. The output drops straight into de novo
binder-design pipelines like **BindCraft**, **RFdiffusion**, or **BoltzGen**.

---

## Installation

With [uv](https://docs.astral.sh/uv/) and `git` installed:

```bash
git clone https://github.com/qbi-hackathon-2026/qbi-hackathon-2026.git
cd qbi-hackathon-2026/bindscout
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
cd qbi-hackathon-2026/bindscout
uv sync
```

You're ready — see [Running Bonsai](#running-bonsai) below.

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
cd qbi-hackathon-2026\bindscout
uv sync
```

You're ready — see [Running Bonsai](#running-bonsai) below.

</details>

---

## Running Bonsai

Start the web app from the `bindscout` folder:

```bash
uv run python -m bindscout.server
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
4. **Ask Bonsai** (optional) — click the chat panel to ask about the loaded
   target or re-run the pipeline with different settings in plain English. The
   assistant needs an Anthropic API key: copy `bindscout/.env.example` to
   `bindscout/.env` and set `ANTHROPIC_API_KEY` (get one at
   [console.anthropic.com](https://console.anthropic.com)). The rest of the app
   works without it.

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

For each target, Bonsai emits files under `bindscout/outputs/<TARGET>/`,
downloadable from the page:

- `trimmed.pdb` — the design-ready trimmed target (author numbering preserved)
- `original.cif` — the chosen assembly (or AlphaFold model)
- `summary.json` — structured result: chosen structure + reasoning, topology, ECD
  ranges, hotspots, epitope patch, avoid set, warnings
- `hotspots.csv` — ranked hotspots with an `in_patch` flag
- an epitope-patch config

---

## Ask Bonsai (chat assistant)

A chat panel in the web app that answers questions about the target you've loaded
and can re-run the pipeline in plain English. Claude orchestrates the **same
deterministic tools** the app already uses, so answers are grounded in real
pipeline output — not guessed. Requires `ANTHROPIC_API_KEY` (see
[Running Bonsai](#running-bonsai)); the rest of the app works without it.

**Try asking:**

- *"Why was this PDB chosen?"* · *"What were the runner-up structures?"*
- *"List the hotspots and the epitope patch."* · *"What's in the avoid set, and why?"*
- *"What are the extracellular-domain ranges?"* · *"Which residues are membrane-proximal?"*
- *"What are the glycosylation sites on EGFR?"* · *"Show the candidate structures."*
- *"Re-trim with a 20-residue membrane buffer."* · *"Run it again without an epitope patch."*
- *"Prefer an apo structure instead."* · *"Prepare insulin instead."* (the viewers refresh)
- *"What does EGFR do?"* — short, general-background answers too.

---

## How Bonsai picks a structure

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
<summary><b>The structural details Bonsai handles correctly</b></summary>

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

## Next steps

- **By-domain / residue-range trimming** — let the assistant trim to a specific
  Pfam domain or residue range (e.g. "trim to just domain 3"), not only the full ECD.
- **Grounded general knowledge** — give the assistant web search so background
  questions cite live sources instead of relying on model memory.
- **Stability & affinity scoring** — integrate a folding/ΔΔG and binding-affinity
  tool so the assistant can answer "will this stay folded?" with a real number
  rather than a caveat.

---

## Development

```bash
cd bindscout
uv sync --extra dev
uv run pytest -m "not network"     # offline unit tests (run in CI)
uv run pytest                       # full suite (hits UniProt/PDBe/RCSB)
```

See [`CONTRIBUTING.md`](CONTRIBUTING.md) for how to contribute.

---

## License

[MIT](LICENSE) © Bonsai contributors.
