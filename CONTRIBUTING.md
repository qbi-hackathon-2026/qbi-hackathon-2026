# Contributing to BindScout

Thanks for your interest in improving BindScout! Bug reports, fixes, and small
focused improvements are all welcome.

## Ways to contribute

- **Report a bug** — open a GitHub Issue describing what you did, what you
  expected, and what happened. Include the protein/accession you searched and
  any error text from the page or server logs.
- **Suggest a fix** — open a Pull Request (see below).
- **Propose a feature** — open an Issue first so we can discuss scope before you
  invest time.

## Development setup

The app is a deterministic Python pipeline (`bindscout/src/bindscout/`, managed
with [uv](https://docs.astral.sh/uv/)) served by a FastAPI app, plus a single
static frontend file (`bindscout/frontend/index.html`). There is no build step
for the frontend.

```bash
cd bindscout
uv sync                                 # creates an isolated env (Python >=3.11)
uv run python -m bindscout.server        # then open http://127.0.0.1:8000/
```

`uv` downloads the right Python automatically, so you don't need a matching
system interpreter. See the [root README](README.md) for the interfaces (web app
/ CLI / MCP).

## Running tests

```bash
cd bindscout
uv sync --extra dev
uv run pytest -m "not network"   # offline unit tests (what CI runs)
uv run pytest                     # full suite (hits UniProt/PDBe/RCSB)
```

CI runs the `not network` subset on every pull request, so make sure that
passes locally before opening a PR. Tests marked `network` require live
internet access to external services (UniProt/PDBe/RCSB/EBI) and may be flaky in
CI.

## Pull request guidelines

- **Branch from `main`** and open the PR against `main`.
- **Keep PRs focused** — one logical change per PR is much easier to review.
- **Match the surrounding style.** The backend is plain typed Python; the
  frontend is a single vanilla-JS/CSS file with no framework or build step.
- **Don't commit secrets, credentials, or generated output.** The pipeline
  writes regenerated files under `bindscout/outputs/` (and caches under
  `bindscout/cache/`) — these are gitignored and must not be added.
- **Describe what and why** in the PR body, and link any related Issue.

## A note on CI and deployment for outside contributors

For security, the deployment workflow only runs **after a maintainer merges to
`main`** — it never runs on pull requests, and it does **not** expose any cloud
credentials to PR builds. Pull requests from forks run only the read-only
checks (tests/lint) with no secrets. This is expected: your PR's CI will not
deploy anything, and that's by design.

## Code of conduct

Be respectful and constructive. Assume good intent, keep discussion focused on
the work, and help make this a welcoming project for researchers and developers
alike.
