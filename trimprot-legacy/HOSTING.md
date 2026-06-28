# TrimProt hosting — work log & remaining steps

> **Status: TABLED.** This is a local planning doc (gitignored, not published).
> It captures where the hosting effort stands so it can be resumed later.
>
> Decisions locked in: **Google Cloud Run** (scale-to-zero), **cost minimized**
> (paid personally for now, cold starts acceptable), **shared cross-user cache**
> as a core feature, **Workload Identity Federation** for CI auth (no stored
> keys), public GitHub repo open to outside contributors.

## Architecture (target)

```
  Browser ──> Frontend (static index.html, served by the backend at /)
                 │   GET /api/search, /api/run, /api/download/*
                 ▼
          FastAPI backend — ONE Cloud Run service (min-instances=0)
                 │
        ┌────────┴─────────┐
        ▼                  ▼
   Cache lookup        Object storage (GCS bucket)   ← Phase 2
   accession → result   trimmed PDB / CIF / AF outputs,
   JSON + filenames     written once, served to all, survive restarts
```

The shared cache is the main cost lever: each protein is computed once, ever,
for all users; subsequent lookups are near-free cache hits. The expensive part
is a *cold* pipeline run (CIF/AlphaFold download + gemmi compute), not idle time.

## What's DONE

- **Phase 0 — repo prep (merged to `main`, PR #3):**
  - `.gitignore` rewritten for the Python backend + static frontend; defensively
    ignores any GCP key files, plus `.claude/`, `CLAUDE.md`, and this `HOSTING.md`.
  - `LICENSE` (MIT, "TrimProt contributors" — independent team, not UCSF).
  - `CONTRIBUTING.md`, PR template.
  - Read-only CI workflow (`.github/workflows/ci.yml`): runs `pytest -m "not
    network"` on PRs (incl. forks) and pushes to main; **no secrets**,
    read-only permissions — safe for untrusted PRs. (Verified: 68 pass.)
  - `requirements.txt` / `requirements-dev.txt` as the single dependency source.

- **Tooling:** `gcloud` CLI installed locally via Homebrew (574.x) at
  `/opt/homebrew/share/google-cloud-sdk/bin`. Not yet on the default shell PATH.

- **Phase 1 code — committed on branch `cloud-run-deploy` (NOT yet merged):**
  - `trimprot/Dockerfile` (at build-context root so the image includes the
    backend *and* the `frontend/` it serves at `/`); starts uvicorn honoring
    Cloud Run's `$PORT`. Start command validated locally (serves `/` and
    `/api/search`, HTTP 200) without Docker.
  - `trimprot/.dockerignore` (keeps the image lean / secret-free).
  - `.github/workflows/deploy.yml`: builds via Cloud Build and deploys to Cloud
    Run **on push to `main` only** (never PRs/forks), auth via **WIF**
    (`id-token: write`, no stored key); `min-instances=0`, `max-instances=3`,
    1Gi / 1 CPU, 300s timeout. Config comes from GitHub **repo variables** (not
    secrets): `GCP_PROJECT_ID`, `GCP_REGION`, `GCP_WIF_PROVIDER`,
    `GCP_DEPLOY_SA`, `CLOUD_RUN_SERVICE`.

## What's LEFT (to resume Phase 1 → live URL)

All remaining Phase-1 steps need a GCP project and are run by a human (they touch
billing / create cloud resources). Sequence:

1. **GCP project setup (you):**
   - `gcloud auth login`
   - Create/choose a project; **enable billing** (via console — needs a billing
     account). Note the **project ID** and pick a **region** (default suggestion:
     `us-central1`).
   - Set a **budget alert** (~$5–10) as the cost guardrail.
2. **Enable APIs + service account (commands to be generated):** Cloud Run, Cloud
   Build, Artifact Registry; create the runtime/deploy service account with
   minimal roles.
3. **Workload Identity Federation (commands to be generated):** create a pool +
   provider trusting **only** `qbi-hackathon-2026/qbi-hackathon-2026` on `main`;
   grant the SA impersonation from that provider.
4. **GitHub repo variables:** add the 5 `vars.*` listed above (Settings → Secrets
   and variables → Actions → **Variables** tab). No secrets needed (WIF).
5. **Ship:** open PR for `cloud-run-deploy` → merge → push-to-main triggers
   `deploy.yml` → live Cloud Run URL. Verify with `curl /api/run?...` and the page.
   - ⚠️ Don't merge `cloud-run-deploy` until steps 1–4 are done, or the first
     deploy run will fail for lack of GCP config (harmless, but noisy).

## Phase 2 (after Phase 1 is live): durable shared cache + GCS outputs

- Today `api.py`'s `_cache` is an in-memory dict — empty on restart, per-instance.
  `run_egfr.py` writes outputs to `trimprot/backend/output/` (incl. AlphaFold
  `AF-*` files), which is ephemeral on Cloud Run. Both must move to GCS so the
  cache is shared + persistent and `/api/download/*` survives redeploys.
- Cheapest primitives: cache as a GCS object per accession
  (`gs://<bucket>/cache/<acc>.json`); outputs under `gs://<bucket>/outputs/`;
  lifecycle rule to expire large files. (Graduate to Firestore only if needed;
  avoid always-on Cloud SQL.)
- New `storage.py` behind `api.py`; the pipeline modules stay unchanged.

## Phase 3 (polish)

Pre-warm common targets (EGFR `P00533`, TNF-α `P01375`), a "popular/recent
proteins" list off the cache index, optional cache-stats endpoint.

## Cost guardrails (set day one — personal account)

- Cloud Run `min-instances=0`, small `max-instances` cap (in `deploy.yml`).
- GCP **budget alert** ($5–10).
- GCS lifecycle rule to expire large output objects (Phase 2).

## Notes / gotchas

- gcloud not on default PATH; new shells may need the Homebrew SDK path added, or
  prefix `PATH=/opt/homebrew/share/google-cloud-sdk/bin:$PATH`.
- Frontend could later move to Vercel for per-branch preview URLs (point its
  `API_BASE` at the Cloud Run URL) — but do **not** try to host the stateful
  backend on Vercel serverless; it fights the filesystem + cache model.
