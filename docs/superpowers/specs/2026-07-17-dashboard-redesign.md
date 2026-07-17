# Dashboard Redesign — Report/Paper Style

**Date:** 2026-07-17
**Status:** Approved design. Supersedes the interactive-explorer parts of
`2026-07-16-dashboard-design.md`.

---

## 1. Why

The first dashboard was an interactive explorer. In use it read as clunky and
"broken": the Part 3 tab exposed three independent filter selectboxes, and **10
of the 18 condition × treatment × sample-type combinations are structurally
empty** (every `healthy` condition and every `none` treatment has no
responder/non-responder label). Users hit dead combos and assumed a bug. The
methods were described in the sidebar but the *results* were not shown
prominently.

**There is no cohort bug.** The required cohort (melanoma + miraclib + PBMC) has
656 subjects (331 responders / 325 non-responders); all 8 meaningful
combinations are large. The interactive filtering was obscuring correct results.

## 2. Decision

Pivot to a **report/paper style**: each tab presents its instruction part's
required analysis as a clean report, fixed on the required filters, with results
front and center. No filter selectboxes that can produce empty cohorts.

## 3. Structure — 3 tabs, one per instruction part

Minimal sidebar: app title + a one-line dataset description only. **All methods
AND results live in the tabs** (fixing "methods on the sidebar, no results").

- **Part 2 · Cell-Type Frequencies**
- **Part 3 · Statistical Analysis**
- **Part 4 · Subset Analysis**

### Part 2 — Summary table (plain)
Framing sentence; a KPI row (samples / subjects / projects); a mean-composition
chart (mean % per population); then the full per-sample frequency table
(`sample, total_count, population, count, percentage`), searchable/scrollable.
No filter selectboxes.

### Part 3 — Statistical Analysis as a paper (fixed: melanoma + miraclib + PBMC)
No filters. Structured report:
1. **Question** — do population frequencies differ between miraclib responders
   and non-responders?
2. **Cohort** — 656 subjects (331 R / 325 NR), 3 timepoints, PBMC.
3. **Methods** — mixed-effects model (primary), Mann-Whitney U cross-check,
   Benjamini-Hochberg FDR, bootstrap CIs; one-line rationale each.
4. **Results** — the plain-language finding (prominent, via `finding_summary`);
   boxplots for all five populations; and the **full statistics table always
   visible** (not in an expander) with paper-quality column names: Δ
   (responder−non, pp), 95% CI, mixed-model p, FDR q, MWU q, effect size,
   significant. cd4_t_cell highlighted.
5. **Interpretation** — cd4_t_cell is a modest candidate (+0.64 pp, q=0.025),
   non-concordant with the conservative test (q=0.062); no other population
   differs; near-zero random-effect variance noted.

### Part 4 — Subset Analysis, degenerate cases surfaced
- **Required baseline subset** (melanoma / miraclib / PBMC, t=0): breakdowns by
  project, response, and sex (tables + small bars), and the **posed answer
  10206.15** as a prominent metric.
- **Data-coverage callout:** a table of baseline (t=0, PBMC) subset sizes across
  condition × treatment, explicitly flagging the structurally **degenerate**
  ones (every `healthy` condition and every `none` treatment → 0
  responder-labelled samples). Makes the degeneracy explicit rather than a
  stumble.

## 4. Removed / changed

- Delete the three independent Part 3 filter selectboxes (source of the
  dead-combo confusion). Part 3 is fixed.
- Sidebar shrinks to title + one line; the methods blurb moves into Part 3's
  Methods section next to the results.
- Rename section files to their parts: `overview.py → part2.py`,
  `response.py → part3.py`, `subset.py → part4.py`; update `app.py` imports and
  tab labels.
- `components.py` gains two pure, tested helpers:
  - `format_stats_table(comparison) -> DataFrame` — the paper-quality display
    table (readable column names, rounded).
  - `subset_coverage(metadata) -> DataFrame` — baseline subset sizes per
    condition × treatment with a `degenerate` flag.
- Keep `data.py` loaders; Part 3 calls them with the fixed required args. Keep
  `finding_summary`, `ensure_database`, the theme, and the self-building DB.

## 5. Testing

- Unit: `format_stats_table` (columns/rounding on a `compare_responders`-shaped
  df) and `subset_coverage` (degenerate flag true for healthy/none, false for
  melanoma/miraclib) — pure functions.
- Keep `ensure_database` / `finding_summary` / cross-thread tests.
- Update the `AppTest` smoke test: asserts no exception and that the Part 3
  finding renders (fixed cohort always produces a result now).

## 6. Out of scope

No changes to `analysis/`, the pipeline, the DB, or the deployment mechanism
(self-building DB, Streamlit Cloud). Python-3.11 deployment note handled
separately.

## 7. Next steps

writing-plans → subagent-driven implementation: (1) components helpers; (2) app
shell + Part 2 + Part 4 tabs + smoke; (3) Part 3 paper.
