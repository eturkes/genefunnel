# GeneFunnel publication acceptance

Fixed review set → adjudicate against this change; one finding-specific re-review.

| ID | Acceptance target | Evidence entry point |
|---|---|---|
| G01 | Exact Pietra identity; score bounds, monotonicity, concavity, plateau; correct novelty scope | `article.Rmd`; `results/algebra.csv` |
| G02 | Explicit estimands, faithful published comparators, fair preprocessing, paired simulation uncertainty | `analyze.R`; `results/simulation*.csv` |
| G03 | Coverage, set size, nuisance scale, and preprocessing stress expose failure regimes | `results/coverage.csv`; `results/simulation*.csv`; `results/ad*.csv` |
| G04 | Public data provenance, immutable fixture identities, analysis-unit and leakage controls | `results/input-manifest.csv`; article Methods |
| G05 | Independent perturbation evaluation; no invented biological ground truth | `results/airway*.csv`; article Results |
| G06 | Disease analysis uses subject-level uncertainty; no causal or diagnostic validation claims | `results/ad*.csv`; article Discussion |
| G07 | Dense/sparse agreement and honest runtime/memory measurement | `results/runtime.csv`; `analyze.R` |
| G08 | Clean-source reproducibility; deterministic outputs, seeds, dependency inventory | `reproduce.R`; `results/replay.csv`; `results/environment.txt` |
| G09 | Appropriate package checks; generated interfaces unchanged | `tools/check-package.R --bioc` |
| G10 | Separate readable manuscript, figures, supplement, claims tied to evidence; external blockers explicit | `article.Rmd`; `article.html`; `REVIEW.md` |
