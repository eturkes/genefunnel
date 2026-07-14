# GeneFunnel execution ledger

Normative specification + completion checkboxes = `PLAN.md`.

## Session 1 baseline - 2026-07-14

- Git: `HEAD 3fd2c23`; `main` clean + one commit ahead of `origin/main` before
  baseline commands.
- Tracked tree: package metadata/namespace; one R wrapper + generated Rcpp R
  binding; one C++ scorer + generated registration + Makevars; one Rd page;
  README/GPL; four excluded bootstrap scripts; RStudio project; AGENTS/PLAN.
  No tests, vignette, CI, release, or benchmark harness. Three compiled objects
  were tracked under `src/`.
- Toolchain: Debian 13; R 4.6.1; Bioconductor 3.23; GCC/G++ 14.2.0; GNU Make
  4.4.1. Project-local packages: Matrix 1.7-5, BiocParallel 1.46.0, Rcpp 1.1.2,
  RcppArmadillo 15.4.0-1, roxygen2 8.0.0, testthat 3.3.2.
- Untouched baseline: `R CMD build .` passed. First `R CMD check --no-manual`
  stopped with one dependency error because BiocParallel/Rcpp/RcppArmadillo
  were absent. After project-local installation, tarball check completed with
  `1 WARNING` (placeholder license) + `1 NOTE` (unqualified `bpparam`).

## Handoff ledger

2026-07-14 | session 1 | commit `HEAD` (`package: establish trustworthy build and test baseline`)
Scope: trustworthy source-only build + installed-package smoke-test baseline.
Changed: ignored/removed native build products; truthful minimum DESCRIPTION;
testthat edition 3 harness; pure-R formula oracle; dense smoke fixtures.
Verified: local testthat = 12 expectations; source tarball manifest excludes local
artifacts; `R CMD build .` passes; installed-tarball `R CMD check --no-manual`
passes with 0 errors, 0 warnings, 1 known note.
Decisions: testthat follows current official Bioconductor preference; scoring R/C++
code stays unchanged in this session; development dependencies stay repo-local.
Remaining: Session 2 - validation, shared membership/coverage preparation,
partial-coverage policy, namespace-qualified backend default.
Risks/blockers: unqualified `bpparam` causes the sole check note and leaves the
implicit default unusable until Session 2; explicit `SerialParam()` is covered.

2026-07-14 | session 2 | commit `HEAD` (`api: lock validation and coverage contract`)
Scope: strict R-level inputs + shared set preparation/coverage + partial-set policy.
Changed: exact base numeric/integer and numeric `Matrix` class boundary; dimensions,
identifiers, finite/non-negative values, and `BiocParallelParam` validation;
`gene_set_coverage()` diagnostics; stable deduplication/matching/indices; aggregate
omission warning; partial scoring; namespace-qualified registered backend default;
native scorer now consumes prepared integer memberships.
Verified: testthat = 96/96 expectations; dense/integer/Matrix + sparse validation
fixtures; generated docs/bindings; source manifest; built-tarball
`R CMD check --no-manual` = `Status: OK`.
Decisions: zero-member declared sets are valid diagnostics (`coverage = NA`) but
unscoreable; zero retained sets return a warned 0-by-sample numeric matrix; native
membership refactor moved forward from Session 3 to eliminate per-sample rematching.
Remaining: Session 3 - sample-specific missing omission + effective size; defensive
native infinity/result checks; scale-aware zero tolerance; oracle fixtures.
Risks/blockers: dense inputs still route through sparse storage; scoring still uses
one task per sample and unused OpenMP flags. Both remain explicitly scheduled.

2026-07-14 | session 3 | commit `HEAD`
Scope: sample-specific missing-value semantics + defensive native numerics.
Changed: native cells omit `NA`/`NaN`, retain zeros, use effective observed size,
return `NA` below two observations, preserve small positive scores, and reject
escaped infinities/non-finite or materially negative arithmetic with context;
canonical dense/sparse/oracle regressions cover the contract.
Verified: testthat = 112/112 expectations; regenerated Rcpp attributes unchanged;
source build + built-tarball `R CMD check --no-manual` = `Status: OK`.
Decisions: an algebraically equivalent below-mean form avoids catastrophic
subtraction while preserving the normative equation; scale-aware tolerance clamps
negative roundoff only, never positive scores; native guards remain defensive
behind strict R validation.
Remaining: Session 4 - distinct dense/sparse paths, representation preservation,
stored sparse missing equivalence, and a testable dispatch seam.
Risks/blockers: dense inputs still route through sparse storage; scoring still uses
one task per sample and unused OpenMP flags. All remain explicitly scheduled.
