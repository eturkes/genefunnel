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

2026-07-14 | session 4 | commit `HEAD`
Scope: representation-preserving dense/sparse scoring dispatch.
Changed: base + `Matrix` dense inputs stay dense through bounded column work;
sparse inputs stay sparse and normalize only bounded chunks; separate native
wrappers share one numerical core; bounded semantic-column validation avoids
dense intermediates and raw structured-storage false positives; dispatch and
dense/sparse equivalence regressions cover zeros, partial coverage, missingness,
and reordered dimensions.
Verified: testthat = 131/131 expectations; regenerated Rcpp/roxygen outputs;
source build + built-tarball `R CMD check --no-manual` = `Status: OK`.
Decisions: storage kind is fixed from the validated top-level input and guarded
at each native boundary; implicit sparse zeros remain random-access observations;
validation reads represented values because raw structured `Matrix` storage can
contain ignored entries or unaggregated duplicates.
Remaining: Session 5 - bounded parallel chunking, deterministic portable backend
coverage, contextual worker errors, and removal of unused OpenMP flags.
Risks/blockers: one BiocParallel task per sample still captures the full input;
sparse random access is semantically correct but remains a profiling candidate.

2026-07-14 | session 5 | commit `HEAD`
Scope: deterministic, worker-bounded BiocParallel execution.
Changed: manager-side `bpiterate()` streams contiguous matrix slices through at
most two tasks per worker (one serial task); namespace worker function captures
no full matrix; membership indices are static worker arguments; explicit task
IDs/ranges restore input order; failures add global chunk/sample context; unused
OpenMP flags removed.
Profile: exploratory local 10,000 x 480, 800 sets x 24, started two-worker SOCK,
minimum of two elapsed runs. Dense 36.75 MiB input: 0.468 s column calls ->
0.319 s four chunks; 3%-dense sparse 1.78 MiB input: 0.450 s -> 0.405 s. This
selects two chunks/worker for limited load balancing; results are context, not a
package performance claim. Baseline `bplapply()` grouped column evaluations into
backend jobs by parameter settings but still invoked one native call per column
and serialized its full-matrix closure to workers.
Verified: testthat = 162/162 expectations, including dense/sparse serial-SOCK
identity, reordered completion assembly, one sample with two workers, live
backend reuse, nonempty worker chunks, and contextual errors; source build +
built-tarball `R CMD check --no-manual` = `Status: OK`; install compile lines
contain no OpenMP flags.
Decisions: bounded task count takes precedence over fixed-width buffers; iterator
dispatch limits live slices to backend capacity while preserving a small task
count for very wide matrices. Full benchmark/memory characterization remains
Session 7 work.
Remaining: Session 6 - randomized oracle/property/invariance suite plus native
memory and undefined-behavior tooling where practical.
Risks/blockers: sparse native random access remains the known profiling target;
no release blocker identified.
