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

2026-07-14 | session 6 | commit `HEAD`
Scope: complete scientific invariant + adversarial regression suite.
Changed: matrix-level pure-R oracle; 32 seeded randomized dense/sparse fixtures;
direct score-bound, permutation, independence, overlap, scale/shift, minimum-size,
extreme/subnormal, signed-zero, integer, structural-zero, output-shape, and
immutability coverage; oracle handles subnormal mean underflow explicitly.
Verified: testthat = 42 tests/342 expectations; installed-package suite passes
under UBSan and combined ASan/UBSan; source build + built-tarball
`R CMD check --no-manual` = `Status: OK`.
Decisions: sanitizer leak detection stays disabled because R-runtime lifetime
allocations are out of package scope; fixed-seed randomized tests calculate their
expected values independently through the R oracle.
Remaining: Session 7 - reproducible synthetic benchmark generator/harness,
profiles, sparse memory evidence, optimization only where measured.
Risks/blockers: Valgrind is unavailable locally; ASan covers native memory access.
Sparse random access remains the next profiled candidate; no release blocker.

2026-07-14 | session 7 | commit `HEAD`
Scope: reproducible performance evidence + semantics-preserving sparse optimization.
Changed: deterministic dense/sparse low/high-overlap full, hotspot, and CI-smoke
presets; isolated timing/allocation/RSS/digest/environment capture; one-pass
catalogue matching; sparse column streaming with analytical implicit zeros;
qualified README sparse claim.
Verified: 44 tests/349 expectations; built-tarball `R CMD check --no-manual` =
`Status: OK`; combined ASan/UBSan suite passes; smoke serial/SOCK digests match.
Exact 20,000 x 200 hotspot comparison against Session 6, three isolated repeats:
all digests unchanged; median speedup 5.5x-11.3x (median 9.3x). Full 20,000 x
600, 3%-stored sparse serial medians = 0.159-0.165 s; 5.81 MB input + 4.92 MB
output; peak scoring increment = 4.88 MB versus 96 MB logical dense payload.
Decisions: base-R profiling + optional GNU time avoid new dependencies; allocation
tracing uses a separate identity-checked pass; sparse traversal cost now follows
stored entries; performance results remain generated/untracked and contextual.
Remaining: Session 8 - metadata, namespace/source hygiene, maintainer scaffolding,
package docs/NEWS/citation, and cross-platform portability preparation.
Risks/blockers: cold two-worker SOCK medians (0.897-1.221 s) exceed serial at
default sizes because process startup/transfer dominates; claim no parallel
speedup. Workflow wiring for the threshold-free smoke entry point remains
Session 10 CI scope.

2026-07-14 | session 8 | commit `HEAD`
Scope: release metadata + generated-interface/source-package hygiene.
Changed: complete author/license/URL/dependency metadata; four current official
leaf `biocViews`; C++17 build contract; `sourceCpp` replaced by the minimal
Rcpp namespace-loading import; generated namespace/registration refreshed;
package help, NEWS, package + thesis citations; obsolete bootstrap scripts
removed; build/ignore rules cover local and platform-native artifacts.
Verified: official R/Bioconductor metadata, licensing, namespace, citation,
NEWS, and C++ guidance reviewed; devel vocabulary confirms all four views are
leaf terms; citation/NEWS/registration/installed-score smoke passes; source
manifest contains only intentional package files; built-tarball `R CMD check
--no-manual` = `Status: OK` with all 44 tests; Linux GCC C++17, Makevars, and
GNU-extension portability checks pass.
Decisions: retain GPL-3-or-later to match source notices; `Rcpp::evalCpp` import
is required to load C callables used by registered Rcpp wrappers; pin C++17 to
avoid compiler-default drift; actual macOS/Windows execution evidence belongs
with Session 10's cross-platform CI rather than a Debian-only claim.
Remaining: Session 9 durable scientific specification, complete function docs,
vignette, and README; Session 8 macOS/Windows verification through Session 10.
Risks/blockers: empirical macOS/Windows build evidence awaits CI; no local
release blocker.

2026-07-14 | session 9 | commit `HEAD`
Scope: durable scientific contract + complete task-oriented user documentation.
Changed: installed `SCIENTIFIC_SPEC.md`; exact scorer/coverage help + examples;
executable BiocStyle vignette; concise README; vignette metadata/dependencies;
documentation NEWS; regenerated Rd files.
Verified: current official Bioconductor vignette guidance reviewed; 44 tests/349
expectations pass; roxygen/Rd checks pass; standalone HTML has zero external
resources; Chromium PDF inspection confirms formula/tables/layout; source build
executes the vignette; built-tarball `R CMD check --no-manual` = `Status: OK`,
including examples, tests, installed docs, and vignette rebuild.
Decisions: HTML math uses native MathML for deterministic offline rendering;
coverage thresholds remain labelled caller examples; README/vignette state
development status and all scientific comparability/downstream limitations.
Remaining: Session 10 - cross-platform CI/current Bioconductor gates, including
Session 8's empirical macOS/Windows portability evidence; then Sessions 11-12.
Risks/blockers: actual macOS/Windows execution awaits CI; no local blocker.

2026-07-14 | session 10 | commit `HEAD`
Scope: current Bioconductor CI + release-quality automation.
Changed: live official release/devel matrix on Linux plus devel macOS/Windows;
SHA-pinned actions and version-separated dependency caches; source-tarball build,
installed tests, packagebuilder `R CMD check`, Git-clone/tarball BiocCheck, clean
target-library install, offline Rd/link/vignette checks, sparse-memory/backend
identity smoke, and Linux ASan/UBSan jobs; RStudio project removed for repository
hygiene; condition construction and iterator state satisfy BiocCheck coding gates.
Verified locally: dynamic matrix resolves Bioconductor 3.23/3.24 with R 4.6.0;
actionlint + ShellCheck clean; clean tarball install runs 44 tests/349 expectations;
official packagebuilder environment gives `R CMD check` `Status: OK`; BiocCheck
Git-clone hygiene = 0/0/0; tarball = 0 errors/0 warnings/5 explained notes; docs
and benchmark smoke pass; serial/SOCK
digests agree; sparse fixtures stay smaller than logical dense payloads; combined
ASan/UBSan installed-package suite passes.
Decisions: CI fails every package-controlled BiocCheck error/warning; development
version numbering and maintainer-owned Support-site enrollment remain explicit
Session 12 gates. Suggested `Coverage` view is a semantic false positive; ORCID,
funding roles, and mailing-list state remain unasserted. Valgrind stays omitted
because it is unavailable locally; the already-stable stronger sanitizer pair is
gating. Current AI provenance policy is encoded for the release audit and new CI
code carries provenance comments.
Remaining: first remote workflow run must establish macOS/Windows portability and
close Session 8/Definition-of-Done cross-platform boxes; then Session 11.
Risks/blockers: GitHub-hosted results cannot exist before maintainer push; full
submission BiocCheck still needs a 0.99-series version and Support-site state.

2026-07-15 | session 10 follow-up | commit `HEAD`
Scope: reconcile and repair the first remote cross-platform CI run.
Changed: dense/sparse kernels classify zeros correctly when a positive subnormal
mean underflows on 64-bit `long double`; direct regression locks the exact cell;
clean-target CI propagates its temporary library to fresh SOCK workers and asserts
their loaded package paths.
Verified: first remote run proved Linux release/devel, Windows devel, and sanitizer
jobs green plus macOS ARM64 compile/install portable; local forced
`-mlong-double-64` suite = 44 tests/350 expectations; clean-tarball SOCK suite,
`R CMD check` (`Status: OK`), documentation, benchmark smoke, and combined
ASan/UBSan pass; BiocCheck Git clone = 0/0/0 and tarball = 0 errors/0 warnings/5
explained notes.
Decisions: positive totals determine whether exact sparse/dense zeros lie below an
underflowed mean; clean-process library propagation belongs to the CI harness,
while workers still prove they execute the just-built tarball installation.
Remaining: maintainer push + remote rerun must make macOS tests/check and the
quality job green; then close Session 8/10 + Definition-of-Done platform boxes and
begin Session 11.
Risks/blockers: post-fix macOS and clean-runner evidence requires the next remote
workflow; no local release blocker remains.

2026-07-15 | session 11 | commit `HEAD`
Scope: versioned synthetic scientific + computational benchmark protocol.
Changed: protocol 1.0.0 manifest; nine controlled scenarios with analytic/paired
assertions; shared environment capture + plain Markdown report generation;
performance report; combined CI smoke; complete provenance/reproduction boundary.
Verified: controlled = 9 scenarios/14 assertions; smoke = eight serial/SOCK cases
with paired digests + controlled suite; package = 44 tests/350 expectations;
documentation checks pass; clean source build + `R CMD check --no-manual` =
`Status: OK`; tarball BiocCheck = 0 errors/0 warnings/5 explained notes.
Decisions: benchmark semantics carry an explicit protocol version; competitor
methods stay absent until a comparative claim justifies pinned APIs/dependencies;
future real data stay external with URL/license/checksum/immutable provenance.
Remaining: maintainer push + green remote rerun closes Session 8/10 platform/CI
gates; Session 12 release-candidate audit follows.
Risks/blockers: remote evidence for the unpushed portability fix remains the sole
pre-audit gate; no local release blocker.
