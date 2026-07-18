# GeneFunnel execution ledger

Scientific contract = `inst/SCIENTIFIC_SPEC.md`; forward research portfolio =
`PLAN.md`; this file = chronological execution evidence.

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

2026-07-15 | session 12 audit 1 | commit `HEAD`
Scope: finite-overflow numerical portability audit + repair.
Changed: bounded-coefficient evaluation avoids weighted-term overflow; a
max-scaled fallback handles overflowing finite sums while retaining errors for
unrepresentable outputs; the independent R oracle and dense/sparse regressions
cover both paths.
Verified: normal and forced 64-bit-`long double` suites pass; forced coverage
includes repeated maxima, maximum-plus-half, equal quarters, and prior subnormal
cases.
Decisions: scaling is an overflow-only numerical representation and leaves the
equation/API unchanged; ordinary arithmetic order stays stable when intermediates
fit.
Remaining: complete local Session 12 contract/release audit and full benchmarks;
remote rerun still gates macOS/Windows/CI closure and submission numbering.
Risks/blockers: post-fix remote evidence still requires maintainer push; no local
numerical blocker remains.

2026-07-15 | session 12 audit 2 | commit `HEAD`
Scope: local release-candidate trace, reproducibility rerun, and provenance gate.
Changed: current-policy `Assisted-by: OpenAI Codex` citations now cover all
contributed non-generated code; check-log assertion ignores nested scratch trees;
README benchmark link survives the source-tarball boundary; local audit boxes and
handoff state reconciled.
Verified: every scientific requirement traced to implementation/tests/docs; only
the intended two functions export; generated native/roxygen files are current;
fresh committed clone gives BiocCheckGitClone 0/0/0, 260 KiB source tarball,
`R CMD check` = `Status: OK`, and clean-library tarball tests on fresh SOCK
workers; offline docs, benchmark smoke, and combined ASan/UBSan pass; tarball
BiocCheck = 0 errors/0 warnings/5 justified notes; package name absent from
current/archive CRAN and current/removed Bioconductor indices. Clean-HEAD
protocol 1.0.0 archive = 14/14 controlled assertions + 24/24 stable full-run
digests; sparse serial peak increment ~4.9 MiB vs 96 MB logical dense payload.
Decisions: keep `0.0.0.9000` until remote platform gates pass and submission is
authorized; keep advisory BiocCheck notes rather than invent Coverage semantics,
ORCID, or funding; parallel speedup remains unclaimed because SOCK overhead wins
at default benchmark sizes.
Remaining: maintainer push + green macOS/Windows/quality rerun; then enable full
submission checks, use `0.99.0`, disclose AI assistance in the issue/PR, and tag
only after review.
Risks/blockers: all remaining gates require remote evidence or maintainer-owned
submission state; local release audit has no known package-controlled blocker.

2026-07-15 | session 12 audit 3 | commit `HEAD`
Scope: remote portability closure + Bioconductor submission-version promotion.
Changed: green cross-platform/CI/Definition-of-Done gates reconciled; package +
NEWS promoted to `0.99.0`; BiocCheck version validation enabled; Support watched
tag retained as an explicit maintainer-owned gate.
Verified: GitHub Actions run `29385049469` at `93e2e4c` passed all seven jobs,
including Linux release/devel, macOS ARM64 devel, Windows devel, clean-tarball
quality, and ASan/UBSan; local `0.99.0` tarball clean-library tests pass; `R CMD
check --no-manual` = `Status: OK`; tarball BiocCheck with version validation =
0 errors/0 warnings/5 justified notes. Current official guidance still mandates
`0.99.0` and issue + PR AI provenance disclosure.
Decisions: Support account exists, while `genefunnel` is absent from its watched
tags; retain only `no-check-bioc-help` until the maintainer changes that external
profile state. Submission issue/PR and release/tag stay maintainer-owned.
Remaining: add the Support watched tag → remove its check suppression; push this
commit + confirm the exact candidate CI; disclose `Assisted-by: OpenAI Codex` in
the Bioconductor issue and PR; tag/release only after review.
Risks/blockers: watched-tag, submission, push, and review state are external;
all locally actionable release gates pass.

2026-07-15 | post-plan audit 1 | commit `HEAD`
Scope: prevent S3 class spoofing across public input boundaries.
Changed: accept plain base matrices/lists/character identifiers or genuine S4
matrix/backend objects only; shared identifier validation retains member
duplicates solely for canonicalization; docs + NEWS state the tightened contract.
Verified: red exploit changed `c(10, 0)` from the correct score `0` to `20` via
a classed member's `[` method; 47 tests/373 expectations pass after rejection;
clean source-tarball install/tests including fresh SOCK workers pass; `R CMD
check --no-manual` = `Status: OK`; fresh-clone BiocCheck = 0/0/0 and tarball
BiocCheck = 0 errors/0 warnings/5 justified notes.
Decisions: validation and scoring must observe one representation; S3 subclasses
remain caller-convertible but are not implicit package inputs. The Support
watched-tag gate remains suppressed/deferred after the correct account/profile
continued to return an empty API tag list; it no longer blocks other work.
Remaining: choose the next post-plan feature or maintenance objective; resume
submission issue/PR disclosure and review/tag gates when desired.
Risks/blockers: this commit needs the normal maintainer push + cross-platform CI;
no local package-controlled blocker is known.

2026-07-16 | submission gate 1 | commit `HEAD`
Scope: enable the Bioconductor Support watched-tag check.
Changed: removed the final BiocCheck suppression after the maintainer added
`genefunnel` to watched tags; reconciled the Session 12 checklist.
Verified: live Support APIs find the maintainer account and return watched tags
`genefunnel,imputefinder`; unsuppressed source-tarball BiocCheck reports the
maintainer registered and package present in watched tags.
Remaining: disclose `Assisted-by: OpenAI Codex` in the submission issue and PR;
create a tag/release only after maintainer review.
Risks/blockers: submission and review state remain maintainer/external; no local
package-controlled blocker is known.

2026-07-17 | frontier plan 1 | commit `HEAD`
Scope: replace the completed release roadmap with a forward research +
scale portfolio.
Changed: authority now separates durable scientific contract, future plan, and
execution evidence; new workstreams cover exact score decomposition,
aggregation gap, compiled provenance, block execution, reliability, and held-out
validation with explicit kill gates.
Verified: live tracked architecture/status/history; cited official Bioconductor
DelayedArray/SparseArray/HDF5Array/container surfaces; cited primary scoring +
pathway benchmark literature; randomized decomposition/Jensen identities;
Markdown render + whitespace checks.
Decisions: freeze the existing scorer; report all-zero balance as undefined;
treat aggregation gap as a hypothesis rather than biology; make scale adapters
explicit and Bioconductor-native; gate public API/dependency work on candidate
release state and measured value.
Remaining: close/defer release lane, then formalize decomposition semantics and
protocol 2.0 design before implementation.
Risks/blockers: public API/dependency expansion may conflict with the
pre-submission candidate; synthetic fidelity is not biological validation;
research spikes may proceed while public surface changes await the candidate
scope decision.

2026-07-17 | frontier execution 1 | commit `HEAD`
Scope: reconcile the first release-lane gate at the frozen candidate baseline.
Changed: `PLAN.md` records exact post-S3-fix cross-platform CI evidence at
`c789a23`.
Verified: GitHub Actions run `29461199816` completed successfully at exact SHA
`c789a237e579b73fdf47f89a4e8aa54282c5e17f`; all seven checks passed: Linux
release/devel, macOS devel, Windows devel, release-quality, sanitizers, and the
matrix resolver. Full local testthat suite passes; `git diff --check` is clean.
Decisions: `c789a23` remains the package candidate/innovation boundary; later
local planning/profile commits neither change package code nor reopen scope.
Remaining: maintainer chooses candidate freeze versus intentional reopening;
issue/PR AI disclosure and reviewed tag/release remain external release gates.
Risks/blockers: A1 research can proceed, but PLAN execution order blocks public
API/dependency work until the release lane is explicitly closed or deferred.

2026-07-18 | frontier execution 2 | commit `HEAD`
Scope: intentionally reopen candidate scope; complete A1 theorem + component
semantics before public/native implementation.
Changed: release lane records the maintainer decision and defers submission/tag
gates; installed `COMPONENTS_SPEC.md` proves the score factorization, fixes cell
semantics + numerical status contract, audits mathematical/ecological/pathway
prior art, and commits interpretation cases; independent share/TV R oracle +
property tests cover missingness, zeros, bounds, symmetry, homogeneity,
concavity, and Pietra/Bulla identities; README/NEWS/core spec link the evidence.
Verified: primary papers + current official singscore/GSVA APIs; full local
testthat suite; strict GFM → MathML render with all seven display equations;
source tarball contains both installed specifications + theorem tests; `R CMD
check --no-manual` = `Status: OK`; offline Rd/link/vignette checks pass;
tarball BiocCheck = 0 errors/0 warnings/6 advisory notes.
Decisions: candidate deliberately reopened; external submission/release waits
for stabilized new scope. `balance` is exactly Bulla's established evenness
index and the normalized Pietra/Hoover complement - neither index nor
normalization is novel. Retain `balance` as the arithmetic field name; restrict
claims to GeneFunnel factorization, implementation/audit value, and held-out
biology. Core score remains authoritative; A1 oracle covers ordinary values.
Remaining: A2 begins with committed safe-region/error + performance protocols
and a scaled/higher-precision oracle, then one-pass prototype/native work.
Risks/blockers: Bulla balance depends on effective support/member identity;
scaled/unavailable/ill-conditioned sidecars are specified but unimplemented;
diagnostic biological value remains an unvalidated hypothesis.

2026-07-18 | frontier execution 3 | commit `HEAD`
Scope: freeze A2 numerical/status + default-path performance evidence before
native/public implementation.
Changed: protocol `A2-1.0.0` fixes the component result schema, availability/
conditioning states, safe-region budget, extreme fixtures, four dense/sparse
performance workloads, five-call warm estimator, 30 balanced randomized pair
orders, 5% one-sided timing bound, allocation/RSS gates, and Linux quiescence
rule; exact Git-snapshot preparer fingerprints source + installed trees;
dependency-free scaled double-double R oracle covers ~2,000-bit exponent spans.
Verified: 57 tests/1,888 expectations; ordinary + overflow/cancellation/
underflow/subnormal identities; prepared-library smoke = exact digests,
environment identity, allocation/RSS capture, and decision pass; malformed/
outside-repo preparation rejected; documentation checks pass; source build +
`R CMD check --no-manual` = `Status: OK`; tarball BiocCheck = 0 errors/0
warnings/6 advisory notes.
Calibration: an identical-binary 30-pair first-call run failed three of four 5%
bounds (upper ratios 1.091-1.144) while host load reached 10-27 on 8 CPUs; an
earlier design also let the RSS sampler contend with timing. These negative
controls invalidated cold timing as the primary endpoint, isolated passive RSS,
and added fail-closed `load/core <= 0.25` checks. Current saturated host is
correctly rejected; smoke timings carry no performance claim.
Decisions: core score remains authoritative; direct binary64 reconstruction is
required only in the pre-fixed safe region; ordinary zero never masquerades as
scaled underflow; Rmpfr's GMP/MPFR dependency surface is unnecessary for the
finite test oracle; cold timing remains recorded context while the warm paired
interval gates default-path regression.
Remaining: implement one-pass dense/sparse component accumulators + additive R
API; verify native outputs/status across representations/backends; run the exact
locked baseline/candidate gate on a clean, quiescent Linux host.
Risks/blockers: no component API exists yet; final warm equivalence evidence
cannot be claimed from the currently saturated host; the native extreme-value
certificate and biological value remain unimplemented/unvalidated.

2026-07-18 | frontier execution 4 | commit `HEAD`
Scope: implement A2's additive component API and one-pass native diagnostics
under the locked numerical/result contract.
Changed: exported `genefunnel_components()` returns the authoritative score plus
aligned sum, penalty, Bulla balance, effective size/fraction, semantic/
availability/conditioning status, and canonical binary-scaled sidecars; dense
and sparse diagnostic kernels read each input chunk once while the score-only
path allocates no diagnostic matrices; native double-double significands with a
separate exponent cover the complete finite binary64 input range without a new
dependency; help, vignette, README, NEWS, and installed specifications now own
the public and interpretive contract.
Verified: independent ordinary/scaled oracles; committed interpretation and
safe-identity cases; sum/penalty overflow, unsafe cancellation, balance
underflow, subnormals, signed zero, missingness, permutations, dense/sparse,
serial/SOCK, malformed native calls, and randomized safe cells; full local and
clean-tarball suites; `R CMD check --no-manual` `Status: OK`; standalone docs;
benchmark protocol smoke; locked A2 baseline/candidate smoke = 16/16 exact
score digests/shapes/environments, byte-identical median R allocations, and RSS
within smoke limits; combined ASan/UBSan and forced 64-bit-`long double` suites;
tarball BiocCheck = 0 errors/0 warnings (advisory notes only).
Decisions: the primary score remains independently calculated and authoritative;
ordinary diagnostics never use infinity or false zero; only out-of-range
nonzero values receive scaled pairs; undefined diagnostics use explicit
`unavailable`; retain package version `0.99.0` while reopened submission scope
stabilizes.
Remaining: run A2-1.0.0's exact 30-pair default-path timing/allocation/RSS gate
on a clean, quiescent Linux host, then close A2 only if all four workload bounds
pass. Current host load materially exceeds the protocol's fail-closed limit, so
smoke evidence carries no equivalence claim. Gate mode correctly rejected the
environment before timing at load/core 1.226 versus the fixed 0.25 ceiling.
Risks/blockers: external sustained CPU load blocks the performance gate; remote
macOS/Windows CI awaits the maintainer push; component biological value remains
an unvalidated H1 hypothesis.

2026-07-18 | frontier execution 5 | commit `HEAD`
Scope: freeze compiled-catalogue state/wire/performance contracts before spike
results.
Changed: `GFCAT-1` specifies copied canonical members/features, exact stale-use
rejection, fully checked portable integer adjacency, domain-separated SHA-256
identity, and a bounded custom wire format; protocol `C-1.0.0` fixes four
dense/sparse overlap workloads, 20 balanced pairs, call-five amortization,
one-sided 15% minimum improvement, and object/wire/RSS budgets.
Verified: official R 4.5+ SHA-256 + serialization contracts and current
Bioconductor 3.23/R 4.6 runtime support; TSV dimensions/seeds/order balance;
Markdown links/whitespace.
Decisions: use standard `tools::sha256sum(bytes=)` rather than a third-party,
system-crypto, or hand-rolled surface; hashes label content while exact checks
carry reuse safety; float-free wire state eliminates missing/signed-zero
ambiguity; public APIs wait for spike gates.
Remaining: implement the internal constructor/scorer/codec + tamper, stale,
serialization, and fresh/reused SOCK tests; then run `C-1.0.0` on a clean,
quiescent host. A2's locked performance gate remains blocked at load/core above
0.25.
Risks/blockers: the strict 15% all-workload gate may honestly reject compiled
scoring; custom parsing needs adversarial bounds; biological H1 remains open.

2026-07-18 | frontier execution 5 amendment | commit `HEAD`
Scope: correct compiled-catalogue protocol feasibility before results.
Changed: protocol `C-1.0.1` raises the feature universe from 20,000 to 50,000,
making its fixed 2,000 x 25 low-overlap memberships exactly disjoint-capable;
all endpoints, thresholds, seeds, orders, sample/set dimensions, and resource
denominators remain unchanged.
Verified: existing deterministic low-overlap constructor constraint
`sets * size <= features`; all four TSV rows; balanced order codes; Markdown +
whitespace.
Decisions: version the correction rather than silently rewriting `C-1.0.0`;
commit before implementation/performance evidence, so no observed result can
influence it.
Remaining: internal spike + runner, then `C-1.0.1` smoke/full gates.
Risks/blockers: strict benefit gate and host quiescence remain unchanged.

2026-07-18 | frontier execution 6 | commit `HEAD`
Scope: implement the internal compiled-catalogue safety/portability spike.
Changed: unexported constructor + separate scorer preserve the named-list API;
`GFCAT-1` uses bounded big-endian framing and standard SHA-256; complete use-time
validation rejects stale features and mutated schema/members/mappings/coverage/
adjacency/fingerprints before native scoring; portable integer adjacency and
custom bytes contain no pointer/environment.
Verified: fixed cross-process SHA vectors; Latin-1/UTF-8/bytes round trips;
corrupt/truncated/trailing/impossible-length and valid-digest mapping tamper;
12 randomized dense/sparse/list identities; fresh + reused two-worker SOCK;
full fixed low/high catalogue object = 181.9 B/membership and wire = 51.0
B/membership (limits 256/192);
full suite = 75 tests/3,387 expectations; source build; installed-tarball
`R CMD check --no-manual` = `Status: OK`; standalone documentation; tarball
BiocCheck = 0 errors/0 warnings (advisory notes only).
Decisions: keep all functions/classes internal until `C-1.0.1` passes; validate
derived state + fingerprints on every use; rebuild mappings during decode;
retain the custom float-free format instead of long-term R serialization.
Remaining: implement/run the exact `C-1.0.1` paired resource/timing harness;
complete A2's separately locked performance gate when host load/core <= 0.25.
Risks/blockers: validation/hash cost may defeat the 15% benefit gate; current
host load remains ineligible; remote supported-OS evidence awaits push.

2026-07-18 | frontier execution 6 amendment | commit `HEAD`
Scope: make the compiled-catalogue machine manifest fixture-complete pre-run.
Changed: protocol `C-1.0.2` adds fixed dense zero/missing fractions, sparse
stored-missing fraction, and unmatched/duplicate strides already described in
prose; no endpoint, threshold, seed, dimension, repeat, or order changed.
Verified: required-field/types/ranges, disjoint feasibility, exact 50,000
canonical memberships, and balanced orders across all four rows.
Decisions: keep scientific workload constants in the hashed manifest rather
than latent runner code; version before any runner/result exists.
Remaining: implement and smoke the `C-1.0.2` harness; full run awaits a clean,
quiescent Linux host.
Risks/blockers: unchanged strict benefit gate + global load blocker.

2026-07-18 | frontier execution 7 | commit `HEAD`
Scope: make compiled-catalogue protocol `C-1.0.2` executable and auditable.
Changed: exact Git-snapshot preparer + SHA-256 source/installed provenance;
machine-validated dense/sparse fixtures; adjacent paired timing plus fresh,
isolated timing/allocation/RSS workers; ready-synchronized resource-only `/proc`
monitor; compilation/validation/codec/resource facts; 20-pair intervals; fixed
adversarial + fresh/reused SOCK audit; complete machine/human evidence outputs.
Verified: all four full fixtures construct exactly 50,000 canonical memberships;
preparer validates baseline `a573c12` + candidate source/installed trees; final
96-process separate-snapshot smoke produces 32 merged observations, 16/16 paired
digest/fact/environment identities, 13/13 correctness assertions, four summaries,
and `all_pass=TRUE`, `performance_claim=FALSE`; every scoring monitor records
5-9 samples; protocol SHA/schema, four full membership counts, CLI/parse/local
links/whitespace pass. Audit found/fixed duplicate schema, quote reporting,
unsynchronized sampling, mislabelled whole-worker RSS, and resource work between
paired timings.
Decisions: smoke applies correctness/orchestration only; gate uses fresh R
processes, fixed pair order + load checks, one-sided log-ratio intervals, and
fails unavailable/undersampled metrics; worker start/first/reused time remains
descriptive.
Remaining: prepare the committed candidate snapshot; run gate only at clean HEAD
with load/core <= 0.25; check C's benchmark box only if every fixed threshold
passes. A2's independent gate remains open under the same host constraint.
Risks/blockers: smoke compiled/list ratios are intentionally non-evidentiary and
show validation overhead can dominate tiny fixtures; final smoke load/core
oscillates across the strict `0.25` boundary (`0.214-0.276`).

2026-07-18 | frontier execution 7 gate attempt | candidate `315ab7f`
Scope: attempt exact compiled-catalogue gate `C-1.0.2` without weakening its
environment or decision rules.
Changed: no implementation/protocol; prepared and verified separate Git-archive
installs for list baseline `a573c12` and compiled candidate `315ab7f` with
SHA-256 source + installed-tree markers; retained ignored raw rejection evidence.
Verified: clean exact HEAD; 13/13 fixed correctness assertions; gate admitted at
load/core `0.14375`; six complete full-workload pairs (36 isolated workers, 12
merged observations) have exact list/compiled output digests; all load checks
through repeat 2 `sparse_low` pass. Before repeat 2 `sparse_high`, load/core rises
to `0.2875`; machine decision = `all_pass=FALSE`,
`performance_claim=FALSE`, reason `quiescence_rejection:before-r2-sparse_high`.
Decisions: partial timings are inadmissible and receive no threshold summary;
the fixed 20-repeat gate restarts from zero on an adequately quiescent host.
Remaining: rerun C's full gate, then check its benchmark box only if every
correctness/timing/resource bound passes; A2's independent gate also remains.
Risks/blockers: this host crossed the ceiling after six of 80 required pairs;
completion needs sustained load/core <= `0.25`, not a momentary quiet start.

2026-07-18 | frontier execution 8 | commit `HEAD`
Scope: ratify Workstream B's aggregation theorem and pre-API design boundary.
Changed: exact weighted Jensen/coordinate/opposing-mass identities + equality,
strictness, bounds, normalization, physical-sum relation, eligibility, default
rejection/explicit common-support missingness, interpretation, and targeted
prior-art audit; independent base-R oracle + randomized/edge-case properties;
PLAN's first six B gates, README, and NEWS now point to the frozen design.
Verified: 256 randomized weighted fixtures plus proportional/nonproportional
equality, maximal complementary cancellation, zero-weight exclusion, scaling,
physical sums, permutations, and zero-denominator behavior; full source suite =
79 tests/4,927 expectations; standalone documentation; exact source build;
installed-tarball `R CMD check --no-manual` = `Status: OK`; tarball BiocCheck =
0 errors/0 warnings (six advisory notes).
Decisions: weighted mean is primary; zero-weight units leave eligibility before
value/support checks; `J/F(mean)` is only the aggregate-score discrepancy
fraction and is `NA` at zero; field term = aggregation gap, never beta
diversity, generic heterogeneity, complementarity, synergy, coherence, or
Simpson effect;
no biological/public novelty claim follows from the exact theorem.
Remaining: freeze the additive audit schema + controlled/held-out validation
protocol before implementation results, then prototype and test the group API;
rerun C and A2 gates only on a host sustaining load/core <= `0.25`.
Risks/blockers: normalization is unstable near zero and member-scale/noise/
preprocessing dependent; biological value remains an unvalidated H2 hypothesis.

2026-07-18 | frontier execution 9 | commit `HEAD`
Scope: freeze Workstream B implementation/validation protocol before prototype
or empirical results.
Changed: protocol `B-1.0.0` fixes internal signature, four-table schema,
validation/missingness/reason/order/numerical contracts; exact 62,208-latent/
124,416-measurement resolution-IV synthetic design; curve, stability, null,
incremental, known-mixture, donor, multiplicity, and fallback gates; eight-file
SHA-256 data manifest pins CellBench RNA mixtures, Kang GSE96583, and Reactome
v97 without package dependencies; CI smoke validates the machine registry.
Verified: registry = 60 unique rows; design = 1,944 core x 32 orthogonal runs x
two independent measurements with every degree-1/2/3 contrast balanced; all
eight downloaded artifacts = 81,600,675 exact bytes + matching SHA-256;
existing dense/sparse serial/SOCK + controlled benchmark smoke; standalone
documentation; exact source build; installed-tarball `R CMD check --no-manual`
= `Status: OK`.
Decisions: prototype stays internal until controlled, cross-protocol mixture,
and real-data stability gates pass; external benchmark data remain ignored and
dependency-free; CEL-seq2 trains while SORT-seq is held out; Kang inference
uses eight paired donors, never cells/pathways as replicates; biological claims
need the separately frozen held-out + Holm-adjusted exact sign-test gate.
Remaining: implement `.aggregation_audit()` exactly, add API-level adversarial/
oracle coverage, then implement and run synthetic, CellBench, and Kang runners
without changing `B-1.0.0` after results.
Risks/blockers: severe dropout/protocol effects may deliberately reject H2;
data-derived mixture sets and pathway identifier mapping require audited stable
ties/alignment; no theorem guarantees curve accuracy or biological utility.

2026-07-18 | frontier execution 9 amendment | commit `HEAD`
Scope: close a group-audit reporting guarantee gap found during first local
implementation fixtures, before any empirical/protocol endpoint result.
Changed: protocol `B-1.0.1` adds `unit` to `removed_members`; unmatched rows use
`NA`, while missing rows enumerate affected active units in declared-member/
unit order. All external bytes, factors, seeds, splits, endpoints, thresholds,
and decisions remain identical to `B-1.0.0`.
Verified: registry = 60 unique rows; exact 124,416-row synthetic expansion and
resolution-IV checks unchanged; all eight cached external artifacts still
match 81,600,675 declared bytes + SHA-256; Markdown/local-link whitespace.
Decisions: version the correction openly because initial red/green prototype
tests had run; the correction repairs the already-ratified theorem contract and
cannot benefit any unseen scientific outcome.
Remaining: update the internal prototype/detail tests to `B-1.0.1`, complete
adversarial/oracle verification, then run empirical gates without protocol
changes.
Risks/blockers: detail volume can grow with member x missing-unit occurrences;
summary removal count remains unique members and must not be inferred from rows.

2026-07-18 | frontier execution 10 | commit `HEAD`
Scope: implement and adversarially verify the frozen internal aggregation audit.
Changed: `.aggregation_audit()` emits the `B-1.0.1` four-table schema with exact
group/set/unit order, explicit normalized masses, zero-weight exclusion,
reject/intersection support facts, authoritative native scores, an independent
scaled opposing-mass gap, residual certification, and fail-closed numerical
states; pure helpers keep every function below 50 lines and avoid non-local
mutation. Added schema/reason/missingness/extreme/dense-sparse fixtures plus 128
random group audits against the independent oracle and direct mass, scale,
physical-sum, equality, and permutation checks. PLAN/NEWS/spec record only an
unexported prototype; empirical promotion gates remain open.
Verified: full source suite = 89 tests/7,083 expectations with zero failures;
standalone Rd/link/vignette documentation; exact source build; installed-tarball
`R CMD check --no-manual` = `Status: OK`; tarball BiocCheck = 0 errors/0
warnings/6 standing advisory notes, with no function-length or `<<-` note.
Decisions: within-group max-scaled effective weights reject positive-mass
underflow; inactive values cannot affect validity; unexpected
native errors propagate while arithmetic range failures become explicit
`numerically_unavailable`; independent formula and native identity must certify
before any metric is returned; prototype remains internal and dependency-free.
Remaining: implement/run the frozen 124,416-measurement synthetic experiment,
then the pinned CellBench mixture and Kang/Reactome validations without changing
`B-1.0.1` in response to outcomes.
Risks/blockers: empirical H2 may fail one or every co-primary threshold; this R
prototype materializes one retained group/set block and is evidence machinery,
not a performance or public-API commitment.

2026-07-18 | frontier execution 10 amendment | commit `HEAD`
Scope: remove synthetic execution degrees of freedom before any empirical row.
Changed: protocol `B-1.0.2` keeps the B-1.0.1 API, factors, external bytes,
endpoints, thresholds, and decisions; it fixes a shared latent RNG distinct from
A/B measurement seeds, exact overlap/subunit/dropout/zero-total/covariate
operators, paired factor-balanced folds, categorical main-effect encoding,
deterministic QR alias handling, type-8 quantiles, and a within-fold paired
prediction-error bootstrap. The data manifest version changes; its eight URLs,
byte counts, and artifact SHA-256 values do not.
Verified: registry = 79 unique rows + matching SHA-256; design = 62,208 paired
latent scenarios/124,416 measurements; ten latent fold counts differ by at most
two and every declared main-effect level by at most seven; all eight cached
artifacts still match 81,600,675 bytes + declared SHA-256; benchmark/controlled
CI smoke and standalone documentation pass.
Decisions: rejected the first latent-ID modulo fold rule before use because it
perfectly confounded fold parity with the OA baseline factor; mix one-indexed
core and factorial-run indices instead. Bootstrap resamples paired held-out
prediction errors without refitting; all ten CV fits remain the primary model
procedure. No synthetic measurement, model fit, or downloaded-data read
preceded this amendment (manifest hash verification reads bytes only).
Remaining: implement, adversarially smoke-test, and commit the B-1.0.2 runner;
only then execute all controlled rows and record every passed/failed endpoint.
Risks/blockers: regression aliasing is expected because maximum weight is fixed
by declared design factors; frozen pivoted QR chooses one prediction-equivalent
coefficient representation, so coefficients are descriptive rather than causal.

2026-07-18 | frontier execution 11 | commit `HEAD`
Scope: implement B-1.0.2 synthetic generation, analysis, and evidence retention
before running its controlled observations.
Changed: base-R latent/archetype/log-normal/overlap/outlier generator; exact
multinomial subunit + member-dropout measurement; independent reference target;
two-group internal audit call; deterministic serial/fork execution; atomic
protocol/Git/scenario-bound checkpoints; fixed treatment matrices, ten paired
QR fits, within-fold bootstrap, endpoints/strata/report/artifact hashes. Full
runner requires a clean commit, installs that tree in an isolated temporary
library, resumes only identity-matching chunks, records scientific `FAIL` as a
valid outcome, and aborts mechanistic/provenance failures. Every helper is <=49
lines; CI smoke covers generator, audit, full 124,416-row fabricated model
orchestration, endpoint schema, and deterministic replay.
Verified: pre-implementation source seam failed because the runner was absent;
four paired strata pass serial/fork/checkpoint identity, shared-target, complex
zero-gap, dropout/subunit, metric-bound, and seed assertions; a fabricated full
design passes ten folds/62,208 predictions, explicit factor contrasts, QR alias
handling, bootstrap, ten endpoints, and 35 strata; dirty-tree execution refuses;
isolated self-install exposes the committed internal audit; complete 89-test/
7,083-expectation source suite, benchmark/controlled/aggregation CI smoke, and
standalone documentation pass.
Decisions: one audit call scores both eligible A/B groups per latent scenario;
workers never own RNG streams because every latent/measurement seed resets;
post-result thresholds are immutable and hard-checked against the registry;
generated checkpoints/results remain ignored and never become sole claim proof.
Remaining: commit this runner while no controlled observation exists, execute
all 124,416 measurements from that clean commit, adversarially audit artifacts,
then record every passed/failed gate without changing B-1.0.2.
Risks/blockers: full runtime and checkpoint volume are not yet measured; severe
dropout may exclude pairs or reject curve/stability/incremental gates by design.

2026-07-18 | frontier execution 12 | commit `HEAD`
Scope: execute, adversarially audit, and compactly retain B-1.0.2's frozen
controlled factorial result.
Changed: tracked exact ten-endpoint result + provenance/failure-envelope note;
PLAN checks only the controlled factorial item; README, NEWS, and aggregation
spec distinguish a synthetic PASS from dropout robustness or biology. The
ignored full bundle contains 62,208 paired latent scenarios, 124,416 measured
rows, 2,000 bootstraps, 486 resumable checkpoints, and all model/audit facts.
Verified: clean source commit
`d1cbf1510af3bc27d6baf14e3433866449645a77`; protocol SHA-256
`aec8fd4e3e49b953e5ca75e0c2059c1d68436409def0dfd3d351b4ad8a49356f`;
data SHA-256
`ca6dedb4957fa6ef3c21649d342996d1d482cf6f3203e4f2f9bab0bbf1a823e8`;
124,331 eligible measurements, 85 zero-observed-unit measurements, and 85
model-excluded pairs. Curve median/q90 error = 0.0102821/0.182020; overall +
three archetype A/B Spearman = 0.929844/0.919627/0.930966/0.936522; null
median/q95 = 0.000149946/0.0182094; median fold RMSE reduction/bootstrap 95%
lower = 0.218438/0.210229. Every threshold passed. Binary replay of all
checkpoints reproduced endpoint/fold/prediction/bootstrap/coefficient tables
within `5e-16`; 14/14 artifact hashes matched.
Decisions: controlled co-primary status = PASS. At 0.5 dropout, q90 curve error
= 0.304640 and median observed/latent gap = 0.153459/0.053335; the marginal
dropout gap rise 0.09762 is comparable with complementarity 0-to-1's 0.10728.
Treat this as the PLAN fallback concern: retain a severe-dropout warning and no
robustness/public/biological claim, regardless of the aggregate gate pass.
Remaining: implement and run the pinned CellBench CEL-seq2 training/SORT-seq
held-out mixture gate, then Kang/Reactome technical and donor gates; keep the
prototype internal.
Risks/blockers: external preprocessing/alignment may fail closed; CellBench or
Kang can still reject H2, and no threshold may move in response.

2026-07-18 | frontier execution 13 | commit `HEAD`
Scope: close external-data execution degrees of freedom before any CellBench or
Kang profile, selected set, audit metric, correlation, contrast, or endpoint.
Changed: protocol `B-1.0.3` leaves API/synthetic fields, completed B-1.0.2
result, eight payloads, thresholds, donors, pathways, and decisions unchanged.
It fixes CellBench all-row CPM, pooled pure-profile means, common-gene order,
zero/zero ratio semantics, stable selection, equal-mixture normalization,
every-group curve gate, 192-condition split/cross-platform units, C sorting,
and type-8 quantiles. Kang now fixes the two-matrix/barcode join, count/symbol
collapse, within-cell-type halves, retained-unit weights, GMT filtering/order,
16 stability units, training direction, held-out zeros, and exact sign tests.
Verified: pre-endpoint schema audit only - CellBench count/metadata orders are
exact and each platform has the fixed 4 x 4 mixed composition/amount grid with
at least three libraries; Kang has 35,635 gene rows and 14,619 + 14,446 matrix
columns. Concatenated raw barcodes contain 313 cross-file duplicates; base R
`make.unique(..., sep = "")` exactly reproduces all 29,065 author metadata IDs.
All eight cached payload sizes/hashes remain unchanged. Registry = 116 unique
keys; protocol/data SHA-256 =
`c95a2eb91c9e3f8027b461e85ee532a621550038b27f20167b0838dd95c2f7ad` /
`c8743b696f3e05fa623d7f3adc19e27b379303e7f0544a88e054aedf528d5a28`.
Decisions: malformed/alignment failures abort; missing fixed-grid metrics are
retained scientific FAIL results. Schema/dimension/condition facts were the
only downloaded contents inspected before this amendment; no external score or
endpoint exists. B-1.0.2 controlled evidence is inherited without rerunning.
Remaining: implement/adversarially test the CellBench runner from this committed
protocol, then execute and record it before implementing Kang.
Risks/blockers: strict every-group CellBench error gate can fail on one narrow
amount/composition/set; Kang may lose cell-type units or defined pathways, which
must fail visibly rather than change eligibility after inspection.

2026-07-18 | frontier execution 14 | commit `HEAD`
Scope: implement B-1.0.3 CellBench mechanics and evidence runner before reading
an external pure profile or endpoint.
Changed: fail-closed CSV/count validation; outlier removal before C-locale
splits; all-row CPM; platform pure-profile means; exact common-gene order;
zero-aware stable pair/control set selection; normalized composition reference
audits; unbounded observed process-control gap; strict 384-group curve summaries;
192-condition odd/even correlations per platform and cross-platform complete
correlation. The clean-tree runner verifies four hashes, installs that commit in
an isolated library, retains every raw/group/condition/endpoint failure, and
writes session + artifact hashes. Shared clean-install/protocol/artifact helpers
are deduplicated from the synthetic runner.
Verified: red source seam failed while the CellBench implementation was absent;
fabricated integer-count CSVs exercise positive/zero, zero/positive, and
zero/zero ranks; 12 sets/672 memberships; 768 observations; 384 curve + 384
condition rows; all five expected endpoints; exact serial audit/scorer use; and
missing-group scientific failure. Negative-count and reordered-metadata
fixtures abort at the intended parser boundaries. Synthetic runner help,
shared file-map smoke, and isolated install/audit/scorer resolution pass; every
new function is <=40 lines.
Decisions: selection reads only training pure profiles; held-out profiles cannot
affect members. Every error stratum must pass rather than letting a global
average rescue it. Undefined values remain rows and fail the fixed grid; input
schema/alignment failures abort. Full downloaded-data execution waits for this
implementation commit, so no CellBench outcome has been inspected.
Remaining: run documentation + complete benchmark/aggregation smoke, commit,
then execute the four pinned files and adversarially audit/record every endpoint.
Risks/blockers: exact data may reject one or both co-primary gates; output must
remain a process-control result, never biological validation.

2026-07-18 | frontier execution 15 | commit `HEAD`
Scope: execute, adversarially audit, replay, and compactly retain B-1.0.3's
pinned CellBench result.
Changed: tracked exact five-endpoint result + negative-evidence anatomy; PLAN,
README, NEWS, and aggregation spec now close public promotion while keeping the
combined external-validation box open for Kang characterization. No protocol,
threshold, selection, missingness, or decision rule changed after inspection.
Verified: clean source commit
`9b3cab65f626200d28a9ce494b79700247eaf0be`; four inputs = 4,072,960 exact
bytes; 13,906 common genes; 12 training-derived sets; 4,848 library/set rows;
all 96 references eligible/defined. Both curve endpoints = `NA`/FAIL because
only 191/384 fixed groups were complete. CEL-seq2/SORT-seq/cross-platform
complete-case Spearman = 0.977979/0.944496/0.909954, but only 98/93/88 of 192
conditions were complete, so all stability endpoints fail. In total 299/384
curve groups fail, including all 288 pair groups; 1,434 observations have exact
zero measured score. Independent selection/formula/grid/endpoint reconstruction
passes within TSV precision; 16/16 artifact hashes match; a second isolated run
reproduces all 12 scientific tables byte-for-byte.
Decisions: CellBench co-primary decision = FAIL. High finite correlations are
descriptive and cannot override completeness. Pair-set finite errors are also
far beyond threshold, so this is not merely a missingness technicality. The
B-1.0.3 export gate is closed; retain theorem + aggregate-then-score and severe-
dropout/small-denominator warnings. Kang remains pre-specified characterization
only and cannot rescue promotion.
Remaining: implement/adversarially test the committed Kang/Reactome runner,
then execute and record technical stability, held-out directions, and exact
donor sign tests without changing B-1.0.3.
Risks/blockers: Kang may fail alignment, unit eligibility, pathway retention,
technical stability, held-out replication, or biological-claim gates; none
alters the already-negative CellBench decision.

2026-07-18 | frontier execution 16 | commit `HEAD`
Scope: implement B-1.0.3 Kang/Reactome mechanics and evidence runner before
calculating a downloaded-data pathway score or endpoint.
Changed: exact `coordinate real general` Matrix Market parsing rejects declared
dimension/entry/trailing-record mismatches and coalesces duplicate coordinates;
the frozen `make.unique(sep = "")` join must equal metadata order. Empty symbols
drop, duplicate symbols sum in first-occurrence order, and the fixed 96-unit
grid assigns C-locale odd/even halves without non-local mutation. Eligible raw
UMI profiles feed full/odd/even audits; summaries retain the fixed 16 technical
correlations, 16 donor contrasts, training direction, four held-out donors,
exact two-sided donor sign tests, Holm adjustment, and six fail-closed endpoints.
The clean-tree runner verifies four hashes, installs that commit in isolation,
and records preprocessing, assignments, catalogue membership, every audit,
environment, report, and artifact hashes.
Verified: the red source seam failed while the Kang implementation was absent.
A 1,280-cell/12-pathway sparse fixture exercises 96 fixed units, 32 eligible
units, 576 audit rows, exact split stability, positive training/held-out
directions, and both adjusted sign tests; all six fabricated endpoints pass.
Deleting one odd donor-condition produces a retained incomplete technical
failure; zero eligible units fail all six endpoints without aborting, and an
exact-zero training median cannot set a direction. Duplicate matrix coordinates
sum; negative counts, integer/trailing Matrix Market records, reordered barcode
expectations, and count/metadata misalignment abort. Complete benchmark,
controlled, synthetic, CellBench, and Kang CI smoke passes; every new function
is <=40 lines.
Decisions: cells remain aggregation inputs only; donor contrasts remain the
inferential units. Every eight-donor contrast must be defined for the biological
claim, and held-out completeness cannot be rescued by pathway/cell counts.
Full downloaded-data execution waits for this implementation commit.
Remaining: commit, execute the four pinned payloads from the clean snapshot,
adversarially reconstruct endpoints, replay, and track the compact result.
Risks/blockers: fixed data may fail parsing, eligibility, stability, held-out,
or sign-test gates; no observed result may alter B-1.0.3.

2026-07-18 | frontier execution 17 | commit `HEAD`
Scope: repair the Kang metadata boundary exposed by the first pinned run before
any pathway score or endpoint was calculated.
Changed: exact required columns + barcode order remain structural gates;
missing author classification labels now follow the frozen registered-cell/
singlet filter and are explicitly excluded. The retained-cell predicate is
total (no `NA` output), and preprocessing evidence counts missing cell-type
labels plus all registered singlet candidates before unit eligibility.
Verified: all 29,065 author metadata row names exactly equal the frozen joined
barcode order; required fields exist; only `cell` has missing values (9 rows),
which are outside the filter. The corrected boundary yields 23,981 registered
singlet candidates. Fabricated missing-label, malformed matrix/barcode,
misalignment, zero-eligibility, endpoint, and evidence-writer controls pass.
Decisions: unknown classification is an excluded biological label, not schema
corruption. This implements the committed `registered ... only` filter without
changing a threshold, donor, pathway, or decision rule; no endpoint was visible.
Remaining: commit the repair, rerun pinned inputs from that clean snapshot, then
audit/replay and track the result.
Risks/blockers: later parsing or scientific gates can still fail independently.

2026-07-18 | frontier execution 18 | commit `HEAD`
Scope: execute, independently audit, replay, and compactly retain the final
B-1.0.3 Kang/Reactome characterization.
Changed: tracked six-endpoint result + support/failure anatomy; PLAN closes the
ordered external-validation item as executed while selecting the theorem-only
fallback. README, NEWS, aggregation protocol/spec, and result interpretation
reject perturbation sensitivity, a biological claim, and public audit promotion.
Verified: clean source `62d42dc3652d810b3c7ea58c70d29986b42504ce`;
four inputs = 77,527,715 exact bytes. The 313 duplicated raw barcodes join all
29,065 metadata rows exactly; nine missing cell labels filter out. Of 23,981
registered singlet candidates, 23,649 cells in 84/96 eligible units remain.
First-occurrence collapse maps 35,635 rows to 32,938 symbols; 1,728/2,868 human
Reactome pathways (59,969 memberships) retain size 8-128. All 82,944 fixed
pathway/view rows are audit-eligible; 80,633 normalized gaps are defined and
2,311 exact-zero aggregates retain `NA`. All 16 split correlations are defined
over 1,621-1,690 common pathways; median/q10 = 0.916608/0.839918, so both
technical endpoints pass. Interferon alpha/beta fixes positive direction and
matches 3/4 held-out donors, but 6/8 positive signs give exact/Holm p
0.2890625/0.578125. Interferon gamma fixes negative direction, matches only 2/4,
and has 4/8 positive signs, p=1. Endpoint vector = pass/pass/pass/fail/fail/fail.
Independent artifact, unit, pathway, score-identity, contrast, binomial/Holm,
and endpoint reconstruction passes: all 19 artifact hashes match, maximum
relative identity residual = 9.30e-15, and rounded-table rank-tie delta <=4.34e-5
cannot approach a gate. A second isolated run reproduces all 12 scientific
tables byte-exact; ordered hash-stream SHA-256 =
`040d2beca3c0bf41314c6c152fdffd666601bd037e023cab9dfd09aaa837168b`.
Decisions: technical repeatability alone cannot establish donor-replicated
sensitivity. Held-out gate = FAIL; biological gate = FAIL. CellBench had already
closed export; Kang independently withholds perturbation/biological claims.
Workstream B ends with the theorem, internal prototype as maintenance permits,
and aggregate-then-score + severe-dropout/small-denominator warnings.
Remaining: continue PLAN at the next unchecked executable item; no Workstream B
empirical gate remains.
Risks/blockers: none for B execution; public promotion requires a newly frozen
future protocol, not reinterpretation of B-1.0.3.

2026-07-18 | frontier execution 19 | candidate `2694d56`
Scope: resume A2's exact default-path performance gate after Workstream B.
Changed: no package/protocol/threshold. Prepared new fingerprinted Git-archive
installs for locked baseline `9b60a3e` and clean current candidate `2694d56`.
Verified: four-repeat smoke = 16/16 paired score digests/shapes/environments
exact; baseline/candidate median manager allocations are byte-identical in all
four workloads; RSS stays inside smoke bounds; machine smoke decision passes.
The first full-gate admission sampled load/core 0.25625 and a fresh restart
sampled 0.30125, both above the frozen 0.25 ceiling. Both rejected before the
first timing; no partial performance estimate exists.
Decisions: retain A2 unchecked and make no equivalence claim from smoke. Keep
the prepared snapshots as rejection provenance; refresh the candidate marker at
the then-current clean HEAD before a sustained-quiescence restart. Crossing the
admission boundary by a small amount is still a rejection.
Remaining: rerun all 30 pairs/workload from zero when every load check can pass.
Proceed meanwhile to Workstream E's non-performance definition/protocol work.
Risks/blockers: external host load remains the sole A2 execution blocker; the
gate can also reject later if load rises after admission.

2026-07-18 | frontier execution 20 | commit `HEAD`
Scope: fix Workstream E's observed-member deletion quantity and interpretation
boundary before any sensitivity API or empirical diagnostic exists.
Changed: new normative pre-API sensitivity specification; direct brute-force R
oracle; property/canonical-case tests; PLAN, scientific-spec, README, and NEWS
cross-links. For each cell with at least three observed members, delta is full
score minus the score after deleting that member. Compact mathematical facts
are largest absolute delta/member, its signed and observed-sum-normalized delta,
median absolute delta, and effective size; exact ties use stable declared-member
order after deduplication, matching, and per-sample missingness omission.
Verified: ordinary-range direct equation checks lock positive/negative deltas,
all-zero/equal/insufficient-support states, midpoint median, duplicate/unmatched/
`NA`/`NaN` filtering, 100 randomized bounds/homogeneity/common-shift cases,
permutation equivariance, and explicit non-additivity. Clean source tarball
build passes; rebuilt-tarball `R CMD check --no-manual` = `Status: OK`.
Decisions: call the quantity sensitivity, not contribution, importance,
jackknife inference, or causal effect. Normalize the selected signed delta by
the full observed sum only when positive; retain sign after absolute-value
selection. No hidden member-by-set-by-sample array. The future prototype must
freeze extreme-number representation and held-out reliability rules first.
Remaining: commit; freeze E's numerical/API and held-out validation protocol,
then implement compact summaries against brute force. Retry A2 separately only
after the host satisfies its frozen quiescence admission rule.
Risks/blockers: algebra alone establishes no predictive reliability. Public
promotion remains conditional on representation invariance and pre-specified
held-out incremental-effect/technical-repeat gates.

2026-07-18 | frontier execution 21 | commit `HEAD`
Scope: freeze E-1.0.0's internal schema, numerical contract, profiling rule,
controlled reliability design, held-out models, and rejection thresholds before
implementing or calculating a package sensitivity diagnostic.
Changed: 90-row byte-pinned machine registry + normative protocol + fail-closed
parser/design validator; PLAN closes the pre-specification item; sensitivity
spec adds current deletion/leading-edge/dispersion/missing-gene/refinement prior
art and a narrow novelty boundary. The internal seven-field matrix schema uses
exact common-power-of-two arbitrary-integer score numerators: all deletion
numerators share one denominator, so sign/order/ties are exact before ordinary
double conversion. A fixed profile gates optimization research. Controlled
Poisson/NB/dropout repeats cross 5,760 latent profiles with 345,600 independently
masked global-absence/sample-missing rows and five abundance/detection mask
mechanisms.
Verified: registry MD5 `09d9429ae18f64ba00b3bd84955ad71d`;
5,760 scenarios; 345,600 feature rows; 5,760 repeat rows; 576 scenarios in each
of ten folds; every mask leaves >=3 members; 172,800 mask seeds are unique and
in R's integer range. An intentional byte mutation aborts before design use.
Complete benchmark/controlled/aggregation/CellBench/Kang/sensitivity-protocol
smoke passes. Clean source build passes; rebuilt-tarball
`R CMD check --no-manual` = `Status: OK`; every new function is <=40 lines.
Decisions: the baseline includes normalized score (existing balance), making
the increment stricter than coverage/set-size/sum/support alone. Feature-loss
and controlled-repeat median held-out RMSE reductions must each reach 10% with
a one-sided bootstrap lower bound >=5%; neither rescues the other. Synthetic
dropout is a stress operator, not a meaningful-zero or real-repeat claim. A
complete E-1.0.0 pass permits only a newly frozen external stage, never export.
Remaining: commit; implement the exact brute compact prototype and exhaustive
schema/oracle/storage/backend tests, then run the fixed profile before any
sorted-prefix/native optimization. A2 remains independently load-blocked.
Risks/blockers: arbitrary-integer brute recomputation may be slow by design;
the profile decides whether optimization work is justified. Either predictive
gate may fail, selecting the no-public-API fallback without changing algebra.

2026-07-18 | frontier execution 22 | commit `HEAD`
Scope: implement E-1.0.0's exact brute compact sensitivity prototype without
creating a public interface or changing the authoritative scoring path.
Changed: dependency-free base-2^26 arbitrary integers; exact binary64 dyadic
decomposition; independently recomputed full/deleted score numerators; exact
delta sign, magnitude order, ties, and midpoint median; error-certified ordinary
double conversion; explicit semantic/raw/normalized statuses; fixed aligned
matrix schema through bounded dense/sparse BiocParallel iteration. PLAN closes
only the compact-prototype item. Shared chunk ordering now names sensitivity
failures accurately.
Verified: focused sensitivity suite passes. Independent small-integer arithmetic
and 200 randomized integer deletion identities agree exactly; >1,000 finite
binary64 values round-trip across normal/subnormal exponents; a hex-literal
fixture proves canonical exact ordering where subtract-rounded scores select the
later member; wide homogeneous exponents, smallest subnormal, maximum finite,
signed zero, randomized direct recomputation, stable duplicate/unmatched/
missing support, integer/dense/sparse Matrix, and serial/two-worker SOCK paths
pass. Every new package function is <=40 source lines. Complete benchmark/
aggregation/sensitivity-protocol smoke passes. Clean source build passes;
rebuilt-tarball `R CMD check --no-manual` = `Status: OK`.
Decisions: the brute implementation remains the correctness oracle and stores no
member cube. A nonzero raw summary outside the frozen relative-error contract
makes all raw/member fields unavailable; normalization remains separately
statused when raw fields are ordinary. Exact zero is positive zero, all-zero
ties choose the first canonical member, and sets with two global matches stay
aligned as sample-level `not_applicable`. No reliability or biological claim
follows from implementation correctness.
Remaining: commit; run the fixed serial E profile before any optimization. If
the frozen >60 s / >50% exact-arithmetic rule admits optimization, validate any
candidate byte-for-byte/status-for-status against this brute oracle first; then
execute the controlled feature-loss/repeat gates. Retry A2 only under its
independent quiescence gate.
Risks/blockers: brute deletion is intentionally quadratic in observed members;
the pending frozen profile, not intuition, determines whether that warrants a
second implementation. Predictive gates may still select the no-public-API
fallback.

2026-07-18 | frontier execution 23 | commit `HEAD`
Scope: make E's fixed profile executable without silently choosing the fixture
or measurement semantics after observing timing.
Changed: byte-pinned supplement `E-P-1.0.0` records that parent E-1.0.0 fixed
dimensions/seeds but omitted constructors/pass boundaries. It now fixes
column-major unit-rate exponential values, independent within-set sampling,
stable identifiers, isolated clean-Git installation, three explicit-GC/
`gcFirst=FALSE` timed calls, separate `Rprof`/`Rprofmem` passes, exact-stack
attribution, five-output identity, and the unchanged strict 60 s/0.50 OR rule.
New runner records raw profiles, allocations, environment, source/install
fingerprints, artifacts, report, and an optimization-only decision.
Verified: profile registry = 26 unique rows, MD5
`ff86e032b7a766be20c5d61d940991ab`, parent MD5 remains
`09d9429ae18f64ba00b3bd84955ad71d`. Deterministic smoke fixture MD5 =
`6670e2905a56b4c7e612d378ae6d4bcb`; isolated downscaled execution gives one
output MD5 `8d3a853be4e5ac0a20fb86a659c05530` across timed/CPU/allocation
passes. The first rehearsal exposed `system.time()`'s implicit GC contaminating
six of eight samples; freezing `gcFirst=FALSE` before the fixed call gives two
exact samples of three. Dirty-tree gate rehearsal rejects before output creation.
Documentation and complete benchmark/controlled/aggregation/sensitivity
protocol smokes pass; every new function is <=35 source lines.
Decisions: record the parent omission rather than pretending its seeds uniquely
defined a fixture. The supplement is prospective for profiling but explicitly
post-implementation; it changes no reliability design or public-promotion rule.
CPU attribution counts any sampled stack with a `.gf_` frame. Explicit GC occurs
before, outside, each measured pass.
Remaining: commit the supplement/runner, then execute the full profile exactly
once from that clean SHA and retain its decision. Optimization may begin only if
the fixed median or exact share crosses its parent threshold.
Risks/blockers: allocation tracing may add substantial runtime/log volume but is
not a timing observation. The descriptive profile is machine-specific and can
authorize optimization research only, never a reliability claim.

2026-07-18 | frontier execution 24 | commit `HEAD`
Scope: execute E-P-1.0.0's fixed exact-brute profile once from its clean runner
commit and apply the frozen optimization-only decision.
Changed: tracked compact Markdown/machine result plus byte-pinned result parser;
PLAN/protocol/spec/README/NEWS now distinguish completed profiling from pending
correctness-gated optimization and empirical reliability work.
Verified: candidate `37120c125efbdaf3e2c3b30e8049465c4874460b` installed
from Git archive SHA-256
`714bac6613456dcc29940be2940821949248d6ee7bb54f10754510fdb3e41082`;
installed manifest MD5 `d85cf5751d014f656f2b980ccffbb371`; fixed fixture MD5
`32fa8d4843d024765a4adfc676793dbe`. Three elapsed calls = 239.364,
197.017, 213.585 s; median = 213.585 s. The separate 10 ms CPU pass has
23,446 samples, 23,432 exact stacks, share 0.999402883221019. Allocation
pass = 607,529 numeric events / 378,561,864 cumulative manager R bytes. All
five results share MD5 `3d9635e779a9ed1eee453a2a04596369`; every core
artifact hash re-verifies. Raw Rprof/Rprofmem SHA-256 values are retained in the
tracked result. Result registry MD5 = `bfe706bc325fffa3d33b985df962ebc7`.
Decisions: elapsed >60 s = TRUE; exact share >0.50 = TRUE; optimization eligible
= TRUE; performance claim = FALSE. Profile/allocation overhead is excluded from
primary timing. Machine-specific cost cannot imply reliability or API value.
The exact brute implementation remains the oracle and fallback.
Remaining: commit the result; implement a simpler exact optimized algorithm,
then require result/status identity across all fixed/random/extreme/storage/
backend tests plus the frozen workload before adoption. Only afterward execute
the E-1.0.0 controlled feature-loss and measurement-repeat gates. A2 remains
independently load-blocked.
Risks/blockers: optimization may fail equivalence or prove too complex to
maintain, in which case retain the brute internal oracle. Even successful
optimization cannot rescue failed predictive gates or justify export.

2026-07-18 | frontier execution 25 | commit `HEAD`
Scope: implement the first exact sorted-prefix candidate while retaining the
profiled brute path as executable oracle and fallback.
Changed: exact member sorting + arbitrary-integer prefix sums; per-deletion
binary search for the strict below-mean boundary; exact removal of the deleted
member from boundary count/sum; stable bottom-up exact merge order for median
deltas; allocation-free fast path for already-trimmed limbs. The internal cell
path uses the candidate; `.gf_exact_deltas()` remains unchanged as brute oracle.
Verified: red test first failed on absent sorted implementation. Candidate exact
objects are identical to brute for canonical/zero/subnormal/maximum/wide-range
fixtures, 250 seeded cells through size 40, and ten size-128 cells; stable equal-
magnitude ordering is explicit. Existing randomized/direct/schema/extreme/
dense/sparse/integer/serial/SOCK tests pass with a mocked brute failure, proving
the API path bypasses brute. Frozen fixture output MD5 remains
`3d9635e779a9ed1eee453a2a04596369`. Exploratory dirty-tree calls moved from
24.749 s before exact merge/trim work to 12.582 s after it, versus the committed
brute median 213.585 s; these are optimization context, not a speed claim.
Clean source build passes; rebuilt-tarball `R CMD check --no-manual` =
`Status: OK`; all new functions remain <=35 lines.
Decisions: sort original finite binary64 values only to establish their exact
order; every threshold, sum, numerator, sign, delta order, and output conversion
remains arbitrary-integer/rational. Stable sorting preserves canonical ties.
The candidate is not adopted in PLAN until a clean archived commit reproduces
the frozen full-workload digest.
Remaining: commit; install that clean SHA in isolation and record fixed-workload
identity. If it passes, close PLAN's optimized-algorithm item and proceed to the
controlled reliability runner.
Risks/blockers: a clean-candidate mismatch rejects adoption despite randomized
success. Further optimization is optional; correctness and maintainability take
precedence over exploratory elapsed time.

2026-07-18 | frontier execution 26 | commit `HEAD`
Scope: close the optimized-algorithm adoption gate from a clean isolated
candidate and preserve a compact machine-checkable identity result.
Changed: tracked E-O-1.0.0 Markdown/TSV adoption record + MD5-pinned validator;
PLAN closes the optimized-algorithm item; protocol/spec/README/NEWS call the
sorted-prefix path adopted while retaining the brute oracle.
Verified: clean candidate `e5c013fcb0f10e2eb4f0431c1cbd09cd40f73954`
archive SHA-256
`11b7892eb9d60bea77ee63301673da0bae5e0a2744718952a88301f6ede79224`,
archive MD5 `ac7142c66b47d2e5dd5976a1656999bb`, installed manifest MD5
`7f30f8246b98f6a1262a55b2f5217c09`. Fixed fixture MD5 remains
`32fa8d4843d024765a4adfc676793dbe`; candidate output MD5
`3d9635e779a9ed1eee453a2a04596369` equals the committed brute output.
The isolated single call took 15.128 s as non-claim context. Result registry MD5
= `baece34f7435fff0086da286021745cf`; complete protocol smoke validates it.
Decisions: adopt exact sorted-prefix internally. Keep brute executable and its
tests mandatory. `performance_claim = FALSE`; a single correctness call neither
replaces the frozen profile nor creates a performance threshold. No empirical
reliability/public-API gate has run.
Remaining: commit; implement E-1.0.0's controlled scenario generation, feature-
loss/repeat observations, fixed held-out models/bootstrap, and fail-closed report
runner, then execute all 345,600 feature-loss + 5,760 repeat rows.
Risks/blockers: the controlled design may be compute-heavy or reject either
co-primary prediction endpoint. A failed endpoint selects the internal/no-public-
API fallback without affecting exact diagnostic correctness.

2026-07-18 | frontier execution 27 | commit `HEAD`
Scope: close E-1.0.0's remaining controlled-execution degrees prospectively,
before constructing a fixed-grid latent scenario or observation.
Changed: byte-pinned supplement `E-C-1.0.0` fixes exact R RNG calls, canonical
identifiers and row order, installed scoring/sensitivity paths, absence-specific
coverage facts, normalized predictors, dummy/scaling/QR conventions, bootstrap
draw order, clean isolated installation, resumable checkpoint identity, failure
semantics, and evidence artifacts. Parent registry bytes remain unchanged.
Verified: 45 unique supplement rows; MD5
`1c119004ee749b4242f1b73ca6fd8c4c`; parent MD5 remains
`09d9429ae18f64ba00b3bd84955ad71d`; fail-closed parser rejects byte drift and
validates the parent link. Documentation and CI smoke include both identities.
Decisions: record rather than conceal that “Bernoulli,” “global schema,” and a
fixed bootstrap seed did not uniquely specify an R run. The supplement changes
no controlled dimension, equation, target, predictor concept, endpoint,
threshold, claim boundary, or promotion rule. Scientific failure remains a
valid completed result; provenance/schema/arithmetic failures abort.
Remaining: commit; implement deterministic scenario/mask/observation mechanics
against this supplement, then implement held-out models/bootstrap/report and
execute the clean full grid.
Risks/blockers: full execution remains compute-heavy and either co-primary gate
may fail. No controlled value has been generated at this checkpoint.

2026-07-18 | frontier execution 28 | commit `HEAD`
Scope: implement E-1.0.0/E-C-1.0.0's deterministic controlled latent,
measurement, mask, and observation mechanics without fitting or inspecting an
endpoint model.
Changed: scenario-local exact RNG reset; four registered profiles with
log-normal perturbation/capture; Poisson/NB A/B counts and post-count Bernoulli
dropout; five independently seeded weighted mask mechanisms; canonical selected
member records; one authoritative installed score batch per scenario; exact
full/partial A diagnostics; paired absence encodings; fixed 43-field feature and
36-field technical rows; fail-closed schema/seed/order/target/status/ratio/
encoding validation.
Verified: four smoke strata span sizes 8/32/128, all four archetypes, both count
families, dropout/no-dropout, and folds 1/2/9/10. Two runs remain object-identical
after unrelated RNG use; combined digest = `0d925728a6e151e9e2f67fe283a6a92a`;
240 feature + 4 technical rows pass every invariant. Explicit global row absence
and sample `NA` omission produce identical installed-package scores and exact
diagnostics. Direct equation scores agree within the fixed ordinary tolerance;
an altered target fails validation. A 31.069 s mechanics-only sweep constructs
all 5,760 latent/A/B states and 172,800 unique seeds/masks; every selected index
set has the registered size/order/domain, covering 2,822,400 selections with no
invalid or zero-total full measurement.
Decisions: compute every unique partial input once, then reuse its values,
score, and diagnostic across encodings while deriving their distinct factual
coverage fields. Declared size remains the original full set size. Repeat B is
outcome-only; all predictors and diagnostics use A. No controlled gate or
reliability claim has been calculated.
Remaining: commit; exhaustively exercise latent/mask mechanics, implement the
ten-fold models, fixed-prediction cluster bootstrap, resumable clean runner, and
report; then execute the complete grid once from its clean candidate.
Risks/blockers: model implementation must preserve the supplement's global
schema/scaling/alias conventions. Full exact sensitivity remains the dominant
compute cost and either co-primary gate may validly fail.
