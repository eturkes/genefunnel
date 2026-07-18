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
