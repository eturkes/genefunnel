# GeneFunnel fidelity and release plan

**Repository:** `eturkes/genefunnel`\
**Observed baseline:** `main` at `e413316` on 2026-07-14\
**Primary goal:** turn the current algorithmic prototype into a scientifically faithful, sparse-safe, tested, documented, reproducible R/Bioconductor package.\
**Source of scientific intent:** Emir Turkes, *Development of Gene Set Enrichment and Imputation Methods for Transcriptomics and Proteomics: Application in the Study of Neurofibrillary Tangle-bearing Neurons in Alzheimer’s Disease* (UCL PhD thesis, 2025), especially Sections 4.3-4.12 and 5.1.

> Codex will not receive the thesis. This file is the normative substitute for every implementation-relevant thesis requirement. Page/section references are provenance only; no task may block on obtaining the thesis.

## 1. Operating rules for Codex

1. Start every session from the repository root by reading `AGENTS.md`, this file, `.agent/roadmap.md` if present, `git status`, and the current diff/log.
2. Reconcile this plan with the live `HEAD` before editing. The baseline inventory below describes `e413316`; skip work already completed, but verify it with tests rather than assuming it is correct.
3. Treat the specification in Section 2 as authoritative over current behavior, old documentation, benchmark-viewer behavior, or convenient implementation shortcuts.
4. Preserve a plain numeric matrix as the primary return value. Any proposal to change the equation, missing-value semantics, coverage semantics, non-negativity requirement, or return type requires explicit maintainer approval.
5. End each cohesive session with:
   - targeted tests for the changed behavior;
   - the broadest practical package checks;
   - updates to `.agent/roadmap.md` and the checkboxes in this plan;
   - one scoped commit with no unrelated edits.
6. Keep generated files generated: edit roxygen sources rather than `.Rd`; regenerate `RcppExports.*` using `Rcpp::compileAttributes()`; regenerate `NAMESPACE` using roxygen tooling.
7. Never fabricate reproduction of a thesis result. The thesis data and text are unavailable to Codex. Deterministic unit fixtures and synthetic benchmarks specified here are sufficient for the package work; any later real-data reproduction must be explicitly labelled and supplied with data provenance.
8. Prefer a green repository at every session boundary. Split a session rather than committing known failing tests or half-migrated interfaces.

## 2. Normative scientific specification

### 2.1 Inputs and output

`genefunnel()` accepts:

- `mat`: a feature-by-sample, non-negative numeric matrix. Supported release targets are base R numeric/integer matrices and appropriate `Matrix` dense/sparse classes. Sparse input must remain sparse throughout the R-to-C++ path.
- `gene_sets`: a named list. Each element is a character vector of feature identifiers. Matching is exact and case-sensitive; GeneFunnel performs no identifier conversion, synonym resolution, or annotation lookup.
- `BPPARAM`: a `BiocParallelParam`; retain the intended default behavior with a namespace-qualified default such as `BiocParallel::bpparam()`, subject to current official Bioconductor guidance.

The result is a numeric matrix with retained gene sets as rows and samples as columns, preserving their input order and names.

Gene sets may overlap. Scoring one set must not depend on any other set. A single sample and a single gene set are valid inputs.

### 2.2 Exact score

For input matrix \(X\), declared gene set \(G\), and sample \(j\):

1. Treat a gene set as a mathematical set: remove duplicate member identifiers while preserving first occurrence.
2. Match members against `rownames(mat)`. Let \(M = G \cap \mathrm{rownames}(X)\).
3. Within sample \(j\), remove values that are `NA` or `NaN`. Let the remaining observed values be \(x_1,\ldots,x_n\).
4. If \(n < 2\), return `NA_real_` for this gene-set/sample cell.
5. Otherwise calculate:

\[
\bar{x}=\frac{1}{n}\sum_{i=1}^{n}x_i
\]

\[
\operatorname{GeneFunnel}(x_1,\ldots,x_n)
=\sum_{i=1}^{n}x_i
-\frac{n}{2(n-1)}\sum_{i=1}^{n}\left|x_i-\bar{x}\right|.
\]

The effective \(n\) is sample-specific after omitting missing values. It is not the declared set size and not merely the globally matched set size.

### 2.3 Value semantics

- **Zero is observed information.** Include explicit and sparse implicit zeros in the sum, mean, deviation, and effective set size.
- **`NA`/`NaN` is absent information.** Exclude it from all terms and reduce the effective set size for that sample.
- **`Inf`/`-Inf` is invalid input.** Reject it with an actionable error; do not silently omit it.
- **Negative finite values are invalid input.** The non-negativity guarantee and activity interpretation require non-negative input. Do not shift, rescale, normalize, or log-transform data automatically.
- Processed values are acceptable only when their scale has a meaningful zero and remains non-negative. Documentation must not recommend an arbitrary positive shift merely to pass validation, because an additive shift changes the score.

The distinction between a missing row, a missing cell, and a zero is fundamental:

- a feature absent from matrix rows is globally unmeasured and affects gene-set coverage;
- an `NA`/`NaN` cell is missing only for that sample and changes sample-specific \(n\);
- a zero is measured and remains in the calculation.

### 2.4 Gene-set coverage

The scorer must support incomplete matrix coverage:

- Score every set with at least two unique members matched to matrix row names.
- Ignore unmatched declared members during scoring.
- Omit sets with fewer than two globally matched unique members from the output, preserving the order of the remaining sets. Emit at most one concise aggregate warning about omitted sets.
- Do not hard-code a coverage threshold such as 100% or 50% into `genefunnel()`.
- Provide an exported helper, provisionally `gene_set_coverage(gene_sets, features)`, that reports at least:
  - gene-set name;
  - unique declared size;
  - matched size;
  - unmatched size;
  - coverage fraction;
  - duplicate-member count;
  - whether at least two members are scoreable.
- Documentation must show user-controlled filtering before scoring:
  - strict transcriptomic policy: coverage `== 1`;
  - low-coverage proteomic example: coverage `>= 0.5` and matched size `>= 2`.

This preserves the thesis intent: coverage policy belongs to the experiment and caller, while the algorithm can score any available subset of size at least two.

### 2.5 Identifier and structure rules

- Matrix row names are required, non-missing, non-empty, and unique. Duplicate row names are ambiguous and must be rejected.
- Gene-set names are required, non-missing, non-empty, and unique.
- Every gene-set element must be a character vector with non-missing, non-empty identifiers.
- Duplicate members inside a set are deduplicated rather than double-counted; report them through coverage diagnostics.
- Column names, when present, are preserved exactly. If absent, the output may have `NULL` column names; do not invent biological sample names.
- Zero-row or zero-column matrices are invalid. Non-numeric data frames and arbitrary objects must not be silently coerced through a potentially dense `as.matrix()` fallback.

### 2.6 Required mathematical and behavioral properties

For valid finite non-negative observed values:

- score is non-negative, modulo floating-point tolerance;
- score is no greater than the observed sum;
- equal values have zero deviation and score equal to their sum;
- all-zero values score zero;
- a set with exactly one non-zero value and all other values zero scores zero;
- reordering members does not change the score;
- multiplying every observed value by \(c\ge0\) multiplies the score by \(c\);
- adding a constant \(c\ge0\) to every observed member increases the score by \(nc\), because the deviation term is translation-invariant;
- adding, removing, modifying, or reordering other samples does not change a given sample’s score;
- adding, removing, modifying, or reordering other gene sets does not change a given set’s score;
- adding irrelevant matrix rows does not change existing scores;
- dense and sparse representations of the same data produce the same result;
- serial and parallel backends produce the same result and ordering.

Numerical cleanup may clamp a theoretically zero value to exact zero only within a scale-aware tolerance. A materially negative result should trigger an internal error rather than be hidden by unconditional clamping.

### 2.7 Canonical examples

These examples must appear in tests and, in compact form, in user documentation:

| Observed member values | Expected result | Reason |
|---|---:|---|
| `c(4, 4, 4)` | `12` | no deviation; score equals sum |
| `c(0, 0, 0)` | `0` | measured inactivity |
| `c(4, 0)` | `0` | maximally deviant two-member set |
| `c(4, 0, 0)` | `0` | one non-zero member |
| `c(1, 2, NA)` | `2` | omit `NA`; score `c(1, 2)` with effective `n = 2` |
| `c(4, NA)` | `NA` | fewer than two observed members |
| sparse implicit `c(4, 0)` | `0` | implicit zero is included, not treated as missing |

For a declared set `c("A", "B", "C")` and matrix rows containing only `A` and `B`, score `A` and `B`, report coverage `2/3`, and leave the caller to decide whether that coverage is acceptable.

### 2.8 Scope and scientific limitations to document honestly

The package computes scores; it does not perform normalization, imputation, differential testing, identifier mapping, gene-set retrieval, or biological interpretation.

Documentation must state:

- Within-set absolute variability is treated as evidence against coherent activity; this is a scientific assumption, not a universal biological truth.
- The sum term retains expression magnitude and therefore also retains effects of input units, library size, preprocessing, gene-set size, and coverage.
- Sample independence does not by itself make scores automatically comparable across unrelated datasets. Cross-dataset comparison requires compatible units, preprocessing, feature identifiers, gene-set definitions, and coverage.
- Gene-set database revisions can change scores; users should record the exact set collection/version.
- The distributional suitability of scores for every downstream statistical method is not guaranteed. Users must validate downstream assumptions.
- Raw non-negative measurements are preferred where scientifically appropriate because earlier filtering/transformation can discard zeros or alter magnitude, but the package does not enforce a pipeline position.

## 3. Baseline gap inventory at `e413316`

Codex must verify each item against live `HEAD` before acting.

| Area | Baseline state | Required end state |
|---|---|---|
| Default backend | `BPPARAM = bpparam()` is unqualified | Namespace-qualified and tested default |
| Sparse handling | all inputs pass through `as.matrix()` before conversion back to sparse | sparse inputs never become full dense matrices |
| Dense handling | dense inputs are converted to sparse | dense inputs use an appropriate dense path |
| Missing values | `NA`/`NaN` contaminates sum/mean/deviation | omit missing values per set/sample; recompute effective `n` |
| Partial coverage | R wrapper requires every declared member to be present | score the matched subset when at least two members exist |
| Matching cost | row-name hash and set indices rebuilt for every sample | validate and precompute membership indices once per call |
| Task granularity | one R/BiocParallel task per column | bounded column chunks or another profiled strategy |
| Duplicate row names | silently overwrite in C++ map | reject as ambiguous |
| Duplicate set members | double-counted | deduplicate with diagnostics |
| Non-finite values | positive infinity/NaN can reach arithmetic | reject infinities; treat `NA`/`NaN` as missing |
| Numerical cleanup | fixed absolute epsilon | scale-aware tolerance with detection of material negatives |
| OpenMP | compiler/linker flags present but no OpenMP region | remove flags unless a measured, portable OpenMP path is actually used |
| Metadata | placeholder title, author, description, and license text | complete installable package metadata |
| Namespace | imports `Rcpp::sourceCpp` despite compiled registration | generated minimal imports and registered native routines |
| Source tree | `.o` and `.so` build products tracked | source-only repository; generated build products ignored |
| Quality | no test suite, vignette, CI, release, or benchmark harness | comprehensive correctness tests, docs, CI, reproducible benchmarks, release candidate |
| Claims | README claims efficient sparse operation despite densification | claims limited to demonstrated behavior |

## 4. Target public API

### 4.1 `genefunnel()`

Provisional signature, kept deliberately small:

```r
genefunnel(
  mat,
  gene_sets,
  BPPARAM = BiocParallel::bpparam()
)
```

Do not add scoring knobs. In particular, do not add a hidden normalization option, deviation weight, minimum coverage threshold, missing-value replacement, or set-size filter. Experimental policy belongs in caller-side preprocessing, aided by `gene_set_coverage()`.

Expected behavior:

- validate inputs before launching workers;
- prepare unique member/index vectors once;
- omit globally unscoreable sets with one aggregate warning;
- dispatch to dense or sparse kernel without changing representation unnecessarily;
- preserve set/sample ordering;
- return only the score matrix as the primary object.

### 4.2 `gene_set_coverage()`

Provisional signature:

```r
gene_set_coverage(gene_sets, features)
```

`features` may be a character vector or row names extracted by the caller. The helper must use exactly the same validation, deduplication, and matching rules as `genefunnel()`, so coverage reports cannot drift from scoring behavior.

Do not calculate a single “recommended” threshold. Return facts; show policy examples in documentation.

### 4.3 Internal boundaries

A reasonable target separation is:

- R input/class validation;
- R gene-set canonicalization and coverage/index preparation;
- R backend/chunk orchestration;
- native dense kernel;
- native sparse kernel;
- pure R reference scorer used only for tests and scientific verification.

Exact filenames are flexible, but avoid a monolithic wrapper and avoid exporting native helpers as public API.

## 5. Milestones and Codex sessions

Session boundaries are suggested cohesive commits, not a fixed limit. Split any session whose verification cannot be completed cleanly.

### Session 1 - Establish a trustworthy build and test baseline

- [x] Record live `HEAD`, repository tree, package/tool versions, current `R CMD build` and `R CMD check` failures in `.agent/roadmap.md`.
- [x] Remove tracked `src/*.o`, `src/*.so`, DLLs, and similar build products; add appropriate ignore rules.
- [x] Correct only the minimum package metadata needed for reliable install/check/test execution. Use real author/creator information already present in the repository; never invent an ORCID.
- [x] Add a supported R package test framework after checking current official Bioconductor guidance.
- [x] Add a pure R reference implementation in test helpers for the exact formula on already matched member values.
- [x] Add green smoke tests for current complete, finite, non-negative dense behavior: shape, names, all-equal case, all-zero case, one-nonzero case.
- [x] Ensure package installation from a built source tarball, not only `load_all()`, reaches the smoke tests.

**Acceptance:** clean source tree; package installs; test harness runs; formula smoke tests pass; no intended scoring behavior has changed yet except removal of accidental build artifacts.

### Session 2 - Lock the R-level contract and coverage policy

- [x] Implement strict validation for matrix classes, dimensions, row names, gene-set names, member identifiers, duplicates, negative values, and infinities.
- [x] Remove arbitrary coercion paths that can silently densify or convert non-numeric objects.
- [x] Implement shared gene-set canonicalization: deduplicate members, preserve order, exact-match row names, precompute integer indices once.
- [x] Implement `gene_set_coverage()` with the fields defined in Section 2.4.
- [x] Change the scorer’s global set policy from “all members present” to “at least two unique matched members”.
- [x] Emit one aggregate warning for sets omitted because fewer than two members match; keep partially represented but scoreable sets.
- [x] Fix and validate the namespace-qualified `BPPARAM` default.
- [x] Add tests for strict 100% and caller-selected 50% coverage filtering, partial-set scoring, duplicate members, duplicate names, malformed identifiers, and output ordering.

**Acceptance:** partial coverage behaves exactly as specified; invalid structures fail before native code/parallel work; no per-sample rematching exists.

### Session 3 - Correct missing-value and numerical semantics

- [x] Refactor native scoring to consume precomputed integer memberships rather than row names and character sets.
- [x] Omit `NA` and `NaN` per set/sample and recompute effective `n`.
- [x] Return `NA_real_` when fewer than two observed values remain in a cell.
- [x] Preserve explicit and implicit zeros in `n`, mean, deviation, and score.
- [x] Reject any infinity that escapes R validation and report useful context.
- [x] Replace the fixed `1e-9` zero clamp with a scale-aware numerical tolerance.
- [x] Detect materially negative or non-finite native results as internal errors.
- [x] Compare every deterministic fixture against the pure R oracle.

**Acceptance:** the canonical examples in Section 2.7 pass; `c(1, 2, NA)` equals scoring `c(1, 2)`; `c(4, NA)` returns `NA`; no missing value poisons an otherwise scoreable set.

### Session 4 - Preserve dense and sparse representations

- [x] Provide separate or safely dispatched dense and sparse native paths.
- [x] Ensure a base dense matrix remains dense and a sparse `Matrix` remains sparse until the bounded native work buffer.
- [x] Ensure sparse implicit zeros are counted as observed members.
- [x] Support stored `NA`/`NaN` entries in sparse matrices according to the same semantics as dense matrices.
- [x] Avoid whole-matrix dense copies in validation, negativity checks, infinity checks, subsetting, worker export, and native conversion.
- [x] Add dense/sparse equivalence tests across ordinary values, all-zero columns, partial coverage, stored missing entries, and reordered rows/columns.
- [x] Add a testable internal dispatch seam so CI can verify that sparse input reaches the sparse kernel without brittle source-text assertions.

**Acceptance:** dense and sparse outputs agree within a justified tolerance; sparse input never invokes the dense kernel; no `as.matrix()` call occurs on the full sparse input path.

### Session 5 - Make parallel execution deterministic and scalable

- [x] Profile the baseline one-task-per-column design and select a bounded chunking/iteration strategy compatible with BiocParallel.
- [x] Precompute memberships once and reuse them in every chunk.
- [x] Avoid capturing/serializing the full matrix once per column or task, especially under SOCK/Snow backends.
- [x] Preserve exact output order independent of task completion order.
- [x] Validate `BPPARAM` and propagate worker errors with chunk/sample context.
- [x] Avoid nested OpenMP/BiocParallel parallelism. Remove unused OpenMP flags unless a separate measured decision justifies them.
- [x] Test `SerialParam()` and a small portable multi-worker backend; conditionally test platform-specific backends only where supported.
- [x] Test one sample, fewer samples than workers, empty worker chunks, and backend reuse.

**Acceptance:** serial and parallel outputs/names are identical; task count is bounded rather than proportional to millions of cells; no unused OpenMP configuration remains.

### Session 6 - Complete the scientific invariant and adversarial test suite

- [x] Add deterministic randomized tests comparing native dense and sparse results with the pure R oracle.
- [x] Test non-negativity and upper bound by observed sum over broad non-negative inputs and set sizes.
- [x] Test member permutation, row permutation, sample permutation, irrelevant-row addition, set addition/removal, and sample addition/removal invariance.
- [x] Test positive homogeneity and constant-shift behavior.
- [x] Test overlapping gene sets without cross-set effects.
- [x] Test minimum sizes globally and per sample after missing-value omission.
- [x] Test extreme but finite values, signed zero, integer input, zero-only sparse matrices, unnamed columns, and empty/degenerate inputs.
- [x] Test that inputs are not mutated.
- [x] Test stable matrix dimensions and names after sets are omitted.
- [x] Run native memory/undefined-behavior tooling available in the current R toolchain where practical.

**Acceptance:** all specification statements in Section 2 have at least one direct test; randomized oracle comparisons find no dense/sparse/native discrepancy.

### Session 7 - Profile and optimize without changing semantics

- [x] Add a reproducible synthetic data generator and benchmark harness outside package tests.
- [x] Measure elapsed time, throughput, allocations/peak memory, output digest, package versions, hardware, and session information.
- [x] Benchmark at least:
  - dense bulk-like data, approximately 20,000 features, 60 samples, and 1,000 sets of size 20;
  - sparse single-cell-like data, approximately 20,000 features and configurable sample counts up to at least 600 in the default benchmark;
  - high-overlap and low-overlap set catalogues;
  - matrices containing zeros and sample-specific missing values;
  - serial and multi-worker execution.
- [x] Profile native hot spots before optimizing. Likely candidates are sparse random access, repeated set traversal, column copying, and worker serialization.
- [x] Consider a sparse-specific algorithm that accounts analytically for implicit-zero deviation rather than materializing dense columns, but adopt it only after oracle validation and measured benefit.
- [x] Keep performance changes isolated from API/scientific changes.
- [x] Add a small non-gating CI benchmark smoke test and keep full benchmarks manual or scheduled.

**Acceptance:** sparse peak memory scales with stored nonzeros plus output/bounded work buffers rather than the full feature-by-sample product; no documented efficiency claim lacks a reproducible benchmark.

### Session 8 - Finish package metadata, hygiene, and portability

- [x] Replace every placeholder in `DESCRIPTION` with accurate title, description, author, contact, GPL-3 license declaration, URLs, bug tracker, dependencies, and validated `biocViews` terms.
- [x] Remove the unnecessary `Rcpp::sourceCpp` import if package registration makes it redundant.
- [x] Regenerate native registration and namespace files.
- [x] Review `dev/` bootstrap scripts: delete obsolete scaffolding or clearly retain it as maintainer-only and exclude it from builds.
- [x] Expand `.Rbuildignore` for agent files, plans, local caches, benchmark outputs, and other non-package content while retaining user-facing source/docs.
- [x] Ensure no generated binary, check directory, tarball, cache, or platform-specific file is tracked.
- [x] Add package-level documentation, `NEWS.md`, and an appropriate citation entry for the thesis/package without inventing a DOI.
- [ ] Verify compiler portability on Linux, macOS, and Windows; avoid GNU-only assumptions unless guarded.

**Acceptance:** source tarball contains only intentional package content; install/check succeeds on all supported platforms; metadata contains no placeholders or false Bioconductor status.

### Session 9 - Write durable scientific and user documentation

- [x] Create a durable tracked scientific specification, for example `inst/SCIENTIFIC_SPEC.md`, containing the normative formula and value/coverage semantics so correctness does not depend on this temporary plan.
- [x] Expand roxygen documentation for `genefunnel()` and `gene_set_coverage()` with exact input classes, matching rules, equation, zeros, missing values, partial coverage, errors, return shape, and examples.
- [x] Build a vignette using current Bioconductor conventions that covers:
  - algorithm and interpretation;
  - canonical examples;
  - dense and sparse usage;
  - strict versus 50% caller-controlled coverage filtering;
  - serial and parallel execution;
  - exact identifier matching/versioning;
  - limitations from Section 2.8;
  - guidance that GeneFunnel performs no normalization or downstream testing.
- [x] Rewrite the README as a concise quick start and link to the vignette/specification.
- [x] State development/Bioconductor installation status accurately.
- [x] Keep examples small, deterministic, offline, and fast enough for package checks.
- [x] Remove or qualify any claim of “efficient sparse matrices”, superiority, absolute activity, or cross-dataset comparability unless directly supported by the final benchmark evidence.

**Acceptance:** a user with no thesis can understand and reproduce every core behavior from package docs alone; rendered documentation and examples pass checks.

### Session 10 - Add CI and current Bioconductor quality gates

- [x] Consult only current official R/Bioconductor guidance at implementation time; standards may have changed since this plan was written.
- [x] Add CI for supported R/Bioconductor combinations and major operating systems.
- [x] Run at minimum package install, unit tests, `R CMD check`, and `BiocCheck` where applicable.
- [x] Include a job that builds the source tarball and tests the installed tarball in a clean library.
- [x] Add native sanitizer/valgrind checks if supported and stable in CI.
- [x] Add documentation build/link checks and a small sparse-memory/dispatch smoke test.
- [x] Cache dependencies safely without committing caches.
- [x] Ensure tests do not require network access, large external data, a thesis file, or many cores.

**Acceptance:** CI is green from a clean checkout; failures are actionable; package checks produce no errors or warnings and no unexplained notes.

### Session 11 - Encode reproducibility and thesis-derived benchmark protocols

- [ ] Add versioned benchmark scripts and a manifest describing seeds, dimensions, set construction, preprocessing, methods, and recorded environment.
- [ ] Encode the thesis-inspired controlled scenarios without needing thesis text or data:
  - equal-value sets;
  - equal sums with low versus high within-set deviation;
  - all-zero sets;
  - one-nonzero sets;
  - sample-specific missingness;
  - partial matrix coverage;
  - sparse single-cell-like matrices;
  - sample/gene-set independence.
- [ ] Optionally add non-gating competitor benchmarks for GSVA Poisson/Gaussian, ssGSEA, PLAGE, and z-score only if dependencies and current APIs are pinned/documented. Keep competitor code outside core package tests.
- [ ] Do not assert the historical thesis runtimes as universal thresholds. Record observed values with hardware/software context.
- [ ] If real thesis datasets are supplied later, add optional download/checksum/provenance scripts rather than bundling them in the package.
- [ ] Produce a machine-readable result summary and human-readable benchmark report that can be regenerated from a clean environment.

**Acceptance:** scientific and computational claims can be audited from repository scripts; nothing claims that unavailable thesis data were rerun.

### Session 12 - Release-candidate audit

- [ ] Re-read the entire specification and trace every requirement to code, tests, and documentation.
- [ ] Run full checks from a clean clone/library on all CI platforms and the built source tarball.
- [ ] Re-run full dense/sparse/parallel benchmarks and archive result metadata.
- [ ] Audit warnings, messages, errors, numerical tolerances, native registration, package size, examples, and license/citation text.
- [ ] Verify no accidental API expansion or undocumented behavior exists.
- [ ] Disclose non-trivial AI assistance/provenance in the Bioconductor issue and PR, with provenance cited in new code, according to current policy.
- [ ] Set the release/submission version according to current Bioconductor policy only when all gates pass.
- [ ] Write release notes describing breaking behavior relative to the prototype: partial coverage, proper missing-value omission, sparse preservation, duplicate handling, and stricter validation.
- [ ] Create a tag/release only after maintainer review; do not claim Bioconductor availability before acceptance.

**Acceptance:** all Definition-of-Done items in Section 9 pass; `git status` is clean; the source tarball is independently installable and scientifically traceable.

## 6. Required test matrix

The final suite should make this table easy to trace to individual test files.

| Category | Required cases |
|---|---|
| Formula | equal values; all zeros; one non-zero; two unequal values; mixed ordinary values; comparison with pure R oracle |
| Missingness | one `NA`; multiple `NA`; `NaN`; fewer than two observed; different missing patterns by sample; stored sparse missing values |
| Zero semantics | explicit dense zero; implicit sparse zero; zero versus `NA` contrast; all-zero row/column |
| Coverage | 100%; 50%; arbitrary partial; zero/one matched member; unmatched members; set order after omission |
| Identifiers | exact matching; reordered rows; duplicate row names rejected; duplicate set names rejected; duplicate members deduplicated; empty/`NA` identifiers rejected |
| Input values | integer and double; negative rejected; `Inf`/`-Inf` rejected; extreme finite values; signed zero |
| Input classes | base dense; supported Matrix dense; supported Matrix sparse; unsupported/coercion-prone object rejected clearly |
| Independence | sample add/remove/reorder/modify; gene-set add/remove/reorder/modify; overlapping sets; irrelevant matrix rows |
| Parallelism | serial versus portable multi-worker backend; one sample; fewer columns than workers; deterministic ordering; worker error propagation |
| Numerics | non-negative within tolerance; score <= sum; exact zero cleanup; materially negative internal result not hidden |
| API | dimensions; row/column names; return is numeric matrix; inputs unmodified; coverage helper agrees with scorer |

## 7. Benchmark and performance acceptance

Performance is secondary to scientific correctness, but the package must fulfill its stated high-throughput/sparse intent.

### Hard requirements

- No whole sparse input is converted to a dense matrix.
- Membership matching is done once per call, not once per sample.
- Work scheduling is bounded/chunked rather than one task per cell for very large sample counts.
- Output allocation of `n_sets * n_samples` is expected; other memory should be bounded and explained.
- Dense and sparse kernels remain semantically identical under the oracle suite.
- Documentation reports measured context instead of universal speed/memory promises.

### Performance goals, not correctness gates

- Native implementation should substantially outperform the pure R oracle on representative data.
- Sparse execution should benefit from sparsity in memory and, after profiling, preferably runtime.
- Parallel execution should show useful scaling on sufficiently large inputs without regressing small inputs through overhead.
- Any optimization that changes floating-point operation order must remain within documented tolerance and preserve all invariants.

## 8. Documentation content Codex must preserve without the thesis

The final package documentation must explicitly convey these thesis-derived points:

1. GeneFunnel is functional class scoring: feature-by-sample input becomes gene-set-by-sample output.
2. It combines total observed activity with a scaled absolute-deviation penalty.
3. Scores are calculated independently for each sample and each gene set.
4. Zeros are retained because they represent measured lack of activity.
5. Missing values are omitted because they represent lack of information, reducing sample-specific set size.
6. A gene set requires at least two observed members in a sample.
7. Partial assay coverage is allowed; the caller chooses an acceptable coverage threshold using diagnostics.
8. Overlapping sets are valid.
9. Non-negative input is required and no automatic preprocessing occurs.
10. The score is non-negative for valid inputs, with equal values scoring their sum and maximally deviant one-nonzero sets scoring zero.
11. Sparse and parallel support are implementation requirements, not changes to the mathematical result.
12. The method’s assumptions and downstream statistical limitations must be stated, not obscured by benchmark claims.

## 9. Definition of done

### Scientific fidelity

- [x] Equation matches Section 2.2 exactly.
- [x] Zeros, `NA`/`NaN`, infinities, negatives, duplicates, and partial coverage match Sections 2.3-2.5 exactly.
- [x] Every property in Section 2.6 is directly tested.
- [x] Pure R, dense native, sparse native, serial, and parallel paths agree.

### Scalability

- [x] Sparse input is never fully densified.
- [x] Gene-set indices are precomputed once.
- [x] Parallel task count and memory are bounded.
- [x] Benchmark scripts and environment metadata are reproducible.

### Package quality

- [x] No placeholder metadata or tracked binaries.
- [x] Generated namespace/native files are current.
- [x] Complete tests, reference documentation, vignette, README, NEWS, citation, and durable scientific spec exist.
- [ ] Source tarball installs and passes current `R CMD check`/Bioconductor checks on supported platforms.
- [ ] CI is green from a clean checkout without network/data/thesis dependencies.

### Trustworthy claims

- [x] README/vignette claims are limited to tested behavior.
- [x] Cross-dataset comparability and downstream statistical suitability are qualified.
- [x] No exact thesis benchmark is claimed as reproduced unless corresponding data and scripts were actually run.
- [x] Release status and Bioconductor availability are stated accurately.

## 10. Explicit non-goals for this plan

- Changing or tuning the GeneFunnel equation.
- Introducing learned weights, gene-specific dynamic ranges, topology, ranks, cross-sample normalization, or cross-set normalization.
- Bundling GO, Reactome, KEGG, MSigDB, or other gene-set catalogues.
- Performing identifier conversion.
- Adding normalization, imputation, differential-expression testing, plotting, a Shiny app, or server-side analysis.
- Proving that the score distribution satisfies every downstream model.
- Reproducing the complete thesis analyses without the required external datasets.
- Adding broad container/infrastructure work unless needed for reproducible package benchmarks or release checks.

## 11. Maintainer decision triggers

Codex should continue autonomously under the defaults above. Stop and ask only when work would require one of these changes:

- alteration of the equation or deviation scaling;
- acceptance of negative values or automatic transformation;
- treating zeros as missing;
- treating missing cells as zero;
- requiring full gene-set coverage or adding a built-in coverage cutoff;
- changing the primary return type from a matrix;
- adding public scoring parameters;
- supporting a new large abstraction such as `DelayedArray`, `SummarizedExperiment`, or GPU execution before the core release is complete;
- publishing author identifiers, affiliations, citation identifiers, or Bioconductor status that are not verifiable from repository/official sources.

## 12. Handoff ledger template

At the end of each session, append a compact entry to `.agent/roadmap.md`:

```text
YYYY-MM-DD | session N | commit <sha>
Scope: <one-line deliverable>
Changed: <key files/interfaces>
Verified: <targeted tests + package checks + benchmark if relevant>
Decisions: <durable scientific/technical choices>
Remaining: <next unchecked plan items>
Risks/blockers: <only real unresolved items>
```

Do not copy large diffs, transient version numbers, or information already obvious from git history. The ledger exists to let a fresh Codex session resume accurately without the thesis or prior chat context.
