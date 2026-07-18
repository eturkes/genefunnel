<!-- Assisted-by: OpenAI Codex. -->

# genefunnel 0.99.0

## New features

- Implements the GeneFunnel score for dense and sparse non-negative matrices,
  with sample-specific omission of `NA` and `NaN` values.
- Adds `gene_set_coverage()` for exact-match coverage diagnostics and
  caller-controlled filtering.
- Adds `genefunnel_components()` for aligned score, observed-sum, penalty,
  balance, effective-size, observed-fraction, status, and scaled numerical
  diagnostics without changing `genefunnel()`.
- Supports deterministic, bounded parallel scoring through BiocParallel.

## User-visible changes from the prototype

- Partially covered gene sets are scored when at least two unique members
  match; sets below that minimum are omitted with one aggregate warning.
- Duplicate members are counted for diagnostics and deduplicated for scoring.
- Invalid structures, negative values, and infinities are rejected before
  parallel or native work begins.
- S3 subclasses that spoof matrices, gene sets, identifiers, or parallel
  parameters are rejected before their methods can alter validation or scoring.
- Dense inputs remain dense and sparse inputs remain sparse through bounded
  scoring chunks.

## Quality

- Adds deterministic formula, invariant, randomized-oracle, dense-sparse, and
  serial-parallel tests.
- Preserves exact-zero classification when positive subnormal means underflow
  on platforms whose C++ `long double` has the range of `double`.
- Preserves representable finite scores when an intermediate sum or weighted
  term would overflow on those platforms.
- Adds a reproducible synthetic benchmark harness with identity-checked output
  digests and environment metadata.
- Locks a dependency-free scaled numerical oracle and paired default-path
  performance protocol before component implementation.
- Freezes the internal observed-member sensitivity schema, dependency-free
  exact dyadic deletion oracle, profile-before-optimization rule, controlled
  feature-loss/measurement-repeat design, held-out models, bootstrap, and
  rejection gates before package diagnostics existed.
- Adds an unexported exact brute-force observed-member sensitivity prototype
  with compact aligned summaries, explicit representability status, exact
  canonical ties, bounded matrix iteration, and dense/sparse serial/SOCK
  identity. Public promotion is governed by the frozen empirical protocol.
- Byte-pins the previously implicit sensitivity-profile fixture and measurement
  passes before profiling, with isolated clean-SHA installation, output identity,
  CPU/allocation attribution, environment capture, and an optimization-only
  decision runner.
- Records the fixed exact-brute profile: median call 213.585 seconds and
  exact-arithmetic stack share 0.999403 both cross their frozen optimization
  triggers. The result makes no performance, reliability, or public-API claim.
- Adopts an exact sorted-prefix sensitivity path while retaining the brute
  oracle. Fixed/extreme/randomized delta objects, the representation/backend
  suite, and the clean frozen-workload digest remain identical.
- Byte-pins controlled sensitivity execution details - R draw calls, row order,
  predictor encoding/scaling, bootstrap order, isolated installation,
  checkpoints, and evidence artifacts - before constructing the fixed grid.
- Adds deterministic controlled sensitivity profile/count/mask generation and
  fixed-schema feature-loss/measurement-repeat observations, with authoritative
  package scores, exact partial-input diagnostics, and fail-closed invariants.
- Adds frozen ten-fold baseline/augmented models, scenario-cluster bootstrap,
  co-primary endpoint/strata summaries, and a clean-archive resumable controlled
  runner that distinguishes scientific failure from execution failure.
- Records the complete negative controlled sensitivity result: both feature-
  loss and measurement-repeat prediction gates failed by wide margins. Retains
  deterministic study-dependent thinning curves as a failure envelope and
  selects the internal-oracle/no-public-API fallback.
- Implements dependency-free scaled double-double native diagnostics for
  overflow, cancellation, and underflow cases while keeping the ordinary score
  authoritative.
- Adds a versioned controlled-science protocol covering analytic score cases,
  sparse/dense identity, partial coverage, missingness, and independence, with
  machine-readable results and regenerated Markdown reports.
- Adds an unexported, fail-closed aggregation-audit prototype with explicit
  weights, support/missingness facts, identity residuals, and independent
  randomized-oracle coverage; empirical promotion gates remain open.
- Freezes the controlled audit runner's shared-latent RNG, measurement/dropout,
  cross-validation, model, quantile, and bootstrap execution rules before any
  synthetic or downloaded-data result.
- Adds a clean-commit, isolated-install synthetic audit runner with deterministic
  fork-safe generation, resumable identity-checked checkpoints, paired model/
  bootstrap summaries, full failed-endpoint retention, and artifact hashes.
- Records the complete 124,416-measurement controlled result: every frozen
  synthetic gate passed, while severe-dropout diagnostics remain comparable to
  the planted complementarity effect and preclude a robustness claim.
- Freezes CellBench and Kang external-data execution details - including
  zero-abundance ranking, condition-level gate units, duplicated barcode joins,
  cell-type eligibility, and exact sign tests - before any external endpoint.
- Adds a fail-closed CellBench runner with verified inputs, CEL-seq2-only set
  selection, isolated package installation, complete fixed-grid error and
  split/cross-platform stability evidence, and retained scientific failures.
- Records the negative CellBench result: both co-primary gates failed because
  pair-set errors were large and exact-zero measured scores left much of the
  fixed condition grid undefined; high complete-case correlations do not
  rescue the pre-specified decision.
- Adds a fail-closed Kang/Reactome runner with exact sparse-matrix/barcode
  alignment, deterministic cell-type splits, raw-UMI pseudobulk profiles,
  fixed 16-group technical stability, held-out donor directions, exact sign
  tests, and complete scientific-failure evidence.
- Records the Kang result: technical split stability passed, but interferon
  gamma repeated in only two of four held-out donors and neither primary
  pathway passed its Holm-adjusted exact donor sign test. The perturbation and
  biological-effect claims therefore fail.

## Documentation

- Adds a durable scientific specification for the equation, value and coverage
  semantics, invariants, computational equivalence, and limitations.
- Records the proved magnitude/balance factorization, prior-art boundary, and
  component API semantics without changing the primary scorer.
- Records the weighted aggregation-gap theorem, equality condition, normalized
  discrepancy, eligibility/missingness policy, and prior-art boundary before
  any group-level API is implemented.
- Fixes the exact observed-member deletion sensitivity, canonical support,
  sign, compact summaries, algebraic properties, and non-causal boundary before
  any public reliability API or empirical diagnostic result exists.
- Adds an executable BiocStyle vignette covering canonical examples,
  caller-controlled coverage, dense/sparse matrices, and parallel execution.
- Expands function help and the README with the complete public contract and
  accurate development status.
