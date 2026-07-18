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

## Documentation

- Adds a durable scientific specification for the equation, value and coverage
  semantics, invariants, computational equivalence, and limitations.
- Records the proved magnitude/balance factorization, prior-art boundary, and
  component API semantics without changing the primary scorer.
- Records the weighted aggregation-gap theorem, equality condition, normalized
  discrepancy, eligibility/missingness policy, and prior-art boundary before
  any group-level API is implemented.
- Adds an executable BiocStyle vignette covering canonical examples,
  caller-controlled coverage, dense/sparse matrices, and parallel execution.
- Expands function help and the README with the complete public contract and
  accurate development status.
