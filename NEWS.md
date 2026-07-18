<!-- Assisted-by: OpenAI Codex. -->

# genefunnel 0.99.0

## New features

- Implements the GeneFunnel score for dense and sparse non-negative matrices,
  with sample-specific omission of `NA` and `NaN` values.
- Adds `gene_set_coverage()` for exact-match coverage diagnostics and
  caller-controlled filtering.
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
- Adds a versioned controlled-science protocol covering analytic score cases,
  sparse/dense identity, partial coverage, missingness, and independence, with
  machine-readable results and regenerated Markdown reports.

## Documentation

- Adds a durable scientific specification for the equation, value and coverage
  semantics, invariants, computational equivalence, and limitations.
- Records the proved magnitude/balance factorization, prior-art boundary, and
  pre-implementation component semantics without changing the public API.
- Adds an executable BiocStyle vignette covering canonical examples,
  caller-controlled coverage, dense/sparse matrices, and parallel execution.
- Expands function help and the README with the complete public contract and
  accurate development status.
