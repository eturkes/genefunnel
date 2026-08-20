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

- Adds a Lean 4.32.2 kernel proof for the executable exact-rational below-mean
  scorer. The proof establishes equation equality and exact score bounds for
  every modeled cell. It covers missingness, insufficient support, and
  negative-value behavior. The proof boundary explicitly excludes native
  binary64, R/FFI, storage, and parallel refinement claims.
- Adds deterministic formula, invariant, randomized-oracle, dense-sparse, and
  serial-parallel tests.
- Preserves exact-zero classification when positive subnormal means underflow
  on platforms whose C++ `long double` has the range of `double`.
- Preserves representable finite scores when an intermediate sum or weighted
  term would overflow on those platforms.
- Adds a reproducible synthetic benchmark harness with identity-checked output
  digests and environment metadata.
- Adds an explicit-version protocol dispatcher and SHA-256 index that freeze
  protocol 1.0.0's runner/manifest/fixture closure without an implicit latest
  alias or new dependency.
- Prospectively freezes design-only validation protocol F-2.0.0 with separate
  bulk RNA, pseudobulk RNA, and proteomics gates. It fixes exact Reactome and
  method bytes, APIs, same-input and native regimes, matched controls, and
  held-out replication. It also fixes minimum effects, multiplicity, null
  tests, exclusions, and RNG namespaces.
- Adds prospective amendment F-2.1.0 without changing the frozen parent. It
  defines cell-type-complete contrasts, dependence-preserving biological-unit
  weights, fixed-panel-only claim scope, and null-permutation eligibility. It
  also pins the Reactome 97 roster and validates assay-input schemas. It
  specifies linear-proteomics eligibility. Metadata-only selection precedes
  opaque-byte closure.
  The amendment validates the acquisition transcript through the first
  held-out access. It includes adversarial contract checks. Its 128
  rows yield a 266-row effective contract and a fixed 5,062-row target roster.
  The hermetic fixture covers 12 sources, 90 tasks, 60 null contrasts, 29
  objects, 37 F-A events, and 99 adversaries. Validation closes typed design
  fields, byte-exact UTF-8 scalars, and underflow-aware numeric parsing. It also
  closes total archive decoder and resource bounds. It closes finite
  acquisition bounds, dependence-group bijection, guarded descriptors, and
  exact S/B/I tree, runtime, and mode evidence. It states the trusted-process boundary. No
  candidate source, archive, or molecular table was selected or accessed.
  Execution requires exact F-S, F-B, F-E, and controlled-runner artifacts.
- Locks a dependency-free scaled numerical oracle and paired default-path
  performance protocol before component implementation.
- Freezes the internal observed-member sensitivity schema, dependency-free
  exact dyadic deletion oracle, and profile-before-optimization rule. It also freezes the
  controlled feature-loss and measurement-repeat design, held-out models,
  bootstrap, and rejection gates before package diagnostics existed.
- Adds an unexported exact brute-force observed-member sensitivity prototype
  with compact aligned summaries, explicit representability status, exact
  canonical ties, bounded matrix iteration, and dense/sparse serial/SOCK
  identity. Public promotion is governed by the frozen empirical protocol.
- Byte-pins the previously implicit sensitivity-profile fixture and measurement
  passes before profiling. The protocol fixes isolated installation, output
  identity, CPU and allocation attribution, environment capture, and an
  optimization-only decision runner.
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
- Records the complete negative controlled sensitivity result. Both
  feature-loss and measurement-repeat prediction gates failed by wide margins. Retains
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
- Records the complete 124,416-measurement controlled result. Every frozen
  synthetic gate passed. However, severe-dropout diagnostics remain comparable
  to the planted complementarity effect and preclude a dropout-stability claim.
- Freezes CellBench and Kang external-data execution details before any external
  endpoint. These details include zero-abundance ranking, condition-level gate
  units, duplicated barcode joins, cell-type eligibility, and exact sign tests.
- Adds a fail-closed CellBench runner with verified inputs, CEL-seq2-only set
  selection, isolated package installation, complete fixed-grid error and
  split/cross-platform stability evidence, and retained scientific failures.
- Records the negative CellBench result. Both co-primary gates failed because
  pair-set errors were large. Exact-zero measured scores also left much of the
  fixed condition grid undefined. High complete-case correlations do not
  rescue the pre-specified decision.
- Adds a fail-closed Kang/Reactome runner with exact sparse-matrix and barcode
  alignment. The runner fixes cell-type splits, raw-UMI pseudobulk profiles,
  fixed 16-group technical stability, held-out donor directions, and exact sign
  tests. It retains complete scientific-failure evidence.
- Records the Kang result. Technical split stability passed, but interferon
  gamma repeated in only two of four held-out donors. Neither primary pathway
  passed its Holm-adjusted exact donor sign test. The perturbation and
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
  sign, compact summaries, and algebraic properties. It also fixes the
  non-causal boundary before any public reliability API or empirical diagnostic
  result.
- Proves the unknown-member boundary. Non-negativity leaves every absent value
  unidentified. One missing member can yield a finite score plateau. Two or
  more missing members can make the score unbounded. Finite caller caps give
  sharp score bounds and conservative deletion-delta enclosures without
  probabilistic meaning.
- Adds an executable BiocStyle vignette with canonical examples,
  caller-controlled coverage policies, dense and sparse matrices, and parallel
  execution. It also covers
  thesis-derived benchmark concepts and evidence boundaries across controlled,
  comparative, empirical, and resource designs.
- Expands function help and the README with the complete public contract and
  accurate development status.
