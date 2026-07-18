<!-- Assisted-by: OpenAI Codex. -->

# GeneFunnel research and scale plan

**Observed baseline:** `main` at `c789a23` on 2026-07-17; package `0.99.0`;
two public functions; pre-submission Bioconductor candidate.

**Mission:** turn one trustworthy score into an auditable gene-set signal
system: explain the score, expose aggregation effects, measure fragility, and
scale to partitioned assays without weakening the original method.

**Plan type:** forward research portfolio. The completed fidelity/release plan
remains recoverable from Git at `c789a23:PLAN.md`.

## 1. Authority and operating model

- [`inst/SCIENTIFIC_SPEC.md`](inst/SCIENTIFIC_SPEC.md) = normative scoring
  contract. It owns the equation, accepted values, missingness, coverage,
  identifiers, invariants, interpretation, and limitations.
- [`inst/AGGREGATION_SPEC.md`](inst/AGGREGATION_SPEC.md) = Workstream B's
  theorem/design contract. It owns the pre-API estimand, eligibility,
  missingness, terminology, and novelty boundary.
- `PLAN.md` = hypotheses, future interfaces, sequencing, evidence gates, and
  kill criteria. Unchecked boxes are prospective work, not release promises.
- [`.agent/roadmap.md`](.agent/roadmap.md) = chronological execution evidence
  and handoffs. Git history preserves the completed release programme.
- Every session starts by reconciling this plan with live `HEAD`, current
  official dependencies, and the latest evidence. Each cohesive result ends
  with targeted verification, the broadest practical package checks, a compact
  roadmap entry, and one scoped commit.
- Exploration may fail. Negative results close a hypothesis honestly; they do
  not justify moving its success criteria.
- Before an experiment starts, a committed/versioned protocol fixes its data or
  workload, primary/secondary endpoints, minimum effect, unit of analysis,
  repeats, uncertainty rule, and multiplicity rule where relevant. A threshold
  changes only with recorded rationale before results are inspected.

## 2. Frozen core

The present scorer is the control condition for every extension.

1. `genefunnel(mat, gene_sets, BPPARAM)` keeps its signature, exact equation,
   value/missingness semantics, ordering, warnings, and plain numeric matrix
   return.
2. Existing list input stays canonical. Dense/sparse and serial/parallel paths
   remain equivalent under the independent oracle and all scientific
   invariants.
3. Zero remains observed; `NA`/`NaN` remains absent; unmatched features remain
   a coverage fact. GeneFunnel adds no automatic normalization, imputation,
   mapping, catalogue retrieval, or coverage threshold.
4. New facts use additive functions or explicit adapters. Rich containers never
   enter through implicit coercion, and the primary scorer never changes return
   class according to input size.
5. The scientific specification continues to govern the frozen scorer and
   coverage helper. Every additive public API receives its own explicit
   contract, or moves to a companion package when it changes package identity.
6. Algebraic diagnostics are descriptive, not causal gene attribution,
   inferential uncertainty, or proof of biological activity.
7. A scale-free component alone does not make separately processed data sets
   comparable. Units, preprocessing, detection, set definition, and coverage
   remain part of the estimand.
8. Signed inputs, learned weights, topology, ranks, or a changed penalty define
   a different method and require a separately named API or companion package.
9. Provenance and diagnostic sidecars never attach to or change the core matrix
   return implicitly.

## 3. Flagship mathematical opportunity

Let an observed non-negative gene-set/sample vector be
\(x=(x_1,\ldots,x_n)\), \(n\ge2\), with total \(T=\sum_i x_i\). When
\(T>0\), define shares \(p_i=x_i/T\), the uniform share vector
\(u_i=1/n\), and
\(\operatorname{TV}(p,u)=\frac12\sum_i|p_i-u_i|\). The existing score
\(F\) factorizes exactly:

\[
F(x)
=T-\frac{n}{2(n-1)}\sum_i|x_i-T/n|
=T\left[1-\frac{n}{n-1}\operatorname{TV}(p,u)\right]
=T B(x).
\]

This exposes four aligned facts already latent in the kernel:

- **magnitude** \(T\): observed total signal;
- **balance** \(B\in[0,1]\): evenness of within-set signal shares;
- **arithmetic penalty** \(T-F\), equal to \(T(1-B)\) when \(T>0\); and
- **effective size** \(n\): sample-specific observed member count.

`Balance` is the working field term; `coherence` is avoided because it commonly
means cross-sample correlation. When \(T=0\), \(F=0\) but \(p\) and \(B\) are
undefined: report `NA`, not “perfect balance.” Mathematically, \(B\) is exactly
[Bulla's evenness index](https://doi.org/10.2307/3545713) and the normalized
complement of the Pietra/Hoover quantity. Neither the index nor its
normalization is novel. Any novelty must concern this exact GeneFunnel
factorization, implementation, and validated biological use. Balance also
depends on effective size, retained member identity, gene-specific abundance/
dynamic range, and measurement noise; it is not automatically comparable
across those changes.

The score is symmetric, positively homogeneous, and concave. For compatible
vectors \(x_j\), non-negative weights \(w_j\) summing to one, projection
\(P=I-\mathbf1\mathbf1^\top/n\), and \(c=n/[2(n-1)]\):

\[
J
=F\!\left(\sum_jw_jx_j\right)-\sum_jw_jF(x_j)
=c\left[\sum_jw_j\lVert Px_j\rVert_1
-\left\lVert P\sum_jw_jx_j\right\rVert_1\right]\ge0.
\]

`J` is the aggregation gap: centered member deviations can cancel during
aggregation, increasing the aggregate score above the weighted unit score. For
example, `c(4, 0)` and `c(0, 4)` each score zero; their mean `c(2, 2)` scores
four. This is an exact
computational fact, not yet a biological heterogeneity measure. The theorem
requires common member support and compatible units; sample-varying missingness
changes the dimension and invalidates the direct comparison. Proportional
profiles are sufficient but not necessary for \(J=0\): equality also holds when
centered coordinate deviations do not oppose in sign. The gap therefore measures
a specific cancellation pattern, not generic heterogeneity.

### Research hypotheses

- **H1 - decomposition:** magnitude, balance, penalty, and effective size make
  score changes more interpretable and flag arithmetic patterns requiring
  confounding investigation.
- **H2 - aggregation:** the aggregation gap quantifies complementary member
  patterns and warns when aggregate-then-score disagrees with
  score-then-aggregate.
- **H3 - reliability:** deterministic deletion/thinning sensitivity predicts
  score error under held-out feature loss and technical-repeat instability.
- **H4 - catalogue reuse:** compiled catalogues amortize repeated matching and
  adjacency construction without changing any score.
- **H5 - block scale:** block source/sink execution bounds scoring buffers
  without changing any score.

## 4. Release lane - protect the candidate first

Exploratory notebooks/benchmarks may proceed, but any DESCRIPTION/NAMESPACE
change, optional or required dependency, S4 registration, public export,
version bump, or contract change waits for an explicit maintainer decision on
the pre-submission candidate scope.

- [x] Record green cross-platform CI at exact innovation baseline `c789a23`:
  [Actions run `29461199816`](https://github.com/eturkes/genefunnel/actions/runs/29461199816)
  passed all seven jobs after the S3-spoofing fix.
- [x] Maintainer intentionally reopened the candidate around new public work on
  2026-07-17.
- [ ] Disclose non-trivial `Assisted-by: OpenAI Codex` provenance in the
  Bioconductor issue and pull request after the reopened scope stabilizes.
- [ ] Create a release/tag only after maintainer review; claim Bioconductor
  availability only after acceptance.

The latter two submission gates are deferred, not waived, while the reopened
scope proceeds through the evidence gates below.

## 5. Workstream A - transparent components

**Outcome:** users can distinguish “more measured signal” from “more evenly
distributed signal” without changing the GeneFunnel score.

### A1 - theorem and semantics

- [x] Complete a prior-art audit covering normalized total variation,
  Hoover/Pietra-style indices, mean absolute deviation, pathway variability,
  and existing gene-set decomposition APIs.
- [x] Add a proof note and independent R oracle for the factorization, bounds,
  homogeneity, concavity, and zero edge case.
- [x] Ratify, test, and document these plan-level component semantics before
  native work:
  - `effective_size` and `observed_sum` are factual even when fewer than two
    members remain;
  - score, penalty, and balance are `NA` when effective size is below two;
  - score and penalty are zero but balance is `NA` for an all-zero scoreable
    cell; and
  - sample-specific observed fraction uses globally matched size as its
    denominator and remains distinct from declared-set coverage.
- [x] Define a representability contract for extreme finite inputs. A score can
  fit in an R double even when its mathematical observed sum/penalty does not;
  representable sum and penalty can also cancel, and a true nonzero balance can
  underflow. The score remains authoritative. Diagnostics must use explicit
  scaled/unavailable/ill-conditioned status, never silently return infinity or
  make the ordinary score path fail.
- [x] Use the exact retained-row universe, aggregate warning, set/sample order,
  and names of `genefunnel()`. Globally unscoreable sets remain available
  through `gene_set_coverage()` rather than appearing as diagnostic-only rows.
- [x] Document that balance measures evenness, not co-expression, regulation,
  pathway truth, or preprocessing invariance.
- [x] Define committed interpretation cases with known mass, balance, and
  support manipulations that the scalar score alone cannot distinguish.

### A2 - one-pass prototype

- [x] Freeze protocol `A2-1.0.0` before implementation: exact result/status
  shape, numerical safe region and tolerance, extreme fixtures, runtime/
  allocation/RSS workloads, paired repeat schedule, estimators, uncertainty,
  and rejection rules.
- [x] Add a dependency-free, platform-independent scaled double-double oracle
  covering total/penalty overflow, unsafe cancellation, balance underflow,
  subnormals, and semantic edge states.
- [x] Prototype `genefunnel_components()` as an additive API returning a named
  list of aligned base matrices for ordinary values plus an aligned, documented
  representation/status for extreme components.
- [x] Extend dense and sparse kernels to emit requested accumulators in one pass;
  the ordinary `genefunnel()` path allocates no diagnostic outputs.
- [x] Verify both identities in a higher-precision/scaled oracle. Require direct
  binary64 reconstruction only in a numerically safe region defined before
  tests; explicitly cover total overflow, subtractive cancellation, and balance
  underflow outside it.
- [x] Cover missingness, partial coverage, all-zero sets, one-nonzero sets,
  signed zero, subnormals/overflow, member permutations, and
  dense/sparse/serial/parallel identity.
- [ ] Measure default-path runtime, allocation, and RSS against a freshly rerun
  baseline. Fix dense/sparse workloads, paired/interleaved repeat count,
  estimator, environment, and confidence/equivalence rule before
  implementation; an upper confidence bound above 5% regression blocks changes
  to the default path.

**Go:** diagnostics match the authoritative score, add no second input pass,
pass the committed interpretation cases, and satisfy the numerical status and
performance contracts.
**Fallback:** if held-out benchmarks show no value beyond score, sum, and
effective size, keep the algebra in the scientific specification and keep the
new API internal.

## 6. Workstream B - aggregation audit

**Outcome:** quantify when weighted aggregation raises GeneFunnel score through
cancellation of different member patterns.

- [x] Prove and property-test the weighted Jensen formula independently of the
  native implementation.
- [x] Audit prior art on Jensen gaps for inequality/diversity indices,
  within-versus-pooled diversity decomposition, mixture/Simpson effects, and
  pseudobulk pathway scoring before making a novelty claim.
- [x] Freeze the primary estimand as a weighted mean
  \(\bar{x}=\sum_jw_jx_j\), \(w_j\ge0\), \(\sum_jw_j=1\), after common
  preprocessing. Homogeneity relates a physical non-negative sum to a scaled
  mean, but preprocessing need not commute with physical pooling. Ignore
  zero-weight units explicitly.
- [x] Specify eligibility before an API: non-negative weights, compatible units
  and preprocessing, identical ordered members, and common observed support.
- [x] For missing data, either reject the group or explicitly recompute all
  terms on a reported common-support intersection. Report removed member
  identities, reject intersections below two, and never turn missing cells into
  zero. Both unit and aggregate scores use the same intersection.
- [x] Define normalized gap as \(J/F(\bar{x})\in[0,1]\) when the aggregate score
  is positive and `NA` otherwise; verify whether that normalization answers the
  intended task before freezing it publicly.
- [x] Prototype a group-level audit returning aggregate score, weighted unit
  score, absolute gap, normalized gap where defined, removed member identities/
  count, and eligibility reasons.
- [x] Verify zero gap for proportional member profiles, positive gap for planted
  complementary profiles, the broader no-opposing-deviation equality case,
  permutation invariance, and non-negativity within numerical tolerance.
- [ ] Run a pre-specified factorial experiment over mixture proportion, library
  depth, dropout, gene-set size, overlap, outliers, weights, and aggregation
  level. Vary gene baselines, dynamic ranges, variances, and correlations;
  separate depth effects from member-pattern cancellation; stratify complexes,
  cascades, regulatory signatures, and mixed-direction sets.
- [ ] Validate first on known artificial mixtures or purified populations, then
  on perturbation data with biological replicates and donor-level inference;
  never treat cells as independent donors.

**Go:** gap meets the protocol's pre-specified curve-error,
replicate-stability, and incremental-effect thresholds under controlled
complementary mixtures, beyond depth, detected-feature count, and aggregate
score alone.
**Fallback:** if depth/dropout drives the gap as strongly as planted
complementarity, or the replication threshold fails, ship only the theorem and
a warning about aggregate-then-score behavior.

## 7. Workstream C - compiled catalogues and provenance

**Outcome:** repeated batches reuse exact matching and sparse adjacency safely;
every result can identify the catalogue and feature universe that produced it.

- [x] Design a logically immutable compiled object containing canonical
  members, exact feature order, integer memberships, coverage facts,
  row-to-set adjacency, schema/formula version, and deterministic content
  fingerprints. Expose no mutators; copy canonical state and fully validity-check
  every use.
- [x] Reject stale reuse on any feature identity/order mismatch. A fingerprint
  supports content identity and quick mismatch detection; it never replaces
  exact safety checks.
- [x] Keep the named-list API authoritative. Prototype compiled scoring behind a
  separate entry point; the frozen `genefunnel()` never accepts this object.
- [ ] Make serialization platform-independent and versioned. Research the hash
  implementation, maintenance, licensing, and supported platforms before
  selecting a dependency. Canonical encoding covers Unicode/encoding, names and
  order, integer width/endianness, missing values, signed zero, serialization
  version, and streamed result order. Fingerprints identify encoded content;
  they neither prove origin nor act as security signatures.
- [ ] Add an optional provenance sidecar with package/formula version,
  canonical catalogue fingerprint, ordered-feature fingerprint, caller-supplied
  catalogue/version and preprocessing labels, coverage-policy label, backend,
  and result digest. Unknown metadata stays unknown; the sidecar never changes
  or attributes origin to the core matrix implicitly.
- [x] Separate a portable R-side adjacency representation from ephemeral native
  per-process caches. Define SOCK worker initialization/reuse and never
  serialize a native pointer.
- [x] Reject malformed/tampered objects before scoring and verify portable
  serialization plus fresh/reused SOCK workers.
- [ ] Benchmark compilation, validation, serialization, worker initialization,
  warm/cold scoring, cumulative repeated-call time, object size, and RSS
  separately with paired/interleaved repeats and uncertainty.
- [x] Before implementation, ratify or replace with rationale the planning
  default: compilation amortizes by call five and the lower confidence bound on
  end-to-end improvement reaches 15%. Set a separate byte-per-membership budget
  for object/RSS growth.

**Go:** outputs/errors match the list path and the pre-specified cumulative-time,
amortization, serialization, and memory budgets pass.
**Fallback:** if matching/adjacency misses its minimum useful effect, retain an
internal prepared representation and ship only the provenance helper.

## 8. Workstream D - blockwise and ecosystem execution

**Outcome:** score inputs and outputs larger than memory through explicit,
transactional adapters.

The preferred ecosystem is Bioconductor-native: [DelayedArray block
processing](https://www.bioconductor.org/packages/devel/bioc/manuals/DelayedArray/man/DelayedArray.pdf)
already supports gridded sparse-aware parallel traversal;
[SparseArray](https://bioconductor.org/packages/devel/bioc/html/SparseArray.html)
provides current multidimensional sparse containers; and
[HDF5Array](https://bioconductor.org/packages/devel/bioc/html/HDF5Array.html)
provides dense/sparse on-disk DelayedArray backends. A later explicit
[SummarizedExperiment](https://bioconductor.org/packages/devel/bioc/html/SummarizedExperiment.html)
adapter can preserve assay/sample metadata. These are prospective optional
dependencies, not current requirements; HDF5 is the reference prototype, not a
pre-approved shipped dependency.

- [ ] Generalize the current column iterator into internal block-source and
  block-sink contracts with explicit viewport IDs and deterministic assembly.
- [ ] Preserve the current validation-before-scoring/error-order contract for
  `genefunnel()`. Only the transactional file-backed API may fuse value
  validation with block traversal; structural names/dimensions validate first
  and native guards remain fail-closed.
- [ ] Begin with all-feature, bounded-column blocks: the current score is not
  additive across feature-row blocks. Arbitrary 2-D/row tiling requires a
  separately proved multi-pass reduction and stays out of the first adapter.
- [ ] Read dense and sparse blocks without realizing the whole input. Select
  column width from backend layout and measured I/O. Convert bounded
  `SparseArray` blocks explicitly to the current `Matrix` boundary or implement
  a validated native bridge; defaults must never realize sparse blocks densely.
- [ ] Use exactly one parallel scheduler. Bound in-flight blocks and the reorder
  buffer; when an outer DelayedArray traversal owns `BPPARAM`, inner scoring is
  serial. Account for worker and callback-result memory multipliers.
- [ ] Add an explicit file-backed scoring function with a same-directory
  temporary sink, completeness/commit manifest, validated close, single-writer
  HDF5 policy, and rename/commit recovery tested per supported OS. Promise
  atomic visibility only where the filesystem guarantees it; never overwrite
  an existing destination implicitly.
- [ ] Keep in-memory and file-backed returns in different entry points. Add an
  explicit `SummarizedExperiment` adapter only after the matrix/on-disk contract
  is stable; assay selection must be named and preprocessing stays external.
  Decide between an explicit optional function and an S4 method, including its
  import/conditional-registration policy, before touching NAMESPACE.
- [ ] Prove block-partition, completion-order, dense/sparse, and backend
  invariance against the oracle. Require within-environment paired digest
  identity and cross-platform numerical equivalence, not universal byte
  identity.
- [ ] Cover implicit sparse zeros, stored explicit/signed zeros, stored
  `NA`/`NaN`, backend-specific sparse extraction, interrupted writes, and fresh
  worker processes.
- [ ] Benchmark end-to-end I/O, validation, scoring, writes, cold/warm workers,
  allocation, and process RSS. Require scoring buffers bounded by
  `O(features * block_columns + sets * block_columns + catalogue)` (or stored
  sparse payload) and account separately for `O(samples)` names/viewports/
  completeness metadata and the dense logical output.
- [ ] Before measurement, fix a larger-than-memory workload, maximum wall time,
  throughput floor, RSS/buffer budget, block widths, and repeats/uncertainty.
- [ ] Use single-cell-shaped synthetic matrices only as computational stress
  tests until Workstream F validates whether zero semantics fit cell-level use.

**Go:** no whole-input/output realization, transactional recovery works, results
match the in-memory path, and the pre-specified larger-than-memory time/memory
gates pass.
**Fallback:** reject any adapter that silently densifies, makes
return type dynamic, or cannot preserve score identity. Keep backends beyond
the HDF5 reference (for example Zarr or TileDB) out until a concrete workload
justifies each maintenance surface.

## 9. Workstream E - reliability and member sensitivity

**Outcome:** report how dependent a score is on particular measured members or
assay coverage, without “correcting” the result.

- [ ] Define exact leave-one-observed-member-out deltas for gene-set/sample
  cells with at least three observed members; call them sensitivity, never
  causal contribution.
- [ ] Prototype compact summaries first: largest absolute delta and member,
  signed delta, delta normalized by observed sum where defined, median absolute
  delta, and effective size. Break largest-delta ties by canonical member order.
  Full member-by-set-by-sample arrays require an explicit sink.
- [ ] Add deterministic/seeded coverage-thinning curves stratified by member
  abundance and detection rate. Distinguish global feature absence, sampled
  deletion, and sample-specific missingness; label assay-derived strata as
  study-composition dependent.
- [ ] Validate an optimized algorithm against brute-force recomputation before
  adopting it; profile before pursuing sorted-prefix or native acceleration.
- [ ] Test whether diagnostics predict held-out full-data score error and
  technical-replicate instability beyond coverage, set size, observed sum, and
  effective size. Characterize biological-replicate variation separately; it
  can be real signal rather than error.
- [ ] Keep unknown-member uncertainty honest: with only non-negativity, an
  unobserved value is unbounded. Caller-supplied limits support deterministic
  identification/sensitivity bounds; probabilistic intervals require a declared
  sampling model and belong to a later inferential project.
- [ ] Pre-specify held-out split, reliability endpoint, baseline model, minimum
  incremental effect, repeats, and uncertainty before calculating diagnostics.

**Go:** summaries are reproducible, representation-invariant, and pass the
pre-specified held-out incremental-effect and technical-repeat thresholds.
**Fallback:** if they do not outperform simple coverage/effective-size facts,
retain targeted adversarial tests and omit the public API.

## 10. Workstream F - external scientific validation

**Outcome:** establish where GeneFunnel, its components, and its aggregation
diagnostics are useful - and where they fail.

Recent benchmarks show that single-sample scores can vary with missing genes,
gene-set size, sample size, and cross-study integration choices
([Toro-Domínguez et al.,
2025](https://academic.oup.com/bib/article/26/6/bbaf684/8382580)); cell-level
methods need context-specific tests for detection-rate/dropout bias ([Noureen et
al., 2022](https://elifesciences.org/articles/71994); [Wang and Thakar,
2024](https://academic.oup.com/nargab/article/6/3/lqae124/7770961)). Sparse
speed alone is not scientific novelty: [PLAID](https://academic.oup.com/bioinformatics/article/41/12/btaf621/8321923)
already demonstrates large sparse catalogue/sample workloads. GeneFunnel's
claim must therefore rest on its estimand, diagnostics, and validated domain.

- [ ] Add a version-selecting harness plus frozen runner/manifest/fixtures for
  executable protocol `1.0.0`; introduce `2.0.0` only when new methods,
  scenarios, or assertions land.
- [ ] Pre-specify in that committed protocol target assays, aggregation level,
  primary/secondary hypotheses and endpoints, minimum effects, independent
  replication, exclusions, multiplicity, unit of inference, preprocessing,
  gene-set versions, coverage policy, competitor versions/APIs, and seeds
  before comparative results are inspected.
- [ ] Start with bulk/pseudobulk transcriptomics and proteomics whose
  non-negative meaningful-zero interpretation can be defended. Treat direct
  single-cell use as a separate hypothesis because technical zeros can violate
  the core semantic assumption.
- [ ] Build factorial simulations that vary magnitude and balance independently,
  then add depth, dropout, partial coverage, missing cells, set size, overlap,
  contamination, outliers, and complementary mixtures. Vary member baselines,
  dynamic ranges, variances, correlations, and biologically distinct set
  archetypes.
- [ ] Compare GeneFunnel with its sum/mean ablations and a small justified set of
  current methods in two regimes where feasible: the same compatible input to
  isolate scoring rules, and each method's documented end-to-end preprocessing
  pipeline to assess use. Competitor agreement is not ground truth.
- [ ] Use licensed external perturbation/negative-control/replicate data with
  immutable IDs, URLs, checksums, retrieval dates, and external storage; include
  known artificial mixtures or purified populations. Split protocol development
  from held-out evaluation.
- [ ] Measure perturbation ranking, technical-replicate stability, sensitivity
  to missing features, donor-level generalization, and resource cost. For
  negative controls, pre-specify the downstream contrast/test, null, covariates,
  replication unit, multiplicity control, and empirical false-discovery
  endpoint; GeneFunnel scores alone have no false-positive rate.
- [ ] Match random control sets on size, abundance, detection, correlation/
  co-expression, and catalogue overlap/degree where the endpoint requires it.
- [ ] Track compact manifests and reports or archive immutable result bundles;
  ignored local `benchmark/results/` files never become the sole evidence for a
  public claim.
- [ ] Record negative results and narrow the supported domain when simple
  baselines win.

**Go:** the protocol's combined primary-task rule, minimum effects, multiplicity
control, and independent replication pass against sum/mean and relevant methods,
with a clear failure envelope.
**Kill/narrow:** absent held-out benefit, retain GeneFunnel as the faithful thesis
implementation and document the diagnostic/aggregation findings without claims
of general method superiority.

## 11. Conditional horizon

| Direction | Trigger | Boundary / kill rule |
|---|---|---|
| Overlap-aware pathway landscapes | Components + stability show catalogue redundancy obscures interpretation | Separate downstream layer; scores remain computationally independent but are not statistically independent/redundancy-free |
| Empirical null calibration | A precise inferential estimand and replicated negative controls exist | Companion analysis; ship only if simulations and real controls meet declared error rates |
| Paired directional programmes | Validated need for up/down signatures on compatible non-negative assays | Separately named contrast of frozen scores; no signed input to the core kernel |
| Partitioned/federated manifests | A real multi-site workflow adopts exact catalogue/preprocessing contracts | Provenance + concatenation only; no privacy or automatic cross-site comparability guarantee |
| GPU or alternate native engine | Profiles show the kernel, not I/O/catalogue/output, dominates target workloads | Stop if transfer/setup misses the pre-specified end-to-end gain or portability regresses |
| Interactive UI/server | Repeated usability evidence shows the diagnostics need guided exploration | Companion surface after APIs stabilize; no server or Shiny dependency in the scorer |

Foundation models, GNNs, automatic multi-omics fusion, learned gene weights, and
bundled catalogue retrieval have no current evidence-backed role. They enter the
portfolio only through a concrete scientific hypothesis and the same held-out,
maintenance, and failure gates.

## 12. Cross-cutting acceptance gates

Every merged milestone must satisfy the applicable rows; one success never
compensates for a failed correctness gate.

| Gate | Required evidence |
|---|---|
| Core compatibility | Existing API/spec unchanged for existing calls; complete oracle/invariant/adversarial suite; controlled assertions + same-run paired digests |
| Numerical | Independent reference; extreme finite, missing, all-zero, dense/sparse, serial/parallel, block-partition tests |
| Scientific | Estimand + failure domain stated; simulation and negative controls; held-out/replicate evidence before biological claims |
| Performance | Fresh fixed baseline; paired/interleaved warm/cold repeats; estimator + uncertainty/equivalence rule; elapsed time, allocation, process RSS, output size, digest, hardware/software context |
| Dependency | Current official docs + alternatives reviewed; license, maintenance, portability, transitive cost, and optionality recorded |
| Provenance | Formula/protocol/data/catalogue versions; immutable IDs/checksums; generated result policy; AI assistance citations |
| Package/release | Generated interfaces current; clean tarball install/tests; `R CMD check`, `BiocCheck`, sanitizers, and supported-OS CI proportional to change |
| Claims | Docs distinguish theorem, synthetic evidence, empirical evidence, hypothesis, and inference; no universal speed/biology/comparability language |

## 13. Execution order

1. Close or explicitly defer the release lane.
2. A1 theorem/semantics → A2 component prototype.
3. In parallel: F protocol design and C compiled-object spike.
4. B aggregation audit → controlled validation; E reliability → held-out masking.
5. C production compiler/provenance only if its benchmark gate passes.
6. D block engine → HDF5 reference adapter → optional experiment container.
7. Run F held-out programme; promote only the directions that survive.
8. Revisit the conditional horizon from measured bottlenecks and scientific
   demand, not novelty pressure.
