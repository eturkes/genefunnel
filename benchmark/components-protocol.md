<!-- Assisted-by: OpenAI Codex. -->

# Component implementation protocol A2-1.0.0

**Locked:** 2026-07-18, before component native/public implementation.
**Runtime baseline:** `9b60a3eb138e5fd267586624ccd8bf51907577e7`.
That commit has no `R/`, `src/`, `DESCRIPTION`, or `NAMESPACE` difference from
the green candidate boundary `c789a23`. Results generated under
`benchmark/results/` remain untracked; the roadmap records compact evidence.

This protocol fixes A2's numerical decisions, output representation, workloads,
endpoints, uncertainty, and rejection rules. A result cannot move these gates.

## Public result shape

`genefunnel_components(mat, gene_sets, BPPARAM)` will reuse the scorer's input
contract and return this fixed-order named list:

1. `score`, `observed_sum`, `penalty`, `balance`: base double matrices;
2. `effective_size`: base integer matrix;
3. `observed_fraction`: base double matrix;
4. `status`: named list of aligned base character matrices `semantic`,
   `observed_sum`, `penalty`, `balance`, and `conditioning`; and
5. `scaled`: named list `observed_sum`, `penalty`, `balance`, each containing
   aligned `mantissa` double and `exponent` integer matrices.

Every matrix has the retained-set/sample dimensions, order, and dimnames of
`genefunnel()`. `score` is its authoritative matrix, never reconstructed.
Status component matrices use `ordinary`, `scaled`, or `unavailable`.
`semantic` uses `scoreable`, `too_few_observed`, or `zero_total`;
`conditioning` uses `safe`, `ill_conditioned`, or `not_applicable`.

An ordinary component is finite. A scaled component has ordinary `NA` plus the
canonical pair `(m, e)`, value `m * 2^e`, with `0.5 <= m < 1`; all other pair
cells are `NA`. Exact defined zero is ordinary zero. Undefined or uncertifiable
components are unavailable with ordinary value and pair both `NA`. Finite valid
input should make `observed_sum` ordinary or scaled; `unavailable` remains a
fail-closed state for diagnostics whose certification cannot be completed.

## Independent numerical oracle

`tests/testthat/helper-scaled-reference.R` is test-only arithmetic: a
double-double (~106-bit) significand plus an explicit binary exponent. It uses
neither the native scorer nor platform `long double`. It calculates penalty
from absolute deviations, score from the non-negative below-mean identity, and
balance by scaled division. This is a bounded reference for the committed
fixtures, not a general arbitrary-precision package or a public dependency.
Current [Rmpfr](https://cran.r-project.org/package=Rmpfr) was reviewed: it
offers arbitrary precision but adds compiled GMP/MPFR system requirements and
another R dependency. That surface is disproportionate to this finite,
test-only oracle.

Let binary64 unit roundoff `u = .Machine$double.eps / 2`, effective size `n`,
and subtraction condition number `kappa = (T + Q) / F`. Define the conservative
implementation allowance

```
g(n) = 64 * (n + 1) * u.
```

A positive-total scoreable cell is `safe` only when all four diagnostics are
ordinary, `F` is normal and positive, `g(n) < 1`, and
`g(n) * kappa <= 2^-20`. It is otherwise `ill_conditioned`. Fewer-than-two and
zero-total cells are `not_applicable`. Within the safe region both direct
binary64 identities must satisfy

```
abs((T - Q) - F) <= 8 * g(n) * max(T, Q, F)
abs((T * B) - F) <= 8 * g(n) * max(T, Q, F).
```

The constant is a deliberately conservative implementation budget, not a
novel floating-point theorem. Native review must account for every accumulator
and stay within it. Outside the safe region, tests require correct status,
availability, normalized scaled pairs, and oracle agreement; direct binary64
subtraction is explicitly not evidence.

Committed adversarial cells include:

| Cell | Required fact |
|---|---|
| `c(.Machine$double.xmax, .Machine$double.xmax / 2)` | score representable; sum scaled; penalty/balance ordinary |
| `c(.Machine$double.xmax, .Machine$double.xmax, 0, 0)` | score representable; sum and penalty scaled |
| `c(2^1023, 2^-1022, 0)` | sum and penalty round equal; positive score ordinary; balance pair approximately `(0.75, -2044)` |
| `c(2^-1074, 2^-1074)` | subnormal sum/score retained; balance one |
| signed zeros, `NA`/`NaN`, one observed member, all zero | semantic table in `COMPONENTS_SPEC.md` retained exactly |

Member permutations and dense/sparse/serial/SOCK representations cover the
same cells. A scaled pair agrees when the normalized exponent matches the
oracle and relative significand error is at most `8 * g(n)`; values straddling
a normalization boundary compare by scaled relative error instead.

## Default-path performance gate

The primary risk is `genefunnel()` with diagnostics unrequested. The comparison
uses freshly installed baseline/candidate source snapshots at their exact Git
SHAs, in separate libraries but sharing one
dependency library. `prepare-components.R` creates each snapshot with
`git archive`, installs it, fingerprints the source archive and complete
installed package tree, and writes a marker that gate mode recomputes and
validates. Gate mode also requires a clean runner repository at the full
candidate SHA; every worker records and pair-checks R plus dependency versions
and paths. The runner pins `components-protocol.tsv` at MD5
`974e3aeab67a8a5c21a01cd0f154a4ed`; any byte change requires a versioned
protocol revision. That table fixes two overlap patterns for each dense/sparse
storage path. Fixtures reuse the existing deterministic
constructors. Scoring uses explicit `SerialParam()` to isolate the common
manager/native path; parallel equivalence remains a correctness gate.

Gate evidence is Linux-only and requires `/proc/loadavg` plus compatible GNU
`time`. At start and before every adjacent pair, one-minute load divided by
detected logical CPUs must be at most `0.25`; otherwise the run aborts and its
partial observations are invalid. `load.tsv` retains every check. This
quiescence threshold was fixed after a pre-implementation identical-binary
calibration on an externally saturated host falsified the unguarded design.

- 30 paired repeats per workload. Each tracked order string contains exactly
  15 baseline-first and 15 candidate-first pairs, randomized before native
  implementation with base-R seed `18072026` and independently across
  workloads. Workloads rotate within each repeat.
- Every observation uses a fresh `R --vanilla` process. Fixture construction,
  package loading, and `gc()` occur before a recorded cold scoring call. Five
  subsequent complete calls form the warm timing batch; elapsed time per call
  is the batch time divided by five.
- Primary estimator per workload = geometric mean candidate/baseline warm
  elapsed ratio. Uncertainty = one-sided 95% Student interval on paired log
  ratios. The upper confidence bound must be at most `1.05` for every workload.
  The separately reported cold-call interval is context, not an equivalence
  gate: pre-implementation self-calibration showed first-call host variance is
  too large for a truthful 5% decision.
- A separate `Rprofmem()` call must match the cold/warm output. Compatible GNU
  `time` passively captures whole-process maximum RSS; no concurrent sampler
  may contend with the timer. Unsupported metrics record `NA`, never zero.
- Manager allocation gate per workload: candidate median <= baseline median +
  `max(65536 bytes, 1% of baseline)`. This and code review enforce that the
  ordinary path allocates no diagnostic result matrices.
- Process maximum RSS gate: candidate median <= baseline median +
  `max(4096 KiB, 10% of baseline)`. It includes the identical fixture/package
  startup and all scoring passes; its role is to catch material default-path
  growth, not attribute native bytes precisely.
- Every paired output digest, dimensions, names, missingness, and byte size must
  agree exactly. Any unavailable metric, failed identity, threshold failure,
  load excursion, or incomplete repeat rejects the candidate; repeats are not
  extended after inspection.

The non-gating smoke mode downscales each matrix and uses four repeats solely to
exercise orchestration. Only the fixed 30-repeat gate mode can support A2.

## Decision

All numerical/correctness gates plus all four runtime confidence bounds and
allocation/RSS gates must pass. Otherwise optimize without changing semantics,
or keep the theorem/specification internal as PLAN's fallback requires.
