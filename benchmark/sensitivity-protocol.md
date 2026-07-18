<!-- Assisted-by: OpenAI Codex. -->

# Observed-member sensitivity protocol E-1.0.0

**Frozen:** 2026-07-18, after the mathematical contract and before any package
sensitivity value, performance profile, controlled target, fitted model, or
endpoint exists. The adjacent
[`machine registry`](sensitivity-protocol.tsv) is authoritative. Changing a
field, seed, split, predictor, threshold, or decision rule requires a new
committed protocol version before replacement results. The validator pins the
registry bytes at MD5 `09d9429ae18f64ba00b3bd84955ad71d`.

This first stage asks three ordered questions:

1. Can a compact internal prototype reproduce exact full-minus-deleted
   sensitivities without changing the scorer?
2. Do the summaries predict error after independently seeded feature loss
   beyond score/support/magnitude facts?
3. Do they predict controlled measurement-repeat score instability beyond the
   same facts?

The algebra can remain correct when both prediction gates fail. E-1.0.0 has no
external biological data and cannot support an export, causal claim,
assay-independent reliability claim, or technical-replicate claim about a real
experiment. Passing only permits a separately frozen external stage.

## Internal prototype

The first implementation is unexported:

```r
.gene_set_sensitivity(mat, gene_sets, BPPARAM = BiocParallel::bpparam())
```

It reuses the scorer's complete matrix/gene-set/backend validation, stable
member deduplication and matching, aggregate omission warning, retained set
universe (`matched_size >= 2`), set/sample order, and names. A set with two
global matches remains aligned even though no cell in it can have three
observed members. Dense and sparse inputs keep their storage class through
bounded column iteration. The ordinary `genefunnel()` path calls no sensitivity
code and allocates no sensitivity output.

The fixed-order result is an unclassed list of aligned base matrices:

1. `largest_member`: character;
2. `largest_absolute_delta`: double;
3. `largest_delta`: double;
4. `largest_delta_over_sum`: double;
5. `median_absolute_delta`: double;
6. `effective_size`: integer; and
7. `status`: character matrices `semantic`, `delta`, and `normalized`.

`semantic` is `defined` or `too_few_observed`. `delta` is `ordinary`,
`unavailable`, or `not_applicable`. `normalized` is `ordinary`, `zero_total`,
`unavailable`, or `not_applicable`.

- Fewer than three observed members: only `effective_size` is populated;
  summary values are `NA`, both availability fields are `not_applicable`.
- At least three all-zero members: deltas are exact positive zero, the first
  canonical member wins, `delta = ordinary`, and
  `normalized = zero_total` with an `NA` normalized value.
- Numerically unavailable raw summaries: member and all three raw summary
  values are `NA`; normalized value is also `NA`/`unavailable`.
- Available raw summaries with an unrepresentable nonzero normalization keep
  their raw fields and use `NA`/`unavailable` only for normalization.

The compact call never constructs or returns a member-by-set-by-sample array.
Member-level output requires a future explicit sink contract.

## Exact brute-force arithmetic

The prototype treats each finite binary64 input as its exact dyadic value. For
one observed vector, choose a common power-of-two exponent `e` and non-negative
arbitrary-length integers `z_i` such that `x_i = z_i * 2^e`. Signed zero maps to
integer zero. The implementation is dependency-free and uses exact limb
addition, comparison, subtraction, and multiplication by a machine-sized
integer. An all-zero vector uses `e = 0` and zero integers.

For any retained vector of size `m`, define exact integers

\[
T_m=\sum_jz_j,\qquad
L_m=\sum_{j:mz_j<T_m}z_j,\qquad
\ell_m=|\{j:mz_j<T_m\}|,
\]

and

\[
N_m=(m-1-\ell_m)T_m+mL_m.
\]

The non-negative below-mean identity gives

\[
F_m(x)=2^e\frac{N_m}{m-1}.
\]

For the full size `n` and deletion `i`, every delta has the same positive
denominator and exponent:

\[
\Delta_i
=2^e\frac{D_i}{(n-1)(n-2)},\qquad
D_i=N_n(n-2)-N_{-i}(n-1).
\]

Thus exact signed integers `D_i` determine sign, absolute ordering, exact ties,
the canonical winning member, and the central order statistics. The first
prototype recomputes every `N_{-i}` independently; this is the required brute
oracle, not an optimized prefix algorithm.

Raw and normalized outputs are rounded from the exact rational values. Each
ordinary nonzero value must be finite/nonzero in binary64 and agree with an
independent scaled reference within relative error
`8 * .Machine$double.eps * max(1, effective_size)`. Exact zero returns positive
zero. A mathematical nonzero that overflows or underflows becomes `NA` with
`unavailable` status, never infinity or false zero. The normalized selected
delta divides the exact `D_i` by the exact observed total and is handled
separately. This prototype supplies no scaled sidecar.

Correctness requires exact agreement with direct normative recomputation for
ordinary fixtures, the committed sign/tie cases, wide binary exponents,
missingness, signed zero, dense/sparse storage, serial/SOCK execution, and
member/set/sample permutations. It also checks that full-cell scores agree
with the authoritative scorer within its existing tolerance. A malformed
internal schema or arithmetic invariant aborts; a representational limit is a
cell status.

## Profile before optimization

The exact R brute path is profiled before any sorted-prefix or native
acceleration. The fixed serial workload has 4,096 features, 12 samples, 48
sets of 128 members, deterministic dense values/members, and three fresh
repeats. Record elapsed time, `Rprof()` samples, `Rprofmem()` allocation, output
digest, and environment. Performance is descriptive.

E-1.0.0 fixed the dimensions and seeds but left the deterministic value/member
constructors and profiling-pass boundaries implicit. This is a genuine
execution-specification omission. Byte-pinned supplement
[`E-P-1.0.0`](sensitivity-profile-protocol.tsv) closes only those degrees of
freedom after the brute implementation and before any fixed profile call. It
uses column-major unit-rate exponential values, independent within-set sampling
without replacement, three GC-separated timed calls, then separate `Rprof()`
and `Rprofmem()` calls whose five outputs must have one digest. An exact sample
is a stack containing a `.gf_` frame. The supplement does not retroactively
claim pre-implementation timing specification and changes no correctness,
controlled-design, endpoint, or promotion gate.

Optimization work becomes eligible only if the median call exceeds 60 seconds
or exact deletion arithmetic consumes more than half the sampled time. Any
optimized candidate must first match the exact brute result/status on every
fixed and randomized cell. No speed result can relax correctness or promote
the API.

The tracked [`E-P-1.0.0 result`](sensitivity-profile-result.md) records median
elapsed 213.585 seconds and exact-stack share 0.999402883221019. Both branches
cross their frozen strict thresholds; optimization research is eligible. The
result sets `performance_claim = FALSE` and changes no reliability gate.

The tracked [`sorted-prefix adoption check`](sensitivity-optimization-result.md)
installs clean candidate `e5c013f` in isolation and reproduces the fixed brute
output MD5 exactly. Exact-object tests remain authoritative; the exploratory
single-call timing is not a performance claim.

The parent controlled registry fixes scientific dimensions, equations, seeds,
predictors, and gates, but terms such as Bernoulli draw, global schema, and
fixed-seed bootstrap do not uniquely determine an R execution. Prospective
supplement [`E-C-1.0.0`](sensitivity-controlled-protocol.tsv) byte-pins only
those implementation choices before any controlled scenario is constructed.
It changes no parent design, target, endpoint, threshold, claim boundary, or
promotion rule.

## Controlled design

The full Cartesian design uses registry order, with the first field varying
fastest:

- member count: 8, 32, 128;
- profile: monotone, U-shaped, wave, dominant;
- dynamic range: 4, 64;
- log-normal member perturbation SD: 0.1, 0.5;
- expected set total: 100, 1,000, 10,000 counts;
- negative-binomial dispersion: 0, 0.1;
- independent post-count dropout: 0, 0.2; and
- profile replicate/fold: 1 through 10.

This gives 5,760 latent scenarios and exactly 576 scenarios per held-out fold.
Base R uses `Mersenne-Twister`, `Inversion`, and `Rejection`; every scenario and
measurement resets `base + scenario_id - 1`, so worker count cannot change a
record.

For member coordinate `t_i = (i - 1)/(n - 1)` and dynamic range `r`, the four
unscaled profiles are

\[
q_i=\begin{cases}
r^{t_i-1/2} & \text{monotone},\\
r^{2|t_i-1/2|-1/2} & \text{U-shaped},\\
r^{\sin(2\pi t_i)/2} & \text{wave},\\
\sqrt r\;\mathbf1(i=1)+r^{-1/2}\;\mathbf1(i>1)
& \text{dominant}.
\end{cases}
\]

The latent seed draws canonical-order `z_i ~ N(0, 1)` and replaces
`q_i` by `q_i * exp(sigma * z_i - sigma^2 / 2)`, then normalizes to shares
`p_i`. A fixed assay-capture pattern `c_i = 2^sin(4*pi*t_i)` produces technical
means

\[
\mu_i=D\frac{p_ic_i}{\sum_jp_jc_j}.
\]

For dispersion `phi = 0`, counts are independent Poisson draws. Otherwise
they are independent negative-binomial draws with mean `mu_i` and
`size = 1 / phi`, so variance is `mu_i + phi * mu_i^2`. Each measurement then
draws one Bernoulli dropout indicator per member in canonical order and sets
selected counts to observed zero. Replicates A and B use independent fixed
seeds but the same latent means. They model controlled measurement repeats,
not biological replicates. Post-count dropout is a stress operator, not
evidence that technical dropout is a scientifically meaningful observed zero
in real single-cell data.

The model-implied positive-detection probability used only to stratify masks is

\[
\pi_i=(1-d)\{1-P(Y_i=0)\},
\]

where the count-zero probability is `exp(-mu_i)` for Poisson and
`(1 + phi * mu_i)^(-1 / phi)` otherwise. It is simulation-model dependent and
must not be presented as an assay-independent gene property.

## Independently seeded feature loss

Every scenario crosses removed fractions 0.125, 0.25, and 0.5; mask mechanisms
`uniform`, `low_abundance`, `high_abundance`, `low_detection`, and
`high_detection`; two mask repeats; and two absence encodings. Removed count is
`max(1, min(n - 3, floor(fraction * n + 0.5)))`.

Masks use base `sample.int(..., replace = FALSE)` after rescaling positive
weights by their maximum:

- uniform: `1`;
- low abundance: `1 / p_i`;
- high abundance: `p_i`;
- low detection: `1 - pi_i + 1/16`; and
- high detection: `pi_i + 1/16`.

Selected indices are sorted back to canonical order. The same mask is encoded
twice:

- `global_absence`: selected feature rows are absent, so declared coverage
  falls while observed fraction is one; and
- `sample_missing`: feature rows remain but selected A cells are `NA`, so
  declared coverage is one while observed fraction falls.

The observed values, scores, and sensitivities must be identical across the
two encodings; their coverage facts must remain distinct. These are artificial
sampled deletions. They are not evidence that a real unmeasured member was zero
or randomly missing. The grid contains 345,600 fixed feature-loss rows.

The executable observation layer resets every registered seed locally,
generates the parent profiles/counts/masks in canonical order, batches the
authoritative installed scorer once per scenario, and calls the installed exact
cell diagnostic only on A/full or A/partial inputs. It computes each unique
partial observation once and duplicates it across the two absence encodings;
only declared coverage and observed fraction differ. Smoke checks also score
and diagnose explicit row absence versus `NA` omission independently.

For full measured replicate A score `F_A`, partial score `F_A^p`, and full
observed sum `T_A`, the feature-loss target is

\[
E_{loss}=\begin{cases}
0 & T_A=0,\\
|F_A^p-F_A|/T_A & T_A>0.
\end{cases}
\]

The target uses the same measurement before the independent mask, not the
latent expectation. Sensitivities and baseline facts use only the partial
input.

For full independent measured replicates A/B, the controlled-repeat target is

\[
E_{repeat}=\begin{cases}
0 & T_A+T_B=0,\\
2|F_A-F_B|/(T_A+T_B) & T_A+T_B>0.
\end{cases}
\]

Sensitivities and baseline facts use A only; B is outcome-only. There are 5,760
fixed repeat rows.

## Held-out models and endpoints

Profile replicate number is the ten-fold split. All masks, both absence
encodings, and A/B measurements for a latent scenario stay in its fold. Each
fold model trains on the other nine profile replicates and predicts the held-out
replicate. No row from a held-out latent profile enters training.

Separate feature-loss and repeat models use these baseline predictors:

- `log2(declared_size)`;
- declared coverage;
- observed fraction;
- `log1p(observed_sum)`;
- effective size;
- score divided by observed sum, or zero when the sum is zero; and
- an all-zero indicator.

For positive totals, score divided by sum is the existing balance component.
The baseline therefore prevents sensitivity from winning merely by proxying
the already-shipped decomposition.

The feature-loss baseline additionally includes absence encoding. The augmented
model adds absolute selected delta divided by observed sum, signed selected
delta divided by observed sum, and median absolute delta divided by observed
sum. Each is zero when the observed sum is zero; the baseline already carries
the exact zero indicator. No latent factor, mask mechanism/fraction, expected
abundance, detection propensity, depth, dispersion, or dropout enters either
model.

A global column schema is fixed before folds. Non-binary continuous predictors
are centered/scaled using training means/SDs and applied unchanged to held-out
rows. A zero-SD training column is dropped from both models for that fold.
Base R `lm.fit` supplies pivoted QR least squares; aliased coefficients are
zero. Predictions are not clipped. Missing/non-finite targets, predictors, or
predictions fail the fixed grid. Every controlled cell must have ordinary raw
deltas and normalized status `ordinary` or `zero_total`; any other status also
fails rather than disappearing.

For each fold and target, calculate held-out RMSE and

\[
R=(\operatorname{RMSE}_{base}-\operatorname{RMSE}_{aug})/
\operatorname{RMSE}_{base}.
\]

Zero baseline RMSE is undefined and fails. The point endpoint is the median of
ten fold reductions. Uncertainty uses 2,000 fixed-prediction cluster-bootstrap
replicates: resample latent scenarios with replacement inside each fold,
retain every associated mask row, recompute fold RMSE reductions, then take
their median. Models are not refit. The one-sided 95% lower bound is the type-8
0.05 quantile. Feature loss and controlled repeats use separate fixed seeds.

Both targets are co-primary intersection gates. Each requires median RMSE
reduction at least 0.10 and bootstrap lower bound at least 0.05. Requiring both
is an intersection-union decision; one target cannot compensate for the other.
The 10% point threshold denotes a material prediction improvement, while the
5% lower bound excludes a trivial or unstable gain; both were fixed without a
package sensitivity calculation or endpoint calibration.
Factor-, fraction-, mechanism-, and encoding-stratified results are retained
descriptively and cannot rescue a failed endpoint.

The executable model layer materializes each global numeric schema before
folding; records every training mean, SD, zero-SD drop, QR rank, coefficient,
and alias; and retains every fixed held-out prediction. Bootstrap execution
aggregates squared errors within canonical scenario clusters, resamples those
clusters in the frozen replicate/fold order, and is algebraically identical to
duplicating every associated mask row. A clean-archive runner checkpoints
scenario ranges atomically and treats scientific failure as a completed result.

## Decision

Correctness, exact schema/status, representation invariance, feature-loss
increment, and controlled-repeat increment must all pass. Failure selects the
PLAN fallback: keep targeted adversarial tests and omit a public sensitivity
API. A complete pass permits only a new protocol for external technical-repeat
and held-out scientific evaluation; it does not itself justify export.

## Result

The complete tracked [`E-1.0.0 result`](sensitivity-controlled-result.md) from
clean candidate `5920ea9` retained all 345,600 feature-loss and 5,760
controlled-repeat rows. Feature-loss median fold RMSE reduction/bootstrap lower
bound was 0.00110097/0.000407703; controlled-repeat reduction/lower bound was
0.0193130/0.0112902. All four values fail their frozen 0.10/0.05 requirements.
The 30 tracked thinning-curve rows are descriptive, study-composition-dependent
artificial deletions and cannot rescue the co-primary failure. The fallback is
therefore selected: retain exact internal adversarial tests, omit the public
API, and make no technical- or biological-replicate reliability claim.
