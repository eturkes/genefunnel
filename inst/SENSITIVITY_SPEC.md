<!-- Assisted-by: OpenAI Codex. -->

# GeneFunnel observed-member sensitivity specification

This document fixed Workstream E's mathematical diagnostic before an
implementation or empirical result existed. It defines deterministic score
changes caused by deleting one currently observed gene-set member. It does not
define causal contributions, sampling uncertainty, or a correction to the
GeneFunnel score.

## Status and scope

The frozen score and input semantics in `SCIENTIFIC_SPEC.md` remain
authoritative. Protocol E-1.0.0 fixed the internal compact schema, exact
dyadic-rational brute oracle, binary64 availability contract, and controlled
held-out gates before implementation. The package now contains that unexported
prototype; this specification adds no public function. The complete frozen
controlled test failed both co-primary prediction gates, selecting the internal
oracle/no-public-API fallback. Any later public interface requires a newly
scoped protocol, external evidence, and a separate scope decision.

The delete-one construction resembles the operation used by the jackknife, but
gene-set members are not assumed to be independent and identically distributed
sampling units. These diagnostics therefore are not jackknife bias estimates,
variance estimates, standard errors, confidence intervals, or influence
functions. That inferential boundary follows the sampling assumptions of the
classical jackknife literature, including
[Efron and Stein (1981)](https://doi.org/10.1214/aos/1176345462).

## Canonical observed members

For each declared gene set:

1. remove duplicate identifiers while preserving first occurrence;
2. retain exact, case-sensitive matches to matrix rows in that order; and
3. for each sample, remove members whose value is `NA` or `NaN`, preserving the
   remaining order.

The resulting ordered identifiers are the cell's canonical observed members,
with finite non-negative values \(x=(x_1,\ldots,x_n)\). Explicit and implicit
zeros remain observed. Globally unmatched members and sample-specific missing
members receive no deletion delta. Duplicate declarations never create
multiple deltas.

## Exact deletion delta

Let \(F_m\) denote the frozen GeneFunnel equation evaluated on \(m\) observed
values. For a cell with \(n\ge3\), define the full-minus-deleted delta for
canonical member \(i\) as

\[
\Delta_i(x)=F_n(x)-F_{n-1}(x_{-i}),
\]

where \(x_{-i}\) removes coordinate \(i\) and changes the score's effective-size
coefficient from \(n/[2(n-1)]\) to \((n-1)/[2(n-2)]\). At least three observed
members are required so every deleted vector remains scoreable.

- \(\Delta_i>0\): deleting member \(i\) lowers the score.
- \(\Delta_i<0\): deleting member \(i\) raises the score.
- \(\Delta_i=0\): the exact score is unchanged by deletion.

The sign describes a deterministic perturbation of this cell only. It does not
say that a gene activates, represses, explains, or causes the score.

## Compact mathematical summaries

For \(n\ge3\), let

\[
i^*=\min\operatorname*{arg\,max}_{i=1,\ldots,n}|\Delta_i|,
\]

where the minimum is canonical member order. The compact mathematical summary
contains these quantities:

- effective size \(n\);
- member identifier \(i^*\);
- largest absolute delta \(|\Delta_{i^*}|\);
- its signed delta \(\Delta_{i^*}\);
- its signed observed-sum normalization \(\Delta_{i^*}/T\), only when
  \(T=\sum_i x_i>0\); and
- \(\operatorname{median}_i |\Delta_i|\), using the ordinary midpoint for an
  even number of members.

The normalization is a scale-free perturbation size, not a fraction of the
score or a compositional allocation. It is `NA` when \(T=0\). Exact ties,
including an all-zero vector, select the first canonical observed member. For
\(n<3\), effective size remains factual and every deletion-summary quantity is
undefined.

Full member-level deltas remain a separate prospective output requiring an
explicit sink; the compact interface must not allocate a hidden
member-by-set-by-sample array.

## Exact properties

For finite non-negative \(x\) with \(n\ge3\):

- **Bounds.** Because \(0\le F_n(x)\le T\) and
  \(0\le F_{n-1}(x_{-i})\le T-x_i\),
  \(-(T-x_i)\le\Delta_i\le T\). Hence
  \(-1\le\Delta_i/T\le1\) when \(T>0\).
- **Positive homogeneity.** For \(a\ge0\),
  \(\Delta_i(ax)=a\Delta_i(x)\). The normalized delta is invariant for
  \(a>0\); at \(a=0\) it is undefined.
- **Common shifts.** For any \(c\) preserving non-negativity,
  \(\Delta_i(x+c\mathbf1)=\Delta_i(x)+c\). Sensitivity therefore retains the
  score's input units and meaningful-zero requirement.
- **Permutation equivariance.** Reordering members reorders their attached
  deltas. Largest absolute and median absolute values are invariant; the
  selected identifier can change only within an exact tie because canonical
  declaration order is the tie-breaker.
- **Missingness distinction.** A missing member is omitted before \(n\) and the
  deltas are defined. Replacing it with zero changes the estimand. Adding a
  globally unmatched declaration changes neither values nor deltas.

The deltas are not additive attributions. In general,
\(\sum_i\Delta_i\ne F_n(x)\), and the changed effective-size coefficient means
\(\Delta_i\) is not the member's measured value or a separable term of the
score.

## Unknown-member identification

An unknown member is a canonically declared identifier whose value is absent
for the cell, either because its feature row is globally unmatched or its cell
is `NA`/`NaN`. A measured zero is observed and is not an unknown. Undeclared
features remain outside the gene-set estimand. The score and deletion diagnostic
omit unknown members; they neither impute them nor treat them as zero.

Non-negativity alone assigns each unknown value the set \([0,\infty)\), so no
unknown coordinate has a finite upper bound. This does not always make the
complete-data score unbounded. Let finite observed values \(x\) have length
\(r\), sum \(T\), and let \(q\) unknown members all take value \(L\), with
\(m=r+q\ge2\). Once \(L\) is large enough,

\[
F_m(x,L\mathbf1_q)=
\frac{q(q-1)}{m-1}L+\frac{r+2q-1}{m-1}T.
\]

For \(q=1\) and \(r\ge1\), the first coefficient is zero: the score plateaus at
\((r+1)T/r\), even though the missing value remains unidentified. For
\(q\ge2\), the coefficient is positive and the complete-data score is
unbounded. A finite score plateau therefore must not be presented as a bound
on the missing member, evidence for zero imputation, or reliability.

Finite caller limits produce deterministic bounds without adding a sampling
claim. For componentwise \(0\le \ell_j\le z_j\le u_j<\infty\), define complete
corners \(y^L=(x,\ell)\) and \(y^U=(x,u)\). The score is coordinatewise
nondecreasing: increasing one coordinate by \(h\) shifts its centered residual
by \(h(1-1/m)\), every other residual by \(-h/m\), and can increase total
absolute deviation by at most \(2h(m-1)/m\). The score increase is therefore
non-negative. Consequently,

\[
F_m(y^L)\le F_m(x,z)\le F_m(y^U)
\]

is the sharp score identification interval over the closed box.

For \(m\ge3\), each potential complete-data deletion delta has the finite outer
enclosure

\[
F_m(y^L)-F_{m-1}(y^U_{-i})
\le \Delta_i(x,z) \le
F_m(y^U)-F_{m-1}(y^L_{-i}).
\]

This enclosure is valid but need not be sharp because its two score extrema can
require different unknown vectors. Member ranking, the selected largest member,
and compact summaries remain unidentified unless they agree throughout the
feasible box. Any future implementation must key limits to canonical member
identifiers, retain their units/preprocessing and source, distinguish global
absence from sample missingness, and expose assumption-dependent bounds rather
than a corrected score. The failed E-1.0.0 reliability gate does not justify
such a public interface.

Finite limits are constraints, not a probability distribution. Probability,
confidence, or credible intervals require an explicit joint value and
missingness model, dependence assumptions, a fitting/validation protocol, and
an inferential target. That belongs to a separately pre-specified project.

## Canonical cases

| Canonical observed values | Deletion deltas | Selected member | Signed / normalized delta | Median absolute delta |
|---|---|---:|---:|---:|
| `A=1, B=2, C=7` | `A=0.5, B=2.5, C=2.5` | `B` | `2.5 / 0.25` | `2.5` |
| `A=0, B=1, C=1` | `A=-1, B=1, C=1` | `A` | `-1 / -0.5` | `1` |
| `A=4, B=4, C=4` | `A=4, B=4, C=4` | `A` | `4 / 1/3` | `4` |
| `A=0, B=0, C=0` | all `0` | `A` | `0 / NA` | `0` |
| `A=1, B=2` | undefined | `NA` | `NA / NA` | `NA` |

For the first row, the full score is 4.5 while the deltas sum to 5.5. For the
second row, deleting the observed zero raises the score from 1 to 2; its
negative delta wins an absolute-value tie by canonical order. These cases
prevent contribution-like interpretations and sign erasure.

## Numerical and interpretation boundary

The exact definition is over real arithmetic. Protocol E-1.0.0 represents each
binary64 input as an exact common-power-of-two integer and compares signed
deletion numerators exactly. Sign, order, and ties therefore do not come from
subtracting rounded scores. A nonzero summary that cannot be represented as an
ordinary finite/nonzero double is explicitly unavailable; it never becomes
infinity or false zero. The prototype has no scaled sidecar. The ordinary
`genefunnel()` path allocates no sensitivity output, and its score remains
authoritative.

The retained brute oracle recomputes every deleted score numerator. After the
frozen profile admitted optimization, the active internal path sorts exact
member magnitudes once and forms exact prefix sums. For deletion `i`, it sets
`T_-i = T - z_i`, binary-searches the largest sorted boundary `q` satisfying
`(n - 1) z_(q) < T_-i`, removes `z_i` from the prefix when it lies below that
boundary, and evaluates the same integer `N_-i`. This changes complexity, not
the rational estimand. The brute function remains executable. A clean isolated
candidate reproduces the frozen brute workload digest exactly; fixed and
randomized tests also require exact delta-object identity.

Sensitivity depends on input units, preprocessing, effective size, retained
member identities, gene-specific abundance and dynamic range, measurement
noise, and the meaningful-zero assumption. It is not causal attribution,
feature importance, biological necessity, assay-independent reliability, or
inferential uncertainty. Whether it predicts held-out feature-loss error or
technical-repeat instability was an empirical question governed by the frozen
protocol, not implied by this algebra. E-1.0.0 answered it negatively in its
controlled design: neither augmented model approached the pre-specified
incremental prediction threshold. Biological-replicate variation remains a
different, untested estimand that may contain real signal rather than error.

## Prior art and novelty boundary

Audit date: 2026-07-18. This targeted, non-systematic audit covers primary
method papers relevant to deletion, gene-set member diagnostics, missing-gene
stability, and subset refinement. It is not an exhaustive review or patent
search.

| Area | Prior art | Boundary here |
|---|---|---|
| Delete-one resampling | [Efron and Stein (1981)](https://doi.org/10.1214/aos/1176345462) | Classical jackknife variance work treats observations as sampling units. Gene-set members are not assigned that role here; no bias, variance, confidence, or influence-function inference follows. |
| Leading-edge members | [Subramanian et al. (2005)](https://doi.org/10.1073/pnas.0506580102) | GSEA's leading edge is a phenotype/ranked-list subset driving an enrichment peak. It is not a deterministic deletion delta for one GeneFunnel cell. |
| Per-sample score dispersion | [Foroutan et al. (2018)](https://doi.org/10.1186/s12859-018-2435-4) | `singscore` reports rank-score dispersion. That is useful API precedent but not score change after deleting a member. |
| Missing-gene benchmarks | [Toro-Dominguez et al. (2025)](https://doi.org/10.1093/bib/bbaf684) | Their multi-method benchmark shows that missing genes can destabilize some single-sample scores. It does not supply a per-cell diagnostic that predicts held-out GeneFunnel error. |
| Gene-set stability/refinement | [Abou Choucha and Pasquier (2026)](https://doi.org/10.1038/s41598-026-48119-9) | MEAST uses SVD scores plus a genetic algorithm to optimize active subsets across cell groups. GeneFunnel sensitivity neither learns/refines a set nor calls selected members stable or high-contributing. |

Delete-one perturbation, member ranking, dispersion, missing-gene robustness,
and gene-set refinement are established ideas. Any claim is restricted to the
exact GeneFunnel full-minus-deleted quantity, its faithful compact
implementation, and predictive value that survives pre-specified held-out
gates. A largest-sensitivity member is not a novel biomarker or biological core.

## Executable evidence

`tests/testthat/helper-sensitivity-reference.R` implements direct brute-force
deletion from the normative equation without package code.
`test-sensitivity-theorem.R` locks the canonical cases, support rules, sign,
tie-break, bounds, homogeneity, common-shift behavior, permutation behavior,
non-additivity, coordinate monotonicity, unknown-member non-identification, and
finite-limit score bounds/delta enclosures before a public implementation
exists.
`test-sensitivity-exact.R` checks the limb arithmetic against independent
integer identities, round-trips binary64 across its exponent domain, proves
exact ordering under subtract-rounded tie drift, and checks score-numerator and
wide-exponent conversion. `test-sensitivity-api.R` locks the compact schema,
statuses, direct randomized recomputation, canonical support, extreme values,
storage representations, and serial/SOCK identity.
`benchmark/sensitivity-protocol.md` and its machine registry fix the internal
schema, exact-arithmetic boundary, controlled masks/repeats, folds, models,
uncertainty, and rejection rules; they were committed before a package
diagnostic was calculated. Profile supplement `E-P-1.0.0` transparently closes
the parent's previously implicit fixture/pass mechanics after implementation
but before a fixed profile call; it changes no empirical gate. Its tracked
result makes exact-oracle-preserving optimization eligible and no reliability
or public-interface claim.
Controlled supplement `E-C-1.0.0` then pins the previously implicit R RNG
calls, predictor encoding/scaling, bootstrap draw order, isolated execution,
checkpoint, and evidence mechanics before any controlled scenario is built. It
changes no scientific design, target, threshold, or promotion rule.
The controlled observation implementation uses the installed authoritative
scorer and exact cell diagnostic, resets every scenario/mask seed, retains all
registered factors, and validates target/status/encoding identities before any
held-out model can consume a row.
The model/runner layer records training-only scaling, zero-SD drops, QR aliases,
fixed predictions, scenario-cluster bootstrap draws, all descriptive strata,
isolated source/install fingerprints, resumable checkpoints, and artifact
hashes. Its smoke uses planted targets and has no empirical meaning. The
complete tracked [`controlled result`](../benchmark/sensitivity-controlled-result.md)
retains all 345,600 feature-loss and 5,760 repeat rows. Feature-loss median
fold reduction/bootstrap lower bound was 0.00110097/0.000407703 and controlled-
repeat reduction/lower bound was 0.0193130/0.0112902; every value failed its
0.10/0.05 requirement. The deterministic thinning curves remain descriptive,
study-composition-dependent artificial deletion evidence and cannot rescue a
failed endpoint. The exact diagnostic stays an internal adversarial-test oracle.
