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
prototype; this specification adds no public function. A later public interface
still requires external evidence and a separate scope decision.

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

Sensitivity depends on input units, preprocessing, effective size, retained
member identities, gene-specific abundance and dynamic range, measurement
noise, and the meaningful-zero assumption. It is not causal attribution,
feature importance, biological necessity, assay-independent reliability, or
inferential uncertainty. Whether it predicts held-out feature-loss error or
technical-repeat instability is an empirical question governed by a frozen
protocol, not implied by this algebra.

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
and non-additivity before an optimized or public implementation exists.
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
but before a fixed profile call; it changes no empirical gate.
