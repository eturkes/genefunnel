<!-- Assisted-by: OpenAI Codex. -->

# Aggregation-gap proof and pre-API contract

**Status:** theorem/design contract ratified 2026-07-18; no public aggregation
function exists. The frozen scorer remains governed by `SCIENTIFIC_SPEC.md`.
This document owns Workstream B's estimand, eligibility, missingness,
terminology, and novelty boundary before a prototype can influence them.

## Estimand

Let active units be $x_k=(x_{k1},\ldots,x_{kn})$, $k=1,\ldots,m$, on one
identical ordered member support, with finite non-negative values and $n\ge2$.
Let finite weights satisfy $w_k>0$ and $\sum_k w_k=1$. Zero-weight units are
excluded before value/support checks and reported as excluded; they cannot alter
the estimand or eligibility.

Define

$$
\bar{x}=\sum_k w_kx_k,
\qquad
F(x)=\mathbf 1^{\mathsf T}x-c\lVert Hx\rVert_1,
\qquad
H=I-\frac{\mathbf1\mathbf1^{\mathsf T}}n,
\qquad
c=\frac{n}{2(n-1)}.
$$

The primary discrepancy is the weighted aggregation gap

$$
J=F(\bar{x})-\sum_k w_kF(x_k).
$$

`Weighted mean` is the only primary aggregation estimand. A future interface may
accept non-negative masses $a_k$ with finite positive total $A$, but it must
report and use $w_k=a_k/A$. It must not silently substitute equal weights.

## Exact theorem and proof

Write $\mu_k=n^{-1}\mathbf1^{\mathsf T}x_k$ and
$d_{ki}=x_{ki}-\mu_k$. Linearity gives
$H\bar{x}=\sum_k w_kHx_k$ and exact cancellation of the total-signal terms:

$$
\begin{aligned}
J
&=c\left[\sum_k w_k\lVert Hx_k\rVert_1
-\left\lVert\sum_k w_kHx_k\right\rVert_1\right]\\
&=c\sum_i\left[\sum_k w_k|d_{ki}|
-\left|\sum_k w_kd_{ki}\right|\right]\ge0.
\end{aligned}
$$

The last inequality is the coordinate-wise weighted triangle inequality, so it
does not rely on native scoring code.

For each member define positive and negative centered mass

$$
P_i=\sum_k w_k\max(d_{ki},0),
\qquad
N_i=\sum_k w_k\max(-d_{ki},0).
$$

Since $\sum_k w_k|d_{ki}|=P_i+N_i$ and
$|\sum_k w_kd_{ki}|=|P_i-N_i|$,

$$
J=2c\sum_i\min(P_i,N_i).
$$

Thus the gap measures one exact arithmetic pattern: cancellation between
oppositely signed, within-unit member deviations from their unit means.

### Equality and strictness

$J=0$ exactly when every coordinate has one weak sign across positive-weight
units: for each $i$, $P_i=0$ or $N_i=0$. The gap is positive exactly when at
least one member is above its unit mean in one active unit and below its unit
mean in another.

Profiles $x_k=a_kq$, $a_k\ge0$, are sufficient for equality because all
$Hx_k=a_kHq$ share coordinate signs. Proportionality is not necessary. For
example, `c(4, 3, 1, 0)` and `c(6, 3, 2, 1)` are not proportional, but their
centered deviations have no opposing coordinate and their gap is zero.

The equal mixture of `c(4, 0)` and `c(0, 4)` is the maximum-discrepancy example:
both unit scores are zero, the mean is `c(2, 2)` with score four, and $J=4$.

## Bounds and normalization

For valid non-negative vectors, $F(x_k)\ge0$. Hence

$$
0\le J=F(\bar{x})-\sum_kw_kF(x_k)\le F(\bar{x}).
$$

The provisional normalized gap is

$$
R=\frac{J}{F(\bar{x})}
=1-\frac{\sum_kw_kF(x_k)}{F(\bar{x})}
$$

when $F(\bar{x})>0$, giving $R\in[0,1]$. When $F(\bar{x})=0$, both $J$ and the
weighted unit score are zero; `R` is undefined and must be `NA`, not zero.

`R` answers a narrow intended question: what fraction of the aggregate score is
the aggregate-versus-weighted-unit discrepancy? It is scale-free under one
common positive multiplier, but unstable near a zero denominator and not a
generic heterogeneity fraction. A prototype must report $F(\bar{x})$,
$\sum_kw_kF(x_k)$, $J$, and `R` together. Public promotion still requires the
pre-specified controlled and held-out validation in `PLAN.md`.

## Eligibility and missingness

Scientific eligibility requires compatible units and common preprocessing.
Structural equality cannot prove that counts, concentrations, normalized
abundances, or separately transformed values are scientifically commensurate;
that responsibility remains explicit caller metadata/policy.

Computational eligibility requires:

- at least one positive weight and at least two retained members;
- finite non-negative active-unit values after the chosen missingness policy;
- identical unique non-missing, non-empty member identifiers in identical order;
- one reported effective weight vector summing mathematically to one; and
- one explicit missingness policy, never imputation or missing-to-zero coercion.

The two allowed policies are:

1. `reject` (safe default): any `NA`/`NaN` in an active unit makes the group
   ineligible and reports the affected members/units.
2. `intersection` (explicit): retain only members observed in every active unit,
   preserve their original order, report every removed member identity/count,
   and reject fewer than two retained members. Recompute $n$, $c$, every unit
   score, $\bar{x}$, and the aggregate score on that same intersection.

Zero remains observed under both policies. Infinite or negative values remain
invalid. Missingness in zero-weight units is irrelevant because those units are
excluded first and reported.

## Mean, physical sum, and preprocessing

Positive homogeneity relates a compatible physical sum to the mean. For masses
$a_k\ge0$, $A=\sum_ka_k>0$, and $w_k=a_k/A$,

$$
F\!\left(\sum_ka_kx_k\right)-\sum_ka_kF(x_k)=AJ.
$$

This identity licenses rescaling only after the $x_k$ are compatible physical
quantities. Normalization, filtering, log transforms, detection thresholds, and
other preprocessing generally do not commute with pooling. `mean of processed
units`, `process a physical pool`, and `process a summed count vector` are
distinct estimands unless an external argument proves equivalence.

## Interpretation boundary

`Aggregation gap` is the field term. `Jensen gap` is mathematically accurate
but too generic for a biological field name. `Beta diversity`, `heterogeneity
score`, `complementarity`, `synergy`, `coherence`, and `Simpson effect` each
claim more than the theorem establishes.

The gap is descriptive, not causal attribution, uncertainty, or evidence that
different genes substitute biologically. It depends on weights, member support,
units, preprocessing, gene-specific dynamic ranges, missingness, and noise.
Cells nested within a donor remain subsamples: a gap can be calculated across
cells, but donor/experimental-unit replication governs inference.

## Prior-art audit

Audit date: 2026-07-18. Targeted, non-systematic primary-source/software audit;
not a patent search or exhaustive review.

| Area | Prior art | Relation and boundary |
|---|---|---|
| Jensen differences | [Burbea and Rao (1982)](https://doi.org/10.1109/TIT.1982.1056497) | Non-negative convex/entropy Jensen differences are established. $J$ is the sign-adjusted weighted Jensen difference of concave $F$; generic Jensen-gap novelty is unavailable. |
| Additive pooled diversity | [Lande (1996)](https://doi.org/10.2307/3545743); [Ricotta (2003)](https://doi.org/10.1023/A:1024539526618) | Pooled diversity minus mean within-unit diversity and concavity-based non-negativity are established alpha/beta partition ideas. GeneFunnel uses raw compatible member activity, not species diversity. |
| Independence criticism | [Jost (2007)](https://doi.org/10.1890/06-1736.1) | Additive beta components can retain dependence on alpha, especially with unequal community weights. `J` must not be sold as independent beta diversity or generic heterogeneity. |
| Pietra/Gini decomposition | [Porro and Zenga (2021)](https://doi.org/10.1007/s10182-021-00397-6); [Dagum (1997)](https://doi.org/10.1007/BF01205777) | Subgroup/source decompositions of mean-absolute-deviation and Gini inequality, including overlap/transvariation terms, are established. Here members are coordinates within compatible profiles; the exact opposing-centered-mass identity has a different estimand. |
| Aggregation reversal | [Simpson (1951)](https://doi.org/10.1111/j.2517-6161.1951.tb00088.x) | Collapsed and stratified associations can disagree. GeneFunnel's deterministic non-negative score difference has a fixed direction and no association or causal adjustment, so calling it Simpson's paradox would mislead. |
| Replicate-aware pseudobulk | [Crowell et al. (2020)](https://doi.org/10.1038/s41467-020-19894-4); [Zimmerman et al. (2021)](https://doi.org/10.1038/s41467-021-21038-1); [Squair et al. (2021)](https://doi.org/10.1038/s41467-021-25960-2) | Pseudobulk methods preserve biological-replicate structure for inference. They do not establish that nonlinear pathway scoring commutes with aggregation; this contract also forbids treating cells as independent donors. |
| Pathway-score benchmarks | [Zhang et al. (2020)](https://doi.org/10.1016/j.csbj.2020.10.007); [PathwayBench v3 software record (2026)](https://doi.org/10.5281/zenodo.19601134) | Single-cell and pseudobulk pathway-score accuracy/stability are active benchmark topics; the current software record explicitly includes aggregation stability. Neither source tests this GeneFunnel identity or supports a biological novelty claim. |

## Executable evidence and next gate

`tests/testthat/helper-aggregation-reference.R` implements the coordinate
formula and opposing-mass formula in base R, independently of native code.
`test-aggregation-theorem.R` checks randomized weighted identities, bounds,
normalization, equality/strictness, zero-weight exclusion, scaling, physical
sums, and unit/member permutations.

A future group-audit API must be additive and must return eligibility reasons,
effective weights/support, aggregate score, weighted unit score, absolute gap,
normalized gap/status, and removed member/unit facts. Its schema and factorial
validation thresholds are frozen in benchmark protocol `B-1.0.1`;
the prototype and its empirical evidence remain prospective.
