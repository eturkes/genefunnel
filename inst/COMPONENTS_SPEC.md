<!-- Assisted-by: OpenAI Codex. -->

# GeneFunnel component contract and proof note

**Status:** ratified pre-implementation contract, 2026-07-17. No component API
exists yet. This document constrains the additive prototype; it does not change
`genefunnel()` or `gene_set_coverage()`.

The frozen scorer remains governed by `SCIENTIFIC_SPEC.md`. Terms here apply to
the exact matched, observed gene-set/sample vector used by that scorer.

## Definitions and factorization

Let $x=(x_1,\ldots,x_n)$ contain the finite, observed, non-negative values for
one retained gene set and sample, with $n\ge2$, total
$T=\sum_i x_i$, mean $\bar{x}=T/n$, and GeneFunnel score

$$
F(x)=T-\frac{n}{2(n-1)}\sum_i|x_i-\bar{x}|.
$$

For $T>0$, define shares $p_i=x_i/T$, uniform shares $u_i=1/n$,
and total-variation distance

$$
D(x)=\operatorname{TV}(p,u)=\frac12\sum_i|p_i-u_i|.
$$

Since $x_i-\bar{x}=T(p_i-1/n)$,

$$
\sum_i|x_i-\bar{x}|=2TD(x)
$$

and therefore

$$
F(x)
=T\left[1-\frac{n}{n-1}D(x)\right]
=TB(x),
\qquad
B(x)=1-\frac{n}{n-1}D(x).
$$

This $B$ is exactly [Bulla's evenness index](https://doi.org/10.2307/3545713)
applied to member shares. `balance` remains the GeneFunnel field name because
its intended role is arithmetic explanation, not ecological interpretation.

The four ordinary components are:

- magnitude: $T$, named `observed_sum`;
- balance: $B$, the complement of finite-size-normalized inequality;
- arithmetic penalty: $Q=T-F=T(1-B)$; and
- effective size: $n$, after sample-specific missing-value omission.

For any probability vector on $n$ coordinates,

$$
D(p,u)=1-\sum_i\min(p_i,1/n)\le 1-1/n.
$$

At least one share is at least $1/n$, which proves the inequality. Equality
holds at a one-point mass. Hence $0\le B\le1$, $0\le F\le T$, equal
positive values have $B=1$, and exactly one positive value has $B=0$.

When $T=0$, non-negativity implies that every observed value is zero. The
authoritative score and penalty are zero, but shares and balance are undefined.
Balance is `NA`, not one.

## Structural properties

For fixed $n$, let
$H=I-\mathbf{1}\mathbf{1}^{\mathsf T}/n$ and
$c=n/[2(n-1)]$. Then

$$
F(x)=\mathbf{1}^{\mathsf T}x-c\lVert Hx\rVert_1.
$$

- **Symmetry:** coordinate permutations leave the sum and norm unchanged.
- **Positive homogeneity:** for $a>0$, $F(ax)=aF(x)$, $T(ax)=aT(x)$,
  $Q(ax)=aQ(x)$, and $B(ax)=B(x)$. At $a=0$, score and penalty are zero
  while balance becomes undefined.
- **Concavity:** the first term is affine and the negative of an L1 norm of a
  linear map is concave. Thus
  $F(wx+(1-w)y)\ge wF(x)+(1-w)F(y)$ for $0\le w\le1$, common coordinates,
  and common $n$.

Concavity does not license mixing vectors with different observed support,
units, preprocessing, or member order. Sample-specific missingness can change
the dimension, so it must be resolved before an aggregation comparison.

## Cell and output semantics

The future component function must reuse the scorer's prepared memberships,
retained-set universe, aggregate omission warning, row/sample order, and names.
Sets with fewer than two globally matched members remain absent and are
reported only by `gene_set_coverage()`.

For every retained set/sample cell:

| Field | Fewer than 2 observed | At least 2, all zero | At least 2, positive total |
|---|---:|---:|---:|
| `effective_size` | factual count | factual count | factual count |
| `observed_sum` | factual sum | `0` | $T$ |
| `observed_fraction` | factual fraction | factual fraction | factual fraction |
| `score` | `NA` | `0` | authoritative $F$ |
| `penalty` | `NA` | `0` | $Q$ |
| `balance` | `NA` | `NA` | $B\in[0,1]$ |

`observed_fraction = effective_size / globally_matched_size`. Its denominator
is the unique declared members matched to matrix rows before cell missingness.
This differs from declared-set coverage,
`globally_matched_size / unique_declared_size`, returned by
`gene_set_coverage()`.

`NA` and `NaN` remain absent; explicit and implicit zeros remain observed.
Components never cause a second interpretation of membership, missingness, or
zero.

## Binary64 representation contract

The core score is authoritative. Diagnostics must neither make the ordinary
`genefunnel()` path fail nor replace a representable core score with a value
reconstructed from components. If the score itself is not representable, the
existing core error remains authoritative.

Finite inputs can still create diagnostics outside binary64:

- $T$ and $Q$ can overflow while $F$ remains representable;
- finite rounded $T$ and $Q$ can become equal although $F>0$, so
  binary64 subtraction cannot reconstruct the score; and
- a true positive $B=F/T$ can underflow.

The additive API must therefore return aligned semantic, availability, and
conditioning status rather than `Inf` or a silent false zero:

- semantic state = `scoreable`, `too_few_observed`, or `zero_total`;
- availability = `ordinary`, `scaled`, or `unavailable`; and
- conditioning = `safe`, `ill_conditioned`, or `not_applicable`.

Every defined ordinary diagnostic is a finite double. A mathematically defined
nonzero value that overflows or underflows has ordinary value `NA`, never zero
or infinity. `scaled` means a sidecar supplies a pair
$(m,e)$ with value $m\,2^e$, $m=0$ or $1/2\le m<1$, and integer $e$.
Undefined and unavailable values use an `NA` pair. `unavailable` means no
trustworthy ordinary or scaled diagnostic could be certified; it never changes
the score.

For $F>0$, direct subtraction has condition number

$$
\kappa_{-}=\frac{T+Q}{F}.
$$

A cell is `ill_conditioned` whenever rounding-error bounds cannot certify both
$T-Q$ and $TB$ as reconstructions of the authoritative score. The A2
protocol must commit a concrete safe region and tolerance before testing native
outputs. Outside that region, scaled identities and status - not direct
binary64 equality - are required.

## Interpretation cases

These cases are committed before implementation:

| Observed values | Score | Sum | Balance | Effective size | Point |
|---|---:|---:|---:|---:|---|
| `c(2, 2)` | 4 | 4 | 1 | 2 | equal positive pair |
| `c(6, 2)` | 4 | 8 | 0.5 | 2 | same score, different magnitude/balance |
| `c(4/3, 4/3, 4/3)` | 4 | 4 | 1 | 3 | same score/sum/balance, different support |
| `c(4, 0)` | 0 | 4 | 0 | 2 | one positive member |
| `c(0, 0)` | 0 | 0 | `NA` | 2 | zero total is not perfect balance |
| `c(1, 2, NA)` from 3 matched members | 2 | 3 | 2/3 | 2 | observed fraction is 2/3 |

The scalar score alone cannot distinguish the first three rows. The components
describe the arithmetic distinction; they do not identify its biological cause.

## Interpretation boundary

Balance means evenness of non-negative within-set shares. It is not
co-expression, correlation, coordinated regulation, pathway truth, causal gene
attribution, uncertainty, or preprocessing invariance. It depends on effective
size, retained member identity, gene-specific abundance and dynamic range,
measurement noise, and the meaningful-zero assumption. Equal balance does not
make separately processed data comparable.

Bulla's index is known to depend on the number of represented categories and on
whether zero-share categories remain in the support. That is a limitation here,
not a reason to delete observed zero-valued members: the frozen GeneFunnel
estimand includes them.

The word `coherence` is intentionally excluded because pathway literature often
uses it for cross-sample correlation or coordinated variation. `variability`
is also reserved for across-unit dispersion unless explicitly qualified.

## Prior-art audit

Audit date: 2026-07-17. This was a targeted, non-systematic audit of primary
mathematical/ecological papers and current official gene-set software
documentation. It establishes terminology and blocks unsupported novelty
claims; it is not a patent search or exhaustive systematic review.

| Area | Prior art | Relation and boundary |
|---|---|---|
| Mean absolute deviation / Pietra | [Eliazar and Sokolov (2010)](https://doi.org/10.1016/j.physa.2009.08.006); [Porro and Zenga (2021)](https://doi.org/10.1007/s10182-021-00397-6) | Pietra's index is mean absolute deviation about the mean divided by twice the mean. For equal member weights it equals $D(x)=\operatorname{TV}(p,u)$, with maximum $(n-1)/n$. |
| Hoover index | [Hoover (1936)](https://doi.org/10.2307/1927875) | The equal-population-share Hoover form is the same half-L1 redistribution fraction. Its economic/geographic interpretation does not transfer automatically to genes. |
| Exact ecological index | [Bulla (1994)](https://doi.org/10.2307/3545713) | With overlap $O=\sum_i\min(p_i,1/n)=1-D$, Bulla's $E=(O-1/n)/(1-1/n)$ is exactly $B=1-D/(1-1/n)$. The index and normalization are established prior art. |
| Ecological distinctions and criticism | [Camargo (1992)](https://doi.org/10.1007/BF00195643); [Camargo (1995)](https://doi.org/10.2307/3546000); [Molinari (1996)](https://doi.org/10.2307/3546352); [Gregorius and Gillet (2021)](https://doi.org/10.1007/s10441-021-09429-9) | Above-uniform dominance mass equals $D$, while Camargo's named pairwise evenness index is different from Bulla's. Published criticism shows that Bulla's value depends on richness/support; ecological names and interpretations must not be conflated or transferred to genes. |
| Rank score plus dispersion | [Foroutan et al. (2018)](https://doi.org/10.1186/s12859-018-2435-4); [current `singscore` API](https://bioconductor.org/packages/release/bioc/manuals/singscore/man/singscore.pdf) | `simpleScore()` returns sample scores and dispersions; `multiScore()` returns score and dispersion matrices and accepts `dispersionFun = mad`. This is a strong API precedent, but it operates on within-sample ranks and its dispersion is not an algebraic factor of GeneFunnel. |
| Enrichment diagnostics | [Hänzelmann et al. (2013)](https://doi.org/10.1186/1471-2105-14-7); [current GSVA API](https://bioconductor.org/packages/release/bioc/manuals/GSVA/man/GSVA.pdf) | GSVA returns gene-set/sample enrichment scores; current `gsvaRanks()` and `gsvaEnrichment()` expose intermediate enrichment data for a selected cell. They do not return this magnitude/balance identity. |
| Pathway variability | [PAGODA](https://doi.org/10.1038/nmeth.3734); [pathVar](https://doi.org/10.7717/peerj.3334) | PAGODA tests coordinated overdispersion across cells; pathVar bins per-gene variability across samples and decomposes pathways into bin counts. Both address across-unit heterogeneity, not within-cell share evenness. |
| SVD pathway activity | [Tomfohr et al. (2005)](https://doi.org/10.1186/1471-2105-6-225) | PLAGE derives pathway profiles from a cross-sample SVD. It is neither sample-independent nor an exact decomposition of this score. |

**Novelty boundary:** neither $D$ nor $B$ is a new inequality/evenness index;
$B$ already has the name Bulla's evenness index. Any future claim is restricted
to the exact GeneFunnel factorization, its implementation and auditing value,
and biological utility demonstrated under pre-specified held-out validation.

## Executable evidence

`tests/testthat/helper-reference.R` contains an ordinary-range independent R
oracle based on shares and total variation. `test-components-theorem.R` checks
the factorization, Pietra/Bulla equivalence, bounds, homogeneity, concavity,
zero semantics, missingness facts, and the committed interpretation cases.
Extreme scaled arithmetic belongs to A2's higher-precision/scaled oracle.
