
# GeneFunnel scientific specification

This specification defines score calculation and interpretation. The method
appears in Sections 4.3-4.12 and 5.1 of Emir Turkes' 2025
University College London PhD thesis, *Development of Gene Set Enrichment and
Imputation Methods for Transcriptomics and Proteomics: Application in the Study
of Neurofibrillary Tangle-bearing Neurons in Alzheimer's Disease*.

## Scope

GeneFunnel converts a non-negative feature-by-sample matrix into a numeric
gene-set-by-sample matrix through functional class scoring. It does not
normalize, impute, map identifiers, retrieve sets, test differences, or
interpret results biologically.

`genefunnel(mat, gene_sets, BPPARAM, min_coverage, min_observed, details)` accepts:

- an unclassed base R numeric/integer matrix or supported numeric dense/sparse
  `Matrix`;
- an unclassed named list of unclassed character gene-set members; and
- a `BiocParallelParam`, defaulting to the registered backend.

The default base numeric matrix preserves retained set and sample order and
names. Missing column names stay `NULL`. Single samples and sets are valid.

With `details = TRUE`, results include observed counts, sample coverage, and
global set coverage. Coverage filters leave the score equation unchanged.

## Score calculation

For input matrix \(X\), declared gene set \(G\), and sample \(j\):

1. Deduplicate \(G\) as a mathematical set, preserving first occurrence.
2. Match identifiers exactly and case-sensitively against `rownames(X)`. Let
   \(M = G \cap \operatorname{rownames}(X)\).
3. For sample \(j\), remove `NA` and `NaN` values from the matched members. Let
   the remaining observed values be \(x_1,\ldots,x_n\).
4. Return `NA_real_` when \(n < \texttt{min_observed}\) or
   \(n / |G| < \texttt{min_coverage}\). The defaults are two and zero.
5. Otherwise calculate

\[
\bar{x}=\frac{1}{n}\sum_{i=1}^{n}x_i
\]

and

\[
\operatorname{GeneFunnel}(x_1,\ldots,x_n)
=\sum_{i=1}^{n}x_i
-\frac{n}{2(n-1)}\sum_{i=1}^{n}\left|x_i-\bar{x}\right|.
\]

The effective \(n\) counts observations in one set/sample cell, not declared
or globally matched members.

## Values and missingness

Explicit and implicit sparse zeros contribute to the observed count, sum,
mean, and deviation.

`NA` and `NaN` reduce the observed count in their sample only.

Infinities and negative values are rejected: interpretation and score bounds
require non-negative input. GeneFunnel never shifts, rescales, normalizes, or
log-transforms input automatically.

Processed data need meaningful zeros and non-negative observations. Arbitrary
positive shifts are unsuitable: adding \(c\) to every observation adds \(nc\)
to the score.

Three cases remain distinct:

- an absent row reduces global coverage;
- an `NA`/`NaN` cell reduces observations in one sample; and
- a zero remains observed in the calculation.

## Coverage and identifiers

Sets need two unique global matches to matrix row names. Scoring ignores
unmatched members; sets below the minimum are omitted with at most one
aggregate warning per call. Retained sets stay in input order.

No fractional threshold applies by default. Set `min_coverage` to require a
sample-level fraction of unique declared members, or `min_observed` for more
than two observations. Failed cells become `NA`; their rows stay.

`gene_set_coverage(gene_sets, features)` uses the scorer's validation,
deduplication, and matching rules. It returns set names, unique declared sizes,
matched and unmatched counts, matched fractions, duplicate counts, and whether
at least two members match.

Empty sets are valid diagnostic inputs with `NA` coverage, but cannot be scored.
Choose experiment-appropriate thresholds, for example:

- strict transcriptomic policy: `coverage == 1`; and
- low-coverage proteomic policy: `coverage >= 0.5` and `matched_size >= 2`.

Row names and set names must be unique, non-missing, and non-empty. Members
must be non-missing, non-empty strings. Matching performs no case folding,
synonym resolution, annotation lookup, or identifier conversion. Overlapping
sets remain independent.

## Worked values

| Observed values | Score | Explanation |
|---|---:|---|
| `c(4, 4, 4)` | `12` | No deviation; score equals the sum. |
| `c(0, 0, 0)` | `0` | All observed values are zero. |
| `c(4, 0)` | `0` | Maximally deviant two-member set. |
| `c(4, 0, 0)` | `0` | Only one member is non-zero. |
| `c(1, 2, NA)` | `2` | Omit `NA`; score the two observed values. |
| `c(4, NA)` | `NA` | Fewer than two observed members remain. |
| sparse implicit `c(4, 0)` | `0` | The implicit zero is observed. |

For declared members `c("A", "B", "C")` with only `A` and `B` present,
GeneFunnel scores the pair and reports coverage `2/3`. The user sets acceptance
thresholds.

## Properties of the score

For finite, observed, non-negative values:

- the score is non-negative within floating-point tolerance and does not
  exceed the observed sum;
- equal values score their sum;
- all-zero values and a set with exactly one non-zero value score zero;
- member order does not affect the score;
- multiplying all observed values by \(c \ge 0\) multiplies the score by \(c\);
- adding \(c \ge 0\) to all \(n\) observed values adds \(nc\) to the score;
- other samples, other gene sets, and irrelevant matrix rows do not affect the
  cell being scored;
- dense and sparse representations agree; and
- serial and parallel execution agree and preserve ordering.

Numerical cleanup may set negative round-off to zero only within a scale-aware
tolerance. Materially negative or non-finite results are internal errors and
must not be hidden.

## Relation to inequality measures

For positive total \(S=\sum_i x_i\), define the Pietra (Hoover) index
\(H=\sum_i|x_i-\bar{x}|/(2S)\). Then
\(\operatorname{GeneFunnel}=S[1-nH/(n-1)]\). The all-zero case returns zero
without evaluating this ratio. This is a finite-size Pietra-penalized total,
not a new inequality index or the pairwise-difference Gini coefficient.

For fixed observed size, the score is coordinatewise nondecreasing and concave.
It can become insensitive to further increases in one dominant member.
For example, \(\operatorname{GeneFunnel}(1,1,M)=3\) for every \(M\ge1\).
These properties do not extend to changing the observed member count.

Within-sample equality is not cross-sample correlation or coexpression.
Equal values need not imply pathway activation. Unequal expected gene
abundances can produce a penalty even when all genes respond together.

Pietra definition: Porro F, Zenga M. *Decompositions by sources and by
subpopulations of the Pietra index: two applications to professional football
teams in Italy*. AStA Advances in Statistical Analysis 107:73–100.
<https://doi.org/10.1007/s10182-021-00397-6>.

## Interpretation and limitations

The score penalizes absolute deviation from the mean. Treating variability as
evidence against coherent activity is a scientific assumption, not a universal
biological rule.

The sum retains effects of units, library size, preprocessing, set size, and
coverage. Independent scoring does not make unrelated datasets comparable.
Comparisons need compatible units, preprocessing, identifiers, set definitions,
and coverage policies.

Database revisions can change scores. Record the exact collection and version,
and validate score distributions against downstream model assumptions.

Raw measurements may be preferable when filtering or transformation removes
zeros or alters magnitude. The package prescribes no pipeline position.

## Dense, sparse, and parallel calculations

Dense input stays dense; sparse input stays sparse, accounting analytically
for implicit zeros. BiocParallel processes bounded column chunks. Storage and
scheduling must preserve the equation, missing-value rules, names, and order.
