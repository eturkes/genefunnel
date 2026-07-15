<!-- Assisted-by: OpenAI Codex. -->

# GeneFunnel scientific specification

This document is the durable public contract for GeneFunnel scoring. It is
derived from Sections 4.3-4.12 and 5.1 of Emir Turkes' 2025 University College
London PhD thesis, *Development of Gene Set Enrichment and Imputation Methods
for Transcriptomics and Proteomics: Application in the Study of
Neurofibrillary Tangle-bearing Neurons in Alzheimer's Disease*.

## Scope

GeneFunnel performs functional class scoring. It converts a non-negative
feature-by-sample matrix into a numeric gene-set-by-sample matrix. It does not
normalize or impute measurements, map identifiers, retrieve gene sets, perform
differential testing, or provide biological interpretation.

`genefunnel(mat, gene_sets, BPPARAM)` accepts:

- an unclassed base R numeric/integer matrix or supported numeric dense/sparse
  `Matrix`;
- an unclassed named list of unclassed character gene-set members; and
- a `BiocParallelParam`, defaulting to the registered backend.

The result is a base numeric matrix. Retained gene sets and samples preserve
their input order and names. Missing column names remain `NULL`. One sample and
one gene set are valid.

## Exact score

For input matrix \(X\), declared gene set \(G\), and sample \(j\):

1. Deduplicate \(G\) as a mathematical set, preserving first occurrence.
2. Match identifiers exactly and case-sensitively against `rownames(X)`. Let
   \(M = G \cap \operatorname{rownames}(X)\).
3. For sample \(j\), remove `NA` and `NaN` values from the matched members. Let
   the remaining observed values be \(x_1,\ldots,x_n\).
4. Return `NA_real_` when \(n < 2\).
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

The effective \(n\) is specific to a gene-set/sample cell after missing-value
omission. It is neither the declared set size nor necessarily the globally
matched set size.

## Values and missingness

Zero is observed information. Explicit dense zeros and implicit sparse zeros
participate in the effective size, sum, mean, and deviation.

`NA` and `NaN` are absent information. They are omitted per sample and reduce
the effective size. Missingness in one sample has no effect on another sample.

Positive and negative infinity are invalid. Negative finite values are also
invalid because the activity interpretation and score bounds require
non-negative input. GeneFunnel rejects these values and never shifts,
rescales, normalizes, or log-transforms input automatically.

Processed data are suitable only when zero remains meaningful and every
observed value is non-negative. An arbitrary positive shift is not a valid
workaround: adding \(c\) to every observed member increases the score by
\(nc\).

Three cases remain distinct:

- a matrix row that is absent is globally unmeasured and affects coverage;
- an `NA`/`NaN` cell is unavailable only in that sample and affects its
  effective size; and
- a zero is measured inactivity and remains in the calculation.

## Coverage and identifiers

GeneFunnel scores every set with at least two unique members matched globally
to matrix row names. Unmatched declared members are ignored during scoring.
Sets below the two-match minimum are omitted from the result in input order,
with at most one aggregate warning per call. The scoring API has no fixed
coverage threshold.

`gene_set_coverage(gene_sets, features)` applies the scorer's exact validation,
deduplication, and matching rules. It returns:

- gene-set name;
- unique declared size;
- matched and unmatched sizes;
- matched fraction of unique declared members;
- duplicate-member count; and
- whether at least two members are globally scoreable.

Coverage is undefined (`NA`) for an empty declared set. Empty sets are valid
diagnostic inputs but cannot be scored. Experimental policy belongs to the
caller. Examples include:

- strict transcriptomic policy: `coverage == 1`; and
- low-coverage proteomic policy: `coverage >= 0.5` and `matched_size >= 2`.

Matrix row names are required, unique, non-missing, and non-empty. Gene-set
names are required, unique, non-missing, and non-empty. Every member identifier
must be a non-missing, non-empty string. GeneFunnel performs no case folding,
synonym resolution, annotation lookup, or identifier conversion. Overlapping
gene sets are valid and remain independent.

## Canonical values

| Observed values | Score | Explanation |
|---|---:|---|
| `c(4, 4, 4)` | `12` | No deviation; score equals the sum. |
| `c(0, 0, 0)` | `0` | All members are measured as inactive. |
| `c(4, 0)` | `0` | Maximally deviant two-member set. |
| `c(4, 0, 0)` | `0` | Only one member is non-zero. |
| `c(1, 2, NA)` | `2` | Omit `NA`; score the two observed values. |
| `c(4, NA)` | `NA` | Fewer than two observed members remain. |
| sparse implicit `c(4, 0)` | `0` | The implicit zero is observed. |

For a declared set `c("A", "B", "C")` and matrix rows containing only `A`
and `B`, GeneFunnel scores the available pair, reports coverage `2/3`, and
leaves the acceptance threshold to the caller.

## Required properties

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

Numerical cleanup may turn a theoretically zero negative round-off result into
exact zero only within a scale-aware tolerance. A materially negative or
non-finite result is an internal error and must not be hidden.

## Scientific interpretation and limits

The score combines total observed activity with a scaled absolute-deviation
penalty. Treating within-set variability as evidence against coherent activity
is a scientific assumption, not a universal biological rule.

The sum term retains input magnitude. Scores therefore retain effects of input
units, library size, preprocessing, gene-set size, and feature coverage.
Sample-wise independence does not make unrelated data sets automatically
comparable. Cross-data-set comparison requires compatible units,
preprocessing, identifiers, gene-set definitions, and coverage policies.

Gene-set database revisions can change scores, so analyses should record the
exact collection and version. GeneFunnel does not guarantee that score
distributions satisfy every downstream statistical model; users must validate
the assumptions of their chosen analyses.

Raw non-negative measurements can be preferable where scientifically
appropriate because prior filtering or transformation can remove zeros or
alter magnitude. The package does not prescribe or enforce a pipeline position.

## Computational equivalence

Dense input follows a dense native path. Sparse input remains sparse and
accounts analytically for implicit zeros. Work is split into bounded column
chunks for BiocParallel execution. These storage and scheduling choices are
implementation details: they must preserve the equation, missing-value rules,
names, and ordering above.
