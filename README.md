
# GeneFunnel

GeneFunnel calculates sample-wise gene-set scores from non-negative
feature-by-sample matrices. Each score preserves observed magnitude while
discounting disagreement among the available members. Dense inputs use a dense
path; sparse inputs remain sparse and retain implicit zeros as observations.

GeneFunnel is under development and is not yet available from Bioconductor.

## Why not use a plain sum?

A member sum cannot distinguish coherent activity from the same total driven by
only one or two features. GeneFunnel can:

```r
library(genefunnel)

same_total <- rbind(
    A = c(coherent = 3, split = 6, concentrated = 12),
    B = c(coherent = 3, split = 6, concentrated = 0),
    C = c(coherent = 3, split = 0, concentrated = 0),
    D = c(coherent = 3, split = 0, concentrated = 0)
)

scores <- genefunnel(
    same_total,
    list(pathway = rownames(same_total)),
    BPPARAM = BiocParallel::SerialParam()
)
rbind(
    observed_sum = colSums(same_total),
    GeneFunnel = unname(scores["pathway", ])
)
#>              coherent split concentrated
#> observed_sum       12    12           12
#> GeneFunnel         12     4            0
```

All three samples have the same observed sum. GeneFunnel gives full credit to
the coherent pattern, partial credit to the split pattern, and zero to the
single-feature pattern. Use this score when total magnitude and within-set
agreement both carry meaning. Prefer a rank-based method when scale invariance
is the intended property. This contrast demonstrates the score's behavior; it
does not establish biological superiority for every analysis.

## Installation

Install the development version from GitHub:

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}
BiocManager::install("eturkes/genefunnel")
```

After acceptance, install the Bioconductor package with:

```r
BiocManager::install("genefunnel")
```

## Coverage-aware workflow

```r
library(genefunnel)

mat <- rbind(
    A = c(sample_1 = 4, sample_2 = 1),
    B = c(sample_1 = 4, sample_2 = 2),
    C = c(sample_1 = 4, sample_2 = NA)
)
gene_sets <- list(
    coherent = c("A", "B", "C"),
    partial = c("A", "B", "D"),
    insufficient = c("A", "Z")
)

coverage <- gene_set_coverage(gene_sets, rownames(mat))
coverage

# Example caller policy: at least 50% coverage and two matched members.
keep <- !is.na(coverage$coverage) &
    coverage$coverage >= 0.5 &
    coverage$matched_size >= 2

genefunnel(
    mat,
    gene_sets[keep],
    BPPARAM = BiocParallel::SerialParam()
)
#>          sample_1 sample_2
#> coherent       12        2
#> partial         8        2
```

The second sample omits `C = NA` and scores the observed pair `c(1, 2)`.
Zeros, including implicit sparse zeros, remain in the calculation. A set/sample
cell with fewer than two observed members returns `NA`.

Sparse matrices and portable parallel backends use the same API:

```r
sparse_mat <- Matrix::Matrix(mat, sparse = TRUE)
backend <- BiocParallel::SnowParam(workers = 2L, type = "SOCK")
scores <- genefunnel(sparse_mat, gene_sets[keep], BPPARAM = backend)
```

Matching is exact and case-sensitive. GeneFunnel performs no normalization,
imputation, identifier mapping, gene-set retrieval, downstream testing, or
biological interpretation. Input must be non-negative with a meaningful zero;
an arbitrary positive shift changes the score and is not a valid workaround.

## Validation evidence

The vignette evaluates the equal-total contrast above as part of the package
build. The regression suite also compares 32 deterministic dense and sparse
fixtures against an independent R reference implementation. Additional tests
verify score bounds, mathematical invariants, missing-value semantics, partial
coverage, storage equivalence, and serial/SOCK ordering. This evidence supports
the documented calculation and implementation; it is not a benchmark of
biological usefulness against other methods.

## Documentation

- [Scientific specification](inst/SCIENTIFIC_SPEC.md) - equation, value and
  coverage semantics, invariants, and limitations.
- [Vignette](vignettes/genefunnel.Rmd) - canonical examples, dense/sparse use,
  coverage policies, and serial/parallel execution.
- [Tests](tests/testthat) - deterministic formula, coverage, dense/sparse,
  serial/parallel, randomized-oracle, and validation checks.

After installation, run `vignette("genefunnel", package = "genefunnel")` or
`citation("genefunnel")`.

## License

GeneFunnel is licensed under GPL-3-or-later. See [COPYING](COPYING).
