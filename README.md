
# GeneFunnel

GeneFunnel scores predefined gene sets in each sample of a non-negative
feature-by-sample matrix. The score starts with the total observed signal and
reduces it when members of the set disagree. Dense and sparse matrices use the
same scoring rules, and implicit sparse zeros count as observations.

GeneFunnel is under development and is not yet available from Bioconductor.

## What the score captures

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

A plain sum gives all three samples a value of 12. GeneFunnel gives 12 to the
coherent pattern, 4 to the split pattern, and 0 to the single-feature pattern.
Use it when both the total signal and agreement within the set matter. Choose a
rank-based method when relative ordering matters more than measurement scale.

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

## Check coverage before scoring

Use `gene_set_coverage()` to see which members are present before choosing a
coverage threshold for your experiment.

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

# Example policy: at least 50% coverage and two matched members.
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

For `sample_2`, GeneFunnel omits `C = NA` and scores the observed pair
`c(1, 2)`. Zeros remain in the calculation. A score is `NA` when fewer than two
members are observed in that sample.

## Sparse matrices and parallel processing

Sparse matrices and portable parallel backends use the same API:

```r
sparse_mat <- Matrix::Matrix(mat, sparse = TRUE)
backend <- BiocParallel::SnowParam(workers = 2L, type = "SOCK")
scores <- genefunnel(sparse_mat, gene_sets[keep], BPPARAM = backend)
```

## Input assumptions

Matching is exact and case-sensitive. GeneFunnel does not normalize or impute
data, map identifiers, retrieve gene sets, run statistical tests, or interpret
the result. Input values must be non-negative, and zero must have a meaningful
interpretation. Adding a constant to make negative data positive changes the
score and is not an appropriate workaround.

## How the package is tested

The test suite compares 32 deterministic dense and sparse examples with an
independent R implementation. It also checks score bounds, missing-value rules,
partial coverage, dense/sparse agreement, and serial/parallel agreement. These
tests verify the documented calculation; they do not measure biological
performance against other methods.

## Further reading

- [Scoring specification](inst/SCIENTIFIC_SPEC.md) - equation, value and
  coverage rules, properties, and limitations.
- [Vignette](vignettes/genefunnel.Rmd) - worked examples, coverage checks,
  dense and sparse matrices, and parallel execution.
- [Tests](tests/testthat) - reference calculations, coverage checks, storage
  agreement, parallel agreement, and input validation.

After installation, run `vignette("genefunnel", package = "genefunnel")` or
`citation("genefunnel")`.

## License

GeneFunnel is licensed under GPL-3-or-later. See [COPYING](COPYING).
