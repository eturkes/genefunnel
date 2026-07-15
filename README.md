<!-- Assisted-by: OpenAI Codex. -->

# GeneFunnel

GeneFunnel calculates sample-wise gene-set scores from non-negative
feature-by-sample matrices. Each score combines observed magnitude with a
scaled absolute-deviation penalty. Dense inputs use a dense path; sparse inputs
remain sparse and retain implicit zeros as observations.

GeneFunnel is under development and is not yet available from Bioconductor.

## Installation

Install the development version from GitHub:

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}
BiocManager::install("eturkes/genefunnel")
```

## Quick start

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

## Documentation

- [Scientific specification](inst/SCIENTIFIC_SPEC.md) - equation, value and
  coverage semantics, invariants, and limitations.
- [Vignette](vignettes/genefunnel.Rmd) - canonical examples, dense/sparse use,
  coverage policies, and serial/parallel execution.
- [Benchmark protocol](https://github.com/eturkes/genefunnel/blob/main/benchmark/README.md)
  - controlled scientific assertions plus reproducible timing, memory,
  output-identity, and environment evidence.

After installation, run `vignette("genefunnel", package = "genefunnel")` or
`citation("genefunnel")`.

## License

GeneFunnel is licensed under GPL-3-or-later. See [COPYING](COPYING).
