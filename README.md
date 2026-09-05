# GeneFunnel

GeneFunnel scores predefined gene sets in each sample of a non-negative
feature-by-sample matrix. The score starts with the total observed signal and
reduces it when members of the set disagree. Dense and sparse matrices use the
same scoring rules, and implicit sparse zeros count as observations.

Use it when you need to distinguish a coherent gene set from one driven by a few
large values. The output is a gene-set-by-sample score matrix, with optional
coverage diagnostics. It does not normalize data, test enrichment, or establish
biological activity. A rank-based method addresses a different question when
relative ordering matters more than measurement scale.

## Installation

The development package is not yet available from Bioconductor. Install it from
GitHub in R. Source installation needs a C++17 compiler; the example figure also
needs `ggplot2`.

```r
install.packages(c("BiocManager", "remotes", "ggplot2"))
BiocManager::install("eturkes/genefunnel")
```

## First result

Run this from a writable working folder after installation. It uses a built-in
synthetic example, downloads no data, and saves the figure shown below.

```r
source(system.file("examples", "quickstart.R", package = "genefunnel",
    mustWork = TRUE))
```

![Three four-gene inputs each total 12. GeneFunnel scores the coherent, split, and concentrated patterns as 12, 4, and 0.](man/figures/quickstart.png)

All three samples total 12. GeneFunnel returns **12, 4, and 0** as the signal
becomes concentrated in fewer members. Zeros are observations, not missing
values. The final score of zero does not prove biological inactivity.

The script prints a table and writes these files under `genefunnel-results/`:

| File | Contents |
|---|---|
| `comparison.png` | The labelled figure above, at 180 dpi |
| `scores.csv` | Each sample's total, score, observed count, and coverage |
| `gene-set-coverage.csv` | Declared and matched set sizes |
| `result.rds` | The full result, for `readRDS()` |

Re-running replaces these example files. The script leaves `x`, `gene_sets`,
`result`, `score_table`, and `score_plot` in your R session. Read the
[example script](inst/examples/quickstart.R) to change the inputs or figure.

## Score your data

Replace `x` with your non-negative feature-by-sample matrix. Give its rows
feature identifiers and its columns sample names. Supply `gene_sets` as a named
list of character vectors, using the same identifiers.

```r
library(genefunnel)

result <- genefunnel(
    x, gene_sets,
    BPPARAM = BiocParallel::SerialParam(),
    details = TRUE
)
result$scores
result$coverage
```

Inspect coverage alongside scores. Missing members can lower a score without
showing a biological difference. Compare samples on a justified common input
scale; GeneFunnel does not correct different measurement scales or set sizes.

### Check coverage before scoring

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

### Enforce sample-level coverage

Global matching does not detect missing observations within a sample. Apply
sample-level thresholds when those gaps affect interpretation:

```r
checked <- genefunnel(
    mat,
    gene_sets[keep],
    BPPARAM = BiocParallel::SerialParam(),
    min_coverage = 0.8,
    min_observed = 2L,
    details = TRUE
)
checked$scores
checked$observed
checked$coverage
checked$set_coverage
```

Coverage divides the observed count by the **unique declared set size**.
Unmatched members and missing cells reduce coverage; zeros count as observations.
Cells below either threshold return `NA`, without removing their gene-set row.
The thresholds do not rescale scores or correct missingness bias.

The defaults are `min_coverage = 0`, `min_observed = 2L`, and `details = FALSE`.
They preserve the original score-matrix API and numerical results.

### Sparse matrices and parallel processing

Sparse matrices and portable parallel backends use the same API:

```r
sparse_mat <- Matrix::Matrix(mat, sparse = TRUE)
backend <- BiocParallel::SnowParam(workers = 2L, type = "SOCK")
scores <- genefunnel(sparse_mat, gene_sets[keep], BPPARAM = backend)
```

### Experiment containers

Pass a `SummarizedExperiment` directly and select the assay to score:

```r
se <- SummarizedExperiment::SummarizedExperiment(assays = list(signal = mat))
checked <- genefunnel(
    se,
    gene_sets[keep],
    assay = "signal",
    BPPARAM = BiocParallel::SerialParam(),
    details = TRUE
)
```

The only assay is selected when `assay` is omitted. Multiple assays require
one exact, unambiguous name. Feature matching uses the container's row names;
the output follows its sample names and order.

The selected assay must be a supported numeric matrix. Sparse assays remain
sparse; delayed and disk-backed assays are rejected rather than materialized.
Other assays and metadata stay unchanged. The return type remains a score
matrix or diagnostic list, not an experiment container.

## Input assumptions

Matching is exact and case-sensitive. GeneFunnel does not normalize or impute
data, map identifiers, retrieve gene sets, run statistical tests, or interpret
the result. Input values must be non-negative, and zero must have a meaningful
interpretation. Adding a constant to make negative data positive changes the
score and is not an appropriate workaround.

## Verification

Tests compare the native implementations with an independent R calculation
across dense and sparse inputs, missing values, partial coverage, numerical
extremes, and serial or parallel execution. They verify the documented
calculation, not biological performance.

From the repository root, run the same build and check used by CI:

```sh
Rscript tools/check-package.R
Rscript tools/check-package.R --bioc
```

The second command also runs BiocCheck. Install the package dependencies first.
Logs remain under `.checks/`. The checks execute examples, tests, and vignettes.
The workflow covers Linux, macOS, and Windows, plus an R-devel Linux job.
Its configuration does not establish that every hosted platform has passed.

## Further reading

- [Scoring specification](inst/SCIENTIFIC_SPEC.md) - equation, value and
  coverage rules, properties, and limitations.
- [Vignette](vignettes/genefunnel.Rmd) - worked examples and real-data studies:
  six-gene scores in GSE125583 and mouse brain pseudobulk profiles in GSE60361.
  Both extracts ship with the package; the analyses are descriptive.
- [Tests](tests/testthat) - reference calculations and implementation checks.

For a rendered local vignette, install with `build_vignettes = TRUE` and its
suggested dependencies, then run `vignette("genefunnel")`. The first-run script
does not need a vignette build. Run `citation("genefunnel")` for citation details.

Regenerate the README figure from the example with
`Rscript tools/render-quickstart.R`. The package tests run the same example and
check its numerical results and exported files.

## License

GeneFunnel is licensed under GPL-3-or-later. See [COPYING](COPYING).
