# Assisted-by: OpenAI Codex.

#    This file is part of GeneFunnel.
#    Copyright (C) 2025-2026  Emir Turkes, UK DRI at UCL
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <https://www.gnu.org/licenses/>.
#
#    Emir Turkes can be contacted at emir.turkes@eturkes.com

#' Calculate GeneFunnel scores
#'
#' Convert a non-negative feature-by-sample matrix into sample-wise gene-set
#' scores. Every retained gene set and sample is calculated independently.
#'
#' @param mat A base numeric or integer matrix, or a numeric dense/sparse
#'   [Matrix::Matrix] object, with features in rows and samples in columns. It
#'   must have at least one row and column. Row names must be unique,
#'   non-missing, and non-empty. Values may be finite and non-negative, `NA`,
#'   or `NaN`; negative and infinite values are rejected.
#' @param gene_sets A non-empty named list of character vectors. Set names must
#'   be unique, non-missing, and non-empty. Member identifiers must be
#'   non-missing and non-empty; duplicate members are deduplicated in first
#'   occurrence order.
#' @param BPPARAM A [BiocParallel::BiocParallelParam] backend. The default is
#'   the currently registered [BiocParallel::bpparam()] backend.
#'
#' @details
#' For the members of one gene set that match matrix rows, `NA` and `NaN` are
#' removed separately in each sample. Given the remaining values
#' \eqn{x_1,\ldots,x_n}, where \eqn{n \ge 2}, the score is
#'
#' \deqn{
#' \sum_{i=1}^{n}x_i
#' -\frac{n}{2(n-1)}\sum_{i=1}^{n}|x_i-\bar{x}|,
#' \qquad
#' \bar{x}=\frac{1}{n}\sum_{i=1}^{n}x_i.
#' }
#'
#' Zeros are observed values: explicit dense zeros and implicit sparse zeros
#' contribute to the effective size, sum, mean, and deviation. Missing values
#' are absent information and reduce the sample-specific effective size. A
#' cell with fewer than two observed members is returned as `NA_real_`.
#'
#' @section Matching and coverage:
#' Matching against `rownames(mat)` is exact and case-sensitive. The function
#' performs no identifier conversion or annotation lookup. Unmatched members
#' are ignored. Sets with at least two unique matched members are scored; all
#' other sets are omitted with at most one aggregate warning. Use
#' [gene_set_coverage()] to inspect coverage and apply an experiment-specific
#' threshold before scoring.
#'
#' @section Interpretation:
#' The deviation penalty treats within-set absolute variability as evidence
#' against coherent activity. The sum term retains the effects of units,
#' library size, preprocessing, gene-set size, and assay coverage. GeneFunnel
#' performs no normalization, imputation, identifier mapping, downstream
#' testing, or biological interpretation. Scores from unrelated data sets are
#' comparable only when units, preprocessing, identifiers, gene-set versions,
#' and coverage are compatible.
#'
#' @return A base numeric matrix with retained gene sets in rows and samples in
#'   columns. Retained set order, set names, sample order, and existing column
#'   names are preserved. When every set is omitted, the result has zero rows
#'   and `ncol(mat)` columns.
#'
#' @references
#' Turkes E (2025). *Development of Gene Set Enrichment and Imputation Methods
#' for Transcriptomics and Proteomics: Application in the Study of
#' Neurofibrillary Tangle-bearing Neurons in Alzheimer's Disease*. PhD thesis,
#' University College London.
#'
#' @seealso [gene_set_coverage()], `vignette("genefunnel")`,
#'   `system.file("SCIENTIFIC_SPEC.md", package = "genefunnel")`
#' @export
#'
#' @examples
#' mat <- rbind(
#'     A = c(equal = 4, zeros = 0, one_nonzero = 4, missing = 1),
#'     B = c(equal = 4, zeros = 0, one_nonzero = 0, missing = 2),
#'     C = c(equal = 4, zeros = 0, one_nonzero = 0, missing = NA)
#' )
#' gene_sets <- list(
#'     three_members = c("A", "B", "C"),
#'     two_members = c("A", "B")
#' )
#'
#' genefunnel(
#'     mat,
#'     gene_sets,
#'     BPPARAM = BiocParallel::SerialParam()
#' )

genefunnel <- function(
    mat,
    gene_sets,
    BPPARAM = BiocParallel::bpparam()
) {
    .validate_score_matrix(mat)
    prepared <- .prepare_gene_sets(gene_sets, rownames(mat))
    .validate_bpparam(BPPARAM)

    retained <- prepared$coverage$scoreable
    if (any(!retained)) {
        .warn_unscoreable_sets(prepared$coverage$gene_set[!retained])
    }

    retained_names <- prepared$coverage$gene_set[retained]
    if (!any(retained)) {
        return(matrix(
            double(),
            nrow = 0L,
            ncol = ncol(mat),
            dimnames = list(character(), colnames(mat))
        ))
    }

    indices <- prepared$indices[retained]
    storage <- .matrix_storage(mat)
    ranges <- .column_chunk_ranges(
        ncol(mat),
        BiocParallel::bpnworkers(BPPARAM)
    )
    # ITER slices on the manager; the worker function never captures `mat`.
    chunks <- BiocParallel::bpiterate(
        .matrix_chunk_iterator(mat, ranges),
        .score_matrix_task,
        gene_indices = indices,
        storage = storage,
        BPPARAM = BPPARAM
    )

    scores <- .assemble_score_chunks(
        chunks,
        n_gene_sets = length(indices),
        n_columns = ncol(mat)
    )
    rownames(scores) <- retained_names
    colnames(scores) <- colnames(mat)
    scores
}

.score_matrix_chunk <- function(mat, gene_indices, storage) {
    switch(
        storage,
        dense = calculateScoresDense(
            .as_dense_numeric_chunk(mat),
            gene_indices
        ),
        sparse = calculateScoresSparse(
            .as_sparse_numeric_chunk(mat),
            gene_indices
        ),
        stop("Internal matrix storage kind is invalid.", call. = FALSE)
    )
}
