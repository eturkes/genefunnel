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
#' Calculate one score for each gene set and sample in a non-negative
#' feature-by-sample matrix. Unlike a plain sum, GeneFunnel gives lower scores
#' when the observed members of a set disagree.
#'
#' @param mat A base numeric or integer matrix without an extra class, or a
#'   numeric dense or sparse [Matrix::Matrix] object. Features must be in rows
#'   and samples in columns. Row names must be unique and non-empty. Values can
#'   be non-negative finite numbers, `NA`, or `NaN`.
#' @param gene_sets A non-empty named base list of character vectors. Set names
#'   and member identifiers must be non-missing and non-empty. Duplicate members
#'   within a set are used once, in their original order.
#' @param BPPARAM A [BiocParallel::BiocParallelParam] backend. The default is
#'   the currently registered [BiocParallel::bpparam()] backend.
#'
#' @details
#' Gene-set members are matched to matrix row names. `NA` and `NaN` values are
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
#' Both explicit and implicit sparse zeros count as observed values. A score is
#' `NA_real_` when fewer than two members remain after missing values are
#' removed.
#'
#' @section Matching and coverage:
#' Matching is exact and case-sensitive. GeneFunnel does not convert identifiers
#' or query annotation resources. Unmatched members are ignored. A set is
#' scored when at least two unique members match. Use [gene_set_coverage()] to
#' inspect coverage and choose an appropriate threshold for your experiment.
#'
#' @section Interpretation:
#' The score is high when a set has both a large total signal and similar member
#' values. It therefore depends on measurement units, preprocessing, set size,
#' and feature coverage. Compare scores across data sets only when these inputs
#' are compatible. GeneFunnel does not normalize, impute, test, or interpret the
#' data.
#'
#' @return A base numeric matrix with scored gene sets in rows and samples in
#'   columns. Gene-set and sample order are preserved. If no set has at least
#'   two matched members, the result has zero rows and `ncol(mat)` columns.
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
    storage <- if (inherits(mat, "sparseMatrix")) "sparse" else "dense"
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
