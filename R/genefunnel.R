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
#'   be non-negative finite numbers, `NA`, or `NaN`. A
#'   [SummarizedExperiment::SummarizedExperiment] or subclass is also accepted
#'   when its selected assay has one of these matrix types.
#' @param gene_sets A non-empty named base list of character vectors. Set names
#'   and member identifiers must be non-missing and non-empty. Duplicate members
#'   within a set are used once, in their original order.
#' @param BPPARAM A [BiocParallel::BiocParallelParam] backend. The default is
#'   the currently registered [BiocParallel::bpparam()] backend.
#' @param min_coverage Minimum fraction of unique declared members observed in
#'   each sample, between zero and one. Unmatched members reduce this fraction.
#'   The default, zero, applies no fractional coverage requirement.
#' @param min_observed Minimum number of observed members per sample. This must
#'   be an integer of at least two. The default is two.
#' @param details Return scores and coverage diagnostics when `TRUE`.
#'   The default, `FALSE`, returns only the score matrix.
#' @param assay One exact assay name for `SummarizedExperiment` input.
#'   The default, `NULL`, selects the only assay. Multiple assays require a
#'   name. Leave this `NULL` for matrix input.
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
#' Each sample must also meet `min_observed` and `min_coverage`. Otherwise its
#' score is `NA`, but the row remains in the result. Coverage uses the unique
#' declared set size, not only the matched size. Zeros count as observations;
#' `NA`, `NaN`, and unmatched identifiers do not. These filters do not rescale
#' scores or correct bias from missing members.
#'
#' @section Experiment containers:
#' The selected assay uses the container's feature and sample names. Sparse
#' assays stay sparse. Other assays and metadata remain unchanged. Delayed or
#' disk-backed assays are not supported and are not silently materialized.
#' The return type is the same as for matrix input; metadata are not copied
#' into the score matrix or diagnostic list.
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
#'   With `details = TRUE`, a list contains `scores`, integer `observed` counts,
#'   per-sample `coverage` fractions, and `set_coverage`. The first three
#'   matrices share dimensions. `set_coverage` contains all declared sets,
#'   including skipped sets, as returned by [gene_set_coverage()].
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
    BPPARAM = BiocParallel::bpparam(),
    min_coverage = 0,
    min_observed = 2L,
    details = FALSE,
    assay = NULL
) {
    .validate_observation_policy(min_coverage, min_observed, details)
    mat <- .score_input_matrix(mat, assay)
    .validate_score_matrix(mat)
    prepared <- .prepare_gene_sets(gene_sets, rownames(mat))
    .validate_bpparam(BPPARAM)

    retained <- prepared$coverage$scoreable
    if (any(!retained)) {
        .warn_unscoreable_sets(prepared$coverage$gene_set[!retained])
    }

    if (!any(retained)) {
        return(.assemble_score_output(NULL, mat, prepared$coverage, details))
    }

    indices <- prepared$indices[retained]
    declared_size <- prepared$coverage$declared_size[retained]
    minimum_observed <- pmax(
        min_observed,
        floor(min_coverage * declared_size)
    )
    # Compare fractions directly: ceiling(0.14 * 50) can incorrectly require 8.
    minimum_observed <- as.integer(minimum_observed +
        (minimum_observed / declared_size < min_coverage))
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
        minimum_observed = minimum_observed,
        return_observed = details,
        BPPARAM = BPPARAM
    )

    .assemble_score_output(chunks, mat, prepared$coverage, details)
}

.assemble_score_output <- function(chunks, mat, set_coverage, details) {
    retained_names <- set_coverage$gene_set[set_coverage$scoreable]
    assemble <- function(field) {
        if (length(retained_names) == 0L) {
            return(matrix(double(), 0L, ncol(mat)))
        }
        .assemble_score_chunks(
            chunks, length(retained_names), ncol(mat), field = field
        )
    }
    scores <- assemble("scores")
    dimnames(scores) <- list(retained_names, colnames(mat))
    observed <- NULL
    if (details) {
        observed <- assemble("observed")
        storage.mode(observed) <- "integer"
        dimnames(observed) <- dimnames(scores)
    }
    .score_output(scores, observed, set_coverage, details)
}

.score_output <- function(scores, observed, set_coverage, details) {
    if (!details) {
        return(scores)
    }
    declared_size <- set_coverage$declared_size[set_coverage$scoreable]
    list(
        scores = scores,
        observed = observed,
        coverage = observed / declared_size,
        set_coverage = set_coverage
    )
}

.score_matrix_chunk <- function(
    mat, gene_indices, storage, minimum_observed = NULL,
    return_observed = FALSE
) {
    switch(
        storage,
        dense = calculateScoresDense(
            .as_dense_numeric_chunk(mat),
            gene_indices,
            minimum_observed,
            return_observed
        ),
        sparse = calculateScoresSparse(
            .as_sparse_numeric_chunk(mat),
            gene_indices,
            minimum_observed,
            return_observed
        ),
        stop("Internal matrix storage kind is invalid.", call. = FALSE)
    )
}
