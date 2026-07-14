#    This file is part of GeneFunnel.
#    Copyright (C) 2025  Emir Turkes, UK DRI at UCL
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
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
#    Emir Turkes can be contacted at emir.turkes@eturkes.com

#' GeneFunnel scoring algorithm
#'
#' This function computes GeneFunnel scores for each sample and gene set, based
#' on a matrix of feature expression and a list of gene sets.
#'
#' @param mat A base numeric/integer matrix or numeric `Matrix` object, with
#'   features in rows and samples in columns. Row names must be unique.
#' @param gene_sets A named list of character vectors containing feature
#'   identifiers.
#' @param BPPARAM A `BiocParallelParam` backend.
#'
#' @return A numeric matrix of gene set scores (gene sets x samples).
#' @export
#' @useDynLib genefunnel, .registration = TRUE
#' @importFrom Rcpp sourceCpp

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

    mat <- .as_sparse_numeric_matrix(mat)
    indices <- prepared$indices[retained]
    result <- BiocParallel::bplapply(
        seq_len(ncol(mat)),
        function(column) {
            calculateScores(mat[, column, drop = FALSE], indices)
        },
        BPPARAM = BPPARAM
    )

    scores <- do.call(cbind, result)
    rownames(scores) <- retained_names
    colnames(scores) <- colnames(mat)
    scores
}
