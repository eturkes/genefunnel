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
#' @param mat A (sparse) matrix of feature expression (features x samples).
#' @param gene_sets A named list of character vectors representing gene sets.
#' @param BPPARAM A BiocParallel backend.
#'
#' @return A numeric matrix of gene set scores (gene sets x samples).
#' @export
#' @useDynLib genefunnel, .registration = TRUE
#' @importFrom Rcpp sourceCpp

genefunnel <- function(mat, gene_sets, BPPARAM = bpparam()) {
  if (!inherits(mat, "sparseMatrix")) {
    mat <- as.matrix(mat)
    storage.mode(mat) <- "numeric"
    mat <- Matrix(mat, sparse = TRUE)
  }

  if (any(mat < 0, na.rm = TRUE)) {
    stop("Input matrix contains negative values. genefunnel() expects all values to be non-negative.")
  }

  result <- bplapply(
    seq_len(ncol(mat)),
    function(i) {
      calculateScores(mat[, i, drop = FALSE], rownames(mat), gene_sets)
    },
    BPPARAM = BPPARAM
  )

  scores <- do.call(cbind, result)
  rownames(scores) <- names(gene_sets)
  colnames(scores) <- colnames(mat)

  return(scores)
}
