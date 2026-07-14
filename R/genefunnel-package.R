#' GeneFunnel gene-set activity scoring
#'
#' @description
#' GeneFunnel converts a non-negative feature-by-sample matrix into a numeric
#' gene-set-by-sample score matrix. Every gene set and sample is scored
#' independently from the available member values.
#'
#' @details
#' [genefunnel()] is the primary scorer. [gene_set_coverage()] reports exact
#' identifier coverage so callers can apply an experiment-specific coverage
#' policy before scoring. Dense and sparse inputs are supported; zeros are
#' observed values, `NA` and `NaN` are omitted per sample, and partially
#' covered sets remain scoreable when at least two unique members match.
#'
#' GeneFunnel does not normalize or impute measurements, map identifiers,
#' retrieve gene sets, perform downstream tests, or supply biological
#' interpretation.
#'
#' See `vignette("genefunnel")` for a worked workflow and
#' `system.file("SCIENTIFIC_SPEC.md", package = "genefunnel")` for the durable
#' mathematical and behavioral contract.
#'
#' @seealso [genefunnel()], [gene_set_coverage()], `citation("genefunnel")`
#' @aliases genefunnel-package NULL
#' @keywords internal
#' @useDynLib genefunnel, .registration = TRUE
#' @importClassesFrom Matrix dMatrix sparseMatrix
#' @importFrom Rcpp evalCpp
"_PACKAGE"
