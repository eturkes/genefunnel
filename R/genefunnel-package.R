#' GeneFunnel gene-set activity scoring
#'
#' @description
#' GeneFunnel scores predefined gene sets in each sample of a non-negative
#' feature-by-sample matrix, combining total signal with member agreement.
#'
#' @details
#' Use [gene_set_coverage()] to inspect matching and [genefunnel()] to score
#' dense or sparse matrices. Zeros remain observed; `NA` and `NaN` are omitted
#' per sample. Set `min_coverage` and `min_observed` to enforce observation
#' requirements, and `details = TRUE` to inspect coverage.
#'
#' GeneFunnel does not normalize, impute, map identifiers, retrieve sets,
#' test hypotheses, or interpret results biologically.
#'
#' See `vignette("genefunnel")` for a worked workflow and
#' `system.file("SCIENTIFIC_SPEC.md", package = "genefunnel")` for the full
#' scoring specification.
#'
#' @seealso [genefunnel()], [gene_set_coverage()], `citation("genefunnel")`
#' @aliases genefunnel-package NULL
#' @keywords internal
#' @useDynLib genefunnel, .registration = TRUE
#' @importClassesFrom Matrix dMatrix sparseMatrix
#' @importFrom Rcpp evalCpp
"_PACKAGE"
