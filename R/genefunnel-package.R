#' GeneFunnel gene-set activity scoring
#'
#' @description
#' GeneFunnel scores predefined gene sets in each sample of a non-negative
#' feature-by-sample matrix. The score reflects both the total observed signal
#' and the agreement among members of the set.
#'
#' @details
#' Use [genefunnel()] to calculate scores. Use [gene_set_coverage()] first when
#' you need to check how many members of each set are present. Dense and sparse
#' matrices are supported. Zeros are treated as observations, while `NA` and
#' `NaN` are omitted separately in each sample.
#' Set `min_coverage` and `min_observed` to enforce sample-level observation
#' requirements, and `details = TRUE` to inspect the effective coverage.
#'
#' GeneFunnel only calculates scores. It does not normalize or impute data, map
#' identifiers, retrieve gene sets, run statistical tests, or interpret the
#' results biologically.
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
