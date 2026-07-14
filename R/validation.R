.validate_identifier_vector <- function(identifiers, label) {
    if (!is.character(identifiers) || !is.null(dim(identifiers))) {
        stop(label, " must be a character vector.", call. = FALSE)
    }
    if (anyNA(identifiers)) {
        stop(label, " contains missing identifiers.", call. = FALSE)
    }
    if (any(!nzchar(identifiers))) {
        stop(label, " contains empty identifiers.", call. = FALSE)
    }
    if (anyDuplicated(identifiers)) {
        stop(label, " contains duplicated identifiers.", call. = FALSE)
    }
    invisible(identifiers)
}

.validate_score_matrix <- function(mat) {
    base_numeric_matrix <- is.matrix(mat) &&
        typeof(mat) %in% c("double", "integer")
    matrix_numeric_matrix <- inherits(mat, "dMatrix")
    if (!base_numeric_matrix && !matrix_numeric_matrix) {
        stop(
            "`mat` must be a base numeric or integer matrix, or a numeric Matrix object.",
            call. = FALSE
        )
    }

    if (nrow(mat) == 0L) {
        stop("`mat` must have at least one row.", call. = FALSE)
    }
    if (ncol(mat) == 0L) {
        stop("`mat` must have at least one column.", call. = FALSE)
    }

    features <- rownames(mat)
    if (is.null(features)) {
        stop("`mat` must have row names.", call. = FALSE)
    }
    .validate_identifier_vector(features, "`mat` row names")

    if (any(is.infinite(mat))) {
        stop(
            "`mat` contains infinite values; only finite values, NA, and NaN are allowed.",
            call. = FALSE
        )
    }
    if (any(mat < 0, na.rm = TRUE)) {
        stop(
            "`mat` contains negative values; all observed values must be non-negative.",
            call. = FALSE
        )
    }

    invisible(mat)
}

.validate_features <- function(features) {
    .validate_identifier_vector(features, "`features`")
}

.validate_bpparam <- function(BPPARAM) {
    if (!methods::is(BPPARAM, "BiocParallelParam")) {
        stop("`BPPARAM` must be a BiocParallelParam object.", call. = FALSE)
    }
    invisible(BPPARAM)
}

.as_sparse_numeric_matrix <- function(mat) {
    if (is.matrix(mat)) {
        mat <- Matrix::Matrix(mat, sparse = TRUE, doDiag = FALSE)
    }
    mat <- methods::as(mat, "generalMatrix")
    methods::as(mat, "CsparseMatrix")
}
