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

    contains_negative <- FALSE
    for (column in seq_len(ncol(mat))) {
        column_mat <- mat[, column, drop = FALSE]
        values <- if (is.matrix(column_mat)) {
            column_mat
        } else {
            methods::slot(column_mat, "x")
        }
        if (.contains_negative_score_value(values)) {
            contains_negative <- TRUE
        }
    }

    if (contains_negative) {
        stop(
            paste0(
                "`mat` contains negative values; all observed values must ",
                "be non-negative."
            ),
            call. = FALSE
        )
    }

    invisible(mat)
}

.contains_negative_score_value <- function(values) {
    block_size <- 65536L
    value_count <- length(values)
    first <- 1
    contains_negative <- FALSE

    while (first <= value_count) {
        last <- min(first + block_size - 1, value_count)
        block <- values[seq.int(first, last)]
        if (any(is.infinite(block))) {
            stop(
                paste0(
                    "`mat` contains infinite values; only finite values, ",
                    "NA, and NaN are allowed."
                ),
                call. = FALSE
            )
        }
        if (!contains_negative && any(block < 0, na.rm = TRUE)) {
            contains_negative <- TRUE
        }
        first <- last + 1
    }

    contains_negative
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

.matrix_storage <- function(mat) {
    if (inherits(mat, "sparseMatrix")) "sparse" else "dense"
}

.as_dense_numeric_chunk <- function(mat) {
    if (inherits(mat, "sparseMatrix")) {
        stop("Internal dense scoring received a sparse matrix.", call. = FALSE)
    }
    if (is.matrix(mat) && typeof(mat) == "double") {
        return(mat)
    }
    if (is.matrix(mat)) {
        return(matrix(
            as.double(mat),
            nrow = nrow(mat),
            ncol = ncol(mat),
            dimnames = dimnames(mat)
        ))
    }
    as.matrix(mat)
}

.as_sparse_numeric_chunk <- function(mat) {
    if (!inherits(mat, "sparseMatrix")) {
        stop("Internal sparse scoring received a dense matrix.", call. = FALSE)
    }
    mat <- methods::as(mat, "generalMatrix")
    methods::as(mat, "CsparseMatrix")
}
