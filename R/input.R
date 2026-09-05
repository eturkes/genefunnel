.score_input_matrix <- function(mat, assay) {
    if (!isS4(mat) || !inherits(mat, "SummarizedExperiment")) {
        if (!is.null(assay)) {
            stop("`assay` must be NULL for matrix input.", call. = FALSE)
        }
        return(mat)
    }
    if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) {
        stop("Install SummarizedExperiment to use experiment input.",
            call. = FALSE)
    }

    assay_count <- length(SummarizedExperiment::assays(
        mat, withDimnames = FALSE
    ))
    if (assay_count == 0L) {
        stop("`mat` must contain at least one assay.", call. = FALSE)
    }
    if (is.null(assay)) {
        if (assay_count != 1L) {
            stop("`assay` must name one assay when `mat` has multiple assays.",
                call. = FALSE)
        }
        index <- 1L
    } else {
        valid <- is.character(assay) && !is.object(assay) &&
            is.null(dim(assay)) && length(assay) == 1L &&
            !is.na(assay) && nzchar(assay)
        if (!valid) {
            stop("`assay` must be NULL or one non-empty assay name.",
                call. = FALSE)
        }
        assay_names <- SummarizedExperiment::assayNames(mat)
        index <- which(!is.na(assay_names) & assay_names == assay)
        if (length(index) != 1L) {
            stop("`assay` must match exactly one assay name.", call. = FALSE)
        }
    }
    SummarizedExperiment::assay(mat, index)
}
