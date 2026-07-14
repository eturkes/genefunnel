.validate_gene_sets <- function(gene_sets) {
    valid_list <- is.list(gene_sets) &&
        !is.data.frame(gene_sets) &&
        is.null(dim(gene_sets))
    if (!valid_list) {
        stop("`gene_sets` must be a named list.", call. = FALSE)
    }
    if (length(gene_sets) == 0L) {
        stop("`gene_sets` must contain at least one gene set.", call. = FALSE)
    }

    set_names <- names(gene_sets)
    if (is.null(set_names)) {
        stop("`gene_sets` must have names.", call. = FALSE)
    }
    if (anyNA(set_names)) {
        stop("`gene_sets` contains a missing name.", call. = FALSE)
    }
    if (any(!nzchar(set_names))) {
        stop("`gene_sets` contains an empty name.", call. = FALSE)
    }
    if (anyDuplicated(set_names)) {
        stop("`gene_sets` contains duplicated names.", call. = FALSE)
    }

    for (index in seq_along(gene_sets)) {
        members <- gene_sets[[index]]
        encoded_name <- encodeString(set_names[[index]], quote = '"')
        label <- sprintf("`gene_sets[[%s]]`", encoded_name)
        if (!is.character(members) || !is.null(dim(members))) {
            stop(label, " must be a character vector.", call. = FALSE)
        }
        if (anyNA(members)) {
            stop(label, " contains missing member identifiers.", call. = FALSE)
        }
        if (any(!nzchar(members))) {
            stop(label, " contains empty member identifiers.", call. = FALSE)
        }
    }

    invisible(gene_sets)
}

.prepare_gene_sets <- function(gene_sets, features) {
    .validate_gene_sets(gene_sets)
    .validate_features(features)

    members <- lapply(gene_sets, function(set) set[!duplicated(set)])
    declared_size <- unname(lengths(members))
    flat_indices <- match(
        unlist(members, use.names = FALSE),
        features,
        nomatch = 0L
    )
    ends <- cumsum(declared_size)
    starts <- ends - declared_size + 1L
    indices <- Map(function(first, last, size) {
        if (size == 0L) {
            return(integer())
        }
        matched <- flat_indices[seq.int(first, last)]
        unname(matched[matched > 0L])
    }, starts, ends, declared_size)
    names(indices) <- names(members)

    matched_size <- unname(lengths(indices))
    duplicate_member_count <- unname(lengths(gene_sets)) - declared_size
    coverage_fraction <- rep.int(NA_real_, length(gene_sets))
    nonempty <- declared_size > 0L
    coverage_fraction[nonempty] <-
        matched_size[nonempty] / declared_size[nonempty]

    coverage <- data.frame(
        gene_set = names(gene_sets),
        declared_size = declared_size,
        matched_size = matched_size,
        unmatched_size = declared_size - matched_size,
        coverage = coverage_fraction,
        duplicate_member_count = duplicate_member_count,
        scoreable = matched_size >= 2L,
        row.names = seq_along(gene_sets),
        check.names = FALSE
    )

    list(members = members, indices = indices, coverage = coverage)
}

.warn_unscoreable_sets <- function(set_names) {
    count <- length(set_names)
    shown_names <- utils::head(set_names, 5L)
    remainder <- count - length(shown_names)
    shown <- paste(encodeString(shown_names, quote = '"'), collapse = ", ")
    suffix <- if (remainder > 0L) {
        sprintf(", and %d more", remainder)
    } else {
        ""
    }
    warning(
        sprintf(
            paste0(
                "Omitting %d gene set%s with fewer than two unique members ",
                "matched to `mat` row names: %s%s."
            ),
            count,
            if (count == 1L) "" else "s",
            shown,
            suffix
        ),
        call. = FALSE
    )
}

#' Gene-set coverage diagnostics
#'
#' Reports exact, case-sensitive overlap between declared gene-set members and
#' available features. Duplicate members are counted, then removed while
#' preserving their first occurrence. The function reports facts only; callers
#' choose a coverage policy appropriate to their experiment.
#'
#' @param gene_sets A named list of character vectors containing feature
#'   identifiers.
#' @param features A character vector of unique available feature identifiers.
#'
#' @return A data frame with one row per gene set and columns `gene_set`,
#'   `declared_size` (unique declared members), `matched_size`,
#'   `unmatched_size`, `coverage`, `duplicate_member_count`, and `scoreable`.
#'   Coverage is `NA` for an empty declared set; a set is scoreable when at
#'   least two unique members match.
#' @export
#'
#' @examples
#' sets <- list(pathway = c("A", "B", "C"))
#' gene_set_coverage(sets, c("A", "B"))
gene_set_coverage <- function(gene_sets, features) {
    .prepare_gene_sets(gene_sets, features)$coverage
}
