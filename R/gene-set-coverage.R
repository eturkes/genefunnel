# Assisted-by: OpenAI Codex.

.validate_gene_sets <- function(gene_sets) {
    valid_list <- is.list(gene_sets) &&
        !is.object(gene_sets) &&
        !is.data.frame(gene_sets) &&
        is.null(dim(gene_sets))
    if (!valid_list) {
        stop("`gene_sets` must be an unclassed named list.", call. = FALSE)
    }
    if (length(gene_sets) == 0L) {
        stop("`gene_sets` must contain at least one gene set.", call. = FALSE)
    }

    set_names <- names(gene_sets)
    if (is.null(set_names)) {
        stop("`gene_sets` must have names.", call. = FALSE)
    }
    .validate_identifier_vector(set_names, "`gene_sets` names")

    for (index in seq_along(gene_sets)) {
        members <- gene_sets[[index]]
        encoded_name <- encodeString(set_names[[index]], quote = '"')
        label <- sprintf("`gene_sets[[%s]]`", encoded_name)
        .validate_identifier_vector(
            members,
            label,
            allow_duplicates = TRUE
        )
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
            "Omitting %d gene set%s with fewer than two unique members ",
            count,
            if (count == 1L) "" else "s"
        ),
        "matched to `mat` row names: ",
        shown,
        suffix,
        ".",
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
#' @param gene_sets A non-empty unclassed named list of unclassed character
#'   vectors. Set names must be unique, non-missing, and non-empty. Member
#'   identifiers must be non-missing and non-empty.
#' @param features An unclassed character vector of unique, non-missing,
#'   non-empty feature identifiers, typically `rownames(mat)`.
#'
#' @details
#' This function uses the same validation, stable member deduplication, and
#' exact matching rules as [genefunnel()]. It performs no identifier mapping or
#' case conversion. An empty character vector is a valid declared set: its
#' coverage is undefined and it is not scoreable.
#'
#' @return A data frame with one row per gene set and columns `gene_set`,
#'   `declared_size` (unique declared members), `matched_size`,
#'   `unmatched_size`, `coverage`, `duplicate_member_count`, and `scoreable`.
#'   Coverage is `NA` for an empty declared set; a set is scoreable when at
#'   least two unique members match.
#'
#' @seealso [genefunnel()]
#' @export
#'
#' @examples
#' sets <- list(
#'     complete = c("A", "B"),
#'     partial = c("A", "B", "C"),
#'     duplicated = c("A", "A", "B"),
#'     insufficient = c("A", "Z")
#' )
#' coverage <- gene_set_coverage(sets, c("A", "B"))
#' coverage
#'
#' # Strict transcriptomic policy: retain complete sets.
#' strict <- !is.na(coverage$coverage) & coverage$coverage == 1
#' names(sets[strict])
#'
#' # Example low-coverage proteomic policy: >= 50% and >= 2 matches.
#' relaxed <- !is.na(coverage$coverage) &
#'     coverage$coverage >= 0.5 &
#'     coverage$matched_size >= 2
#' names(sets[relaxed])
gene_set_coverage <- function(gene_sets, features) {
    .prepare_gene_sets(gene_sets, features)$coverage
}
