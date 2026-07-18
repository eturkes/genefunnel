# Assisted-by: OpenAI Codex.

.catalogue_stop <- function(detail) {
    stop("`catalogue` is malformed: ", detail, ".", call. = FALSE)
}

.catalogue_coverage <- function(members, member_index, duplicate_member_count) {
    declared_size <- unname(lengths(members))
    matched_size <- unname(vapply(
        member_index,
        function(index) sum(index > 0L),
        integer(1)
    ))
    coverage_fraction <- rep.int(NA_real_, length(members))
    nonempty <- declared_size > 0L
    coverage_fraction[nonempty] <-
        matched_size[nonempty] / declared_size[nonempty]
    data.frame(
        gene_set = unname(names(members)),
        declared_size = declared_size,
        matched_size = matched_size,
        unmatched_size = declared_size - matched_size,
        coverage = coverage_fraction,
        duplicate_member_count = duplicate_member_count,
        scoreable = matched_size >= 2L,
        row.names = seq_along(members),
        check.names = FALSE
    )
}

.catalogue_adjacency <- function(indices, feature_count) {
    matched_count <- sum(lengths(indices))
    if (matched_count > .Machine$integer.max) {
        stop(
            "Compiled catalogue has too many matched memberships.",
            call. = FALSE
        )
    }
    rows <- unlist(indices, use.names = FALSE)
    sets <- rep.int(seq_along(indices), lengths(indices))
    counts <- tabulate(rows, nbins = feature_count)
    offsets <- as.integer(c(0, cumsum(counts)))
    if (matched_count == 0L) {
        set_indices <- integer()
    } else {
        ordering <- order(rows, sets, method = "radix")
        set_indices <- as.integer(sets[ordering])
    }
    list(offsets = offsets, set_indices = set_indices)
}

.build_compiled_catalogue <- function(
    members,
    duplicate_member_count,
    features
) {
    .validate_gene_sets(members)
    .validate_features(features)
    if (length(features) > .Machine$integer.max) {
        stop("Compiled catalogue has too many features.", call. = FALSE)
    }
    if (!is.integer(duplicate_member_count) ||
        !is.null(attributes(duplicate_member_count)) ||
        length(duplicate_member_count) != length(members) ||
        anyNA(duplicate_member_count) || any(duplicate_member_count < 0L)) {
        stop("Compiled catalogue duplicate counts are invalid.", call. = FALSE)
    }
    if (any(vapply(members, anyDuplicated, integer(1)) > 0L)) {
        stop(
            "Compiled catalogue members must already be canonical.",
            call. = FALSE
        )
    }
    member_index <- lapply(members, function(set) {
        unname(match(set, features, nomatch = 0L))
    })
    names(member_index) <- names(members)
    indices <- lapply(member_index, function(index) {
        unname(index[index > 0L])
    })
    names(indices) <- names(members)
    coverage <- .catalogue_coverage(
        members,
        member_index,
        duplicate_member_count
    )
    row_adjacency <- .catalogue_adjacency(indices, length(features))
    state <- list(
        schema_version = .GF_CATALOGUE_SCHEMA_VERSION,
        formula_version = .GF_FORMULA_VERSION,
        features = features,
        members = members,
        duplicate_member_count = duplicate_member_count,
        member_index = member_index,
        indices = indices,
        coverage = coverage,
        row_adjacency = row_adjacency
    )
    state$fingerprints <- .catalogue_fingerprints(state)
    structure(state, class = "genefunnel_catalogue")
}

.compile_gene_sets <- function(gene_sets, features) {
    .validate_gene_sets(gene_sets)
    .validate_features(features)
    members <- lapply(gene_sets, function(set) {
        unname(set[!duplicated(set)])
    })
    names(members) <- unname(names(gene_sets))
    duplicate_member_count <- as.integer(lengths(gene_sets) - lengths(members))
    .build_compiled_catalogue(
        members,
        unname(duplicate_member_count),
        unname(features)
    )
}

.catalogue_integer_vector <- function(value, length, minimum, maximum) {
    is.integer(value) &&
        !is.object(value) &&
        is.null(attributes(value)) &&
        base::length(value) == length &&
        !anyNA(value) &&
        all(value >= minimum) &&
        all(value <= maximum)
}

.catalogue_named_list <- function(value, expected_names) {
    is.list(value) &&
        !is.object(value) &&
        identical(names(value), expected_names)
}

.catalogue_validate_identifiers <- function(state) {
    identifier_error <- tryCatch({
        .validate_features(state$features)
        .validate_gene_sets(state$members)
        NULL
    }, error = identity)
    if (!is.null(identifier_error)) {
        .catalogue_stop(conditionMessage(identifier_error))
    }
    if (!is.null(names(state$features)) ||
        any(vapply(
            state$members,
            function(set) !is.null(names(set)),
            logical(1)
        )) ||
        any(vapply(state$members, anyDuplicated, integer(1)) > 0L)) {
        .catalogue_stop(
            "canonical identifiers have unexpected names or duplicates"
        )
    }
    invisible(NULL)
}

.catalogue_validate_memberships <- function(state) {
    set_names <- names(state$members)
    set_count <- length(state$members)
    feature_count <- length(state$features)
    valid_duplicates <- .catalogue_integer_vector(
        state$duplicate_member_count,
        set_count,
        0L,
        .Machine$integer.max
    )
    if (!valid_duplicates) {
        .catalogue_stop("duplicate-member counts are invalid")
    }
    if (!.catalogue_named_list(state$member_index, set_names) ||
        !.catalogue_named_list(state$indices, set_names)) {
        .catalogue_stop("membership lists are invalid")
    }

    expected_indices <- vector("list", set_count)
    names(expected_indices) <- set_names
    for (index in seq_len(set_count)) {
        members <- state$members[[index]]
        mapping <- state$member_index[[index]]
        valid_mapping <- .catalogue_integer_vector(
            mapping,
            length(members),
            0L,
            feature_count
        )
        if (!valid_mapping) {
            .catalogue_stop(sprintf("member mapping %d is invalid", index))
        }
        matched <- mapping > 0L
        if (any(matched) &&
            !all(members[matched] == state$features[mapping[matched]])) {
            .catalogue_stop(sprintf(
                "member mapping %d disagrees with identifiers",
                index
            ))
        }
        expected_indices[[index]] <- unname(mapping[matched])
        if (!identical(state$indices[[index]], expected_indices[[index]])) {
            .catalogue_stop(sprintf(
                "matched indices %d are not derived from mappings",
                index
            ))
        }
    }
    expected_indices
}

.catalogue_validate_derived <- function(state, expected_indices) {
    expected_coverage <- .catalogue_coverage(
        state$members,
        state$member_index,
        state$duplicate_member_count
    )
    if (!identical(state$coverage, expected_coverage)) {
        .catalogue_stop("coverage is not derived from canonical state")
    }
    expected_adjacency <- .catalogue_adjacency(
        expected_indices,
        length(state$features)
    )
    valid_adjacency <- is.list(state$row_adjacency) &&
        !is.object(state$row_adjacency) &&
        identical(names(state$row_adjacency), c("offsets", "set_indices")) &&
        identical(state$row_adjacency, expected_adjacency)
    if (!valid_adjacency) {
        .catalogue_stop("row-to-set adjacency is invalid")
    }
    invisible(NULL)
}

.catalogue_validate_fingerprints <- function(state) {
    expected_fingerprint_names <- c(
        "algorithm", "encoding", "catalogue", "features", "content"
    )
    valid_fingerprints <- is.character(state$fingerprints) &&
        !is.object(state$fingerprints) &&
        identical(names(state$fingerprints), expected_fingerprint_names) &&
        identical(length(state$fingerprints), 5L) &&
        identical(
            unname(state$fingerprints[seq_len(2L)]),
            c("SHA-256", "GFCAT-1")
        ) &&
        all(grepl(
            "^[0-9a-f]{64}$",
            unname(state$fingerprints[seq.int(3L, 5L)])
        ))
    if (!valid_fingerprints ||
        !identical(state$fingerprints, .catalogue_fingerprints(state))) {
        .catalogue_stop("fingerprints do not identify the canonical state")
    }
    invisible(NULL)
}

.catalogue_state <- function(catalogue, features = NULL) {
    expected_fields <- c(
        "schema_version", "formula_version", "features", "members",
        "duplicate_member_count", "member_index", "indices", "coverage",
        "row_adjacency", "fingerprints"
    )
    valid_object <- is.list(catalogue) &&
        !isS4(catalogue) &&
        identical(class(catalogue), "genefunnel_catalogue") &&
        identical(names(attributes(catalogue)), c("names", "class"))
    if (!valid_object) {
        .catalogue_stop("it must be an intact genefunnel_catalogue object")
    }
    state <- unclass(catalogue)
    if (!identical(names(state), expected_fields)) {
        .catalogue_stop("field names or order are invalid")
    }
    if (!identical(state$schema_version, .GF_CATALOGUE_SCHEMA_VERSION)) {
        .catalogue_stop("schema version is unsupported")
    }
    if (!identical(state$formula_version, .GF_FORMULA_VERSION)) {
        .catalogue_stop("formula version is unsupported")
    }

    .catalogue_validate_identifiers(state)
    expected_indices <- .catalogue_validate_memberships(state)
    .catalogue_validate_derived(state, expected_indices)
    .catalogue_validate_fingerprints(state)

    if (!is.null(features)) {
        .validate_features(features)
        if (!identical(unname(features), state$features)) {
            stop(
                "`catalogue` feature identities and order do not match ",
                "`mat` row names.",
                call. = FALSE
            )
        }
    }
    state
}

.validate_compiled_catalogue <- function(catalogue, features = NULL) {
    invisible(.catalogue_state(catalogue, features))
}

.genefunnel_compiled <- function(
    mat,
    catalogue,
    BPPARAM = BiocParallel::bpparam()
) {
    .validate_score_matrix(mat)
    state <- .catalogue_state(catalogue, rownames(mat))
    .validate_bpparam(BPPARAM)

    retained <- state$coverage$scoreable
    if (any(!retained)) {
        .warn_unscoreable_sets(state$coverage$gene_set[!retained])
    }
    retained_names <- state$coverage$gene_set[retained]
    if (!any(retained)) {
        return(matrix(
            double(),
            nrow = 0L,
            ncol = ncol(mat),
            dimnames = list(character(), colnames(mat))
        ))
    }

    indices <- state$indices[retained]
    storage <- .matrix_storage(mat)
    ranges <- .column_chunk_ranges(
        ncol(mat),
        BiocParallel::bpnworkers(BPPARAM)
    )
    chunks <- BiocParallel::bpiterate(
        .matrix_chunk_iterator(mat, ranges),
        .score_matrix_task,
        gene_indices = indices,
        storage = storage,
        BPPARAM = BPPARAM
    )
    scores <- .assemble_score_chunks(
        chunks,
        n_gene_sets = length(indices),
        n_columns = ncol(mat)
    )
    rownames(scores) <- retained_names
    colnames(scores) <- colnames(mat)
    scores
}
