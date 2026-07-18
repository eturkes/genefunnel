# Assisted-by: OpenAI Codex.

.validate_aggregation_matrix_structure <- function(mat) {
    base_numeric_matrix <- is.matrix(mat) &&
        !is.object(mat) &&
        typeof(mat) %in% c("double", "integer")
    matrix_numeric_matrix <- isS4(mat) && inherits(mat, "dMatrix")
    if (!base_numeric_matrix && !matrix_numeric_matrix) {
        stop(
            "`mat` must be an unclassed base numeric or integer matrix, ",
            "or a numeric Matrix object.",
            call. = FALSE
        )
    }
    if (nrow(mat) == 0L) {
        stop("`mat` must have at least one row.", call. = FALSE)
    }
    if (ncol(mat) == 0L) {
        stop("`mat` must have at least one column.", call. = FALSE)
    }
    if (is.null(rownames(mat))) {
        stop("`mat` must have row names.", call. = FALSE)
    }
    .validate_identifier_vector(rownames(mat), "`mat` row names")
    if (is.null(colnames(mat))) {
        stop(
            "Aggregation `mat` must have unique column names for unit IDs.",
            call. = FALSE
        )
    }
    .validate_identifier_vector(colnames(mat), "`mat` column names")
    invisible(mat)
}

.validate_aggregation_named_vector <- function(value, units, type, label) {
    valid_type <- switch(
        type,
        character = is.character(value),
        numeric = typeof(value) %in% c("double", "integer"),
        FALSE
    )
    valid_attributes <- identical(names(attributes(value)), "names")
    if (!valid_type || is.object(value) || !is.null(dim(value)) ||
        length(value) != length(units) || !valid_attributes ||
        !identical(names(value), units)) {
        stop(
            label,
            " must be an unclassed vector named exactly and in order by ",
            "`mat` column names.",
            call. = FALSE
        )
    }
    invisible(value)
}

.aggregation_effective_weights <- function(groups, weights) {
    group_order <- unique(unname(groups))
    input_weight <- as.double(weights)
    effective_weight <- rep.int(NA_real_, length(weights))
    for (group in group_order) {
        selected <- groups == group
        maximum <- max(input_weight[selected])
        if (maximum == 0) {
            next
        }
        scaled <- input_weight[selected] / maximum
        effective <- scaled / sum(scaled)
        if (any(input_weight[selected] > 0 & effective == 0)) {
            stop(
                "Positive `weights` cannot underflow during normalization.",
                call. = FALSE
            )
        }
        effective_weight[selected] <- effective
    }
    effective_weight
}

.validate_aggregation_inputs <- function(
    mat,
    gene_sets,
    groups,
    weights,
    missing
) {
    .validate_aggregation_matrix_structure(mat)
    units <- colnames(mat)
    .validate_aggregation_named_vector(groups, units, "character", "`groups`")
    .validate_identifier_vector(groups, "`groups`", allow_duplicates = TRUE)
    .validate_aggregation_named_vector(weights, units, "numeric", "`weights`")
    if (anyNA(weights) || any(!is.finite(weights))) {
        stop("`weights` must contain only finite values.", call. = FALSE)
    }
    if (any(weights < 0)) {
        stop("`weights` must be non-negative.", call. = FALSE)
    }
    .aggregation_effective_weights(groups, weights)
    if (!is.character(missing) || is.object(missing) ||
        !is.null(attributes(missing)) || length(missing) != 1L ||
        is.na(missing) || !missing %in% c("reject", "intersection")) {
        stop(
            "`missing` must be exactly \"reject\" or \"intersection\".",
            call. = FALSE
        )
    }

    prepared <- .prepare_gene_sets(gene_sets, rownames(mat))
    active <- weights > 0
    if (any(active)) {
        .validate_score_matrix(mat[, active, drop = FALSE])
    }
    prepared
}

.aggregation_weight_table <- function(groups, weights, units) {
    input_weight <- as.double(weights)
    effective_weight <- .aggregation_effective_weights(groups, weights)
    data.frame(
        group = unname(groups),
        unit = unname(units),
        input_weight = input_weight,
        effective_weight = effective_weight,
        active = !is.na(effective_weight) & input_weight > 0,
        row.names = seq_along(units),
        check.names = FALSE,
        stringsAsFactors = FALSE
    )
}

.aggregation_dense_block <- function(mat, rows, columns) {
    block <- mat[rows, columns, drop = FALSE]
    if (inherits(block, "Matrix")) {
        block <- as.matrix(block)
    }
    if (typeof(block) == "integer") {
        block <- matrix(
            as.double(block),
            nrow = nrow(block),
            ncol = ncol(block),
            dimnames = dimnames(block)
        )
    }
    block
}

.aggregation_gap_formula <- function(block, weights) {
    member_count <- nrow(block)
    value_scale <- max(block)
    if (value_scale == 0) {
        scaled <- block
        aggregate <- rep.int(0, member_count)
        gap <- 0
        gap_scaled <- 0
    } else {
        scaled <- block / value_scale
        aggregate_scaled <- rowSums(sweep(scaled, 2L, weights, "*"))
        aggregate <- aggregate_scaled * value_scale
        unit_scaled <- t(scaled)
        centered <- sweep(unit_scaled, 1L, rowMeans(unit_scaled), "-")
        positive_mass <- colSums(sweep(
            pmax(centered, 0), 1L, weights, "*"
        ))
        negative_mass <- colSums(sweep(
            pmax(-centered, 0), 1L, weights, "*"
        ))
        coefficient <- member_count / (2 * (member_count - 1L))
        gap_scaled <- 2 * coefficient *
            sum(pmin(positive_mass, negative_mass))
        gap <- gap_scaled * value_scale
    }
    list(
        aggregate = aggregate,
        gap = gap,
        gap_scaled = gap_scaled,
        value_scale = value_scale
    )
}

.aggregation_native_scores <- function(block, aggregate) {
    member_count <- nrow(block)
    unit_count <- ncol(block)
    outcome <- tryCatch(
        list(
            value = calculateScoresDense(
                cbind(block, aggregate),
                list(seq_len(member_count))
            ),
            error = NULL
        ),
        error = function(condition) list(value = NULL, error = condition)
    )
    if (!is.null(outcome$error)) {
        if (!startsWith(
            conditionMessage(outcome$error),
            "Internal GeneFunnel arithmetic produced"
        )) {
            stop(outcome$error)
        }
        return(NULL)
    }
    native <- outcome$value
    if (is.null(native) || !is.matrix(native) ||
        !identical(dim(native), c(1L, unit_count + 1L)) ||
        anyNA(native) || any(!is.finite(native))) {
        return(NULL)
    }
    list(
        aggregate_score = native[[unit_count + 1L]],
        unit_scores = as.double(native[1L, seq_len(unit_count)])
    )
}

.aggregation_certify_scores <- function(formula, native, weights) {
    member_count <- length(formula$aggregate)
    unit_count <- length(weights)
    weighted_unit_score <- sum(weights * native$unit_scores)
    identity_residual <- native$aggregate_score - weighted_unit_score -
        formula$gap
    tolerance <- 256 * .Machine$double.eps *
        max(1, unit_count, member_count) *
        max(
            native$aggregate_score,
            weighted_unit_score,
            formula$gap,
            formula$value_scale
        )
    score_tolerance <- 256 * .Machine$double.eps *
        max(1, unit_count, member_count) *
        max(native$aggregate_score, weighted_unit_score, formula$gap)
    certification_tolerance <- min(tolerance, score_tolerance)
    underflowed_gap <- formula$gap_scaled > 0 && formula$gap == 0
    valid <- all(is.finite(c(
        native$aggregate_score,
        weighted_unit_score,
        formula$gap,
        identity_residual,
        tolerance,
        certification_tolerance
    ))) && !underflowed_gap && formula$gap >= 0 &&
        formula$gap <= native$aggregate_score + certification_tolerance &&
        abs(identity_residual) <= certification_tolerance
    if (!valid) {
        return(NULL)
    }
    list(
        aggregate_score = native$aggregate_score,
        unit_scores = native$unit_scores,
        weighted_unit_score = weighted_unit_score,
        aggregation_gap = formula$gap,
        identity_residual = identity_residual
    )
}

.aggregation_score_block <- function(block, weights) {
    formula <- .aggregation_gap_formula(block, weights)
    native <- .aggregation_native_scores(block, formula$aggregate)
    if (is.null(native)) {
        return(NULL)
    }
    scored <- .aggregation_certify_scores(formula, native, weights)
    if (is.null(scored)) {
        return(NULL)
    }
    normalized_gap <- if (scored$aggregate_score > 0) {
        # Preserve the mathematical [0, 1] range when two independently
        # rounded identities differ only inside the locked allowance.
        if (scored$aggregation_gap > scored$aggregate_score) {
            1
        } else {
            scored$aggregation_gap / scored$aggregate_score
        }
    } else {
        NA_real_
    }
    if (!is.na(normalized_gap) &&
        (!is.finite(normalized_gap) || normalized_gap < 0 ||
            normalized_gap > 1)) {
        return(NULL)
    }
    list(
        aggregate_score = scored$aggregate_score,
        unit_scores = scored$unit_scores,
        weighted_unit_score = scored$weighted_unit_score,
        aggregation_gap = scored$aggregation_gap,
        normalized_gap = normalized_gap,
        normalized_status = if (scored$aggregate_score > 0) {
            "defined"
        } else {
            "zero_aggregate"
        },
        identity_residual = scored$identity_residual
    )
}

.aggregation_empty_removed_members <- function() {
    data.frame(
        group = character(),
        gene_set = character(),
        member = character(),
        unit = character(),
        reason = character(),
        stringsAsFactors = FALSE
    )
}

.aggregation_empty_unit_scores <- function() {
    data.frame(
        group = character(),
        gene_set = character(),
        unit = character(),
        effective_weight = double(),
        score = double(),
        stringsAsFactors = FALSE
    )
}

.aggregation_bind_rows <- function(chunks, empty = NULL) {
    result <- if (length(chunks) == 0L) {
        empty
    } else {
        do.call(rbind, chunks)
    }
    rownames(result) <- seq_len(nrow(result))
    result
}

.aggregation_removed_rows <- function(
    group, set, members, mapping, missing_values = NULL,
    active_units = integer(), units = character()
) {
    chunks <- vector("list", length(members))
    matched_index <- 0L
    for (member_index in seq_along(members)) {
        if (mapping[[member_index]] == 0L) {
            chunks[[member_index]] <- data.frame(
                group = group, gene_set = set,
                member = members[[member_index]], unit = NA_character_,
                reason = "unmatched", stringsAsFactors = FALSE
            )
            next
        }
        matched_index <- matched_index + 1L
        if (is.null(missing_values)) {
            next
        }
        affected <- which(missing_values[matched_index, ])
        if (length(affected) > 0L) {
            chunks[[member_index]] <- data.frame(
                group = rep.int(group, length(affected)),
                gene_set = rep.int(set, length(affected)),
                member = rep.int(members[[member_index]], length(affected)),
                unit = units[active_units[affected]],
                reason = rep.int("missing", length(affected)),
                stringsAsFactors = FALSE
            )
        }
    }
    present <- !vapply(chunks, is.null, logical(1L))
    .aggregation_bind_rows(
        chunks[present],
        .aggregation_empty_removed_members()
    )
}

.aggregation_set_support <- function(
    mat, members, mapping, active_units, missing, units, group, set
) {
    matched_rows <- mapping[mapping > 0L]
    result <- list(
        reason = NULL,
        block = NULL,
        retained_size = as.integer(length(matched_rows)),
        removed_member_count = as.integer(sum(mapping == 0L)),
        removed = .aggregation_removed_rows(group, set, members, mapping)
    )
    if (length(active_units) == 0L) {
        result$reason <- "no_positive_weight"
        return(result)
    }
    if (length(matched_rows) < 2L) {
        result$reason <- "too_few_matched_members"
        return(result)
    }
    block <- .aggregation_dense_block(mat, matched_rows, active_units)
    missing_values <- is.na(block)
    missing_rows <- rowSums(missing_values) > 0L
    result$removed <- .aggregation_removed_rows(
        group, set, members, mapping, missing_values, active_units, units
    )
    result$removed_member_count <- as.integer(
        result$removed_member_count + sum(missing_rows)
    )
    if (any(missing_rows)) {
        if (missing == "reject") {
            result$reason <- "active_missing_values"
            return(result)
        }
        block <- block[!missing_rows, , drop = FALSE]
        result$retained_size <- as.integer(nrow(block))
    }
    if (nrow(block) < 2L) {
        result$reason <- "too_few_common_members"
        return(result)
    }
    result$block <- block
    result
}

.aggregation_evaluate_set <- function(
    mat, members, mapping, active_units, effective, missing, units, group, set
) {
    result <- .aggregation_set_support(
        mat, members, mapping, active_units, missing, units, group, set
    )
    if (!is.null(result$reason)) {
        return(result)
    }
    result$scored <- .aggregation_score_block(result$block, effective)
    result$block <- NULL
    if (is.null(result$scored)) {
        result$reason <- "numerically_unavailable"
    } else {
        result$reason <- "eligible"
    }
    result
}

.aggregation_summary_row <- function(
    group, set, coverage, result, active_count, excluded_count
) {
    eligible <- identical(result$reason, "eligible")
    value <- function(name) {
        if (eligible) result$scored[[name]] else NA_real_
    }
    data.frame(
        group = group, gene_set = set, eligible = eligible,
        reason = result$reason,
        aggregate_score = value("aggregate_score"),
        weighted_unit_score = value("weighted_unit_score"),
        aggregation_gap = value("aggregation_gap"),
        normalized_gap = value("normalized_gap"),
        normalized_status = if (eligible) {
            result$scored$normalized_status
        } else {
            "ineligible"
        },
        identity_residual = value("identity_residual"),
        declared_size = as.integer(coverage$declared_size),
        matched_size = as.integer(coverage$matched_size),
        retained_size = result$retained_size,
        active_unit_count = as.integer(active_count),
        excluded_unit_count = as.integer(excluded_count),
        removed_member_count = result$removed_member_count,
        check.names = FALSE, stringsAsFactors = FALSE
    )
}

.aggregation_unit_score_rows <- function(
    group, set, active_units, units, effective, scored
) {
    count <- length(active_units)
    data.frame(
        group = rep.int(group, count),
        gene_set = rep.int(set, count),
        unit = units[active_units],
        effective_weight = effective,
        score = scored$unit_scores,
        stringsAsFactors = FALSE
    )
}

.aggregation_audit_group <- function(
    mat, prepared, mappings, groups, weight_table, missing, units, group
) {
    group_units <- which(groups == group)
    active_units <- group_units[weight_table$active[group_units]]
    effective <- weight_table$effective_weight[active_units]
    active_count <- length(active_units)
    excluded_count <- length(group_units) - active_count
    set_count <- length(prepared$members)
    summary_chunks <- vector("list", set_count)
    removed_chunks <- list()
    unit_score_chunks <- list()
    for (set_index in seq_len(set_count)) {
        set <- names(prepared$members)[[set_index]]
        result <- .aggregation_evaluate_set(
            mat, prepared$members[[set_index]], mappings[[set_index]],
            active_units, effective, missing, units, group, set
        )
        summary_chunks[[set_index]] <- .aggregation_summary_row(
            group, set, prepared$coverage[set_index, , drop = FALSE],
            result, active_count, excluded_count
        )
        if (nrow(result$removed) > 0L) {
            removed_chunks[[length(removed_chunks) + 1L]] <- result$removed
        }
        if (identical(result$reason, "eligible")) {
            unit_score_chunks[[length(unit_score_chunks) + 1L]] <-
                .aggregation_unit_score_rows(
                    group, set, active_units, units, effective, result$scored
                )
        }
    }
    list(
        summary = .aggregation_bind_rows(summary_chunks),
        removed = .aggregation_bind_rows(
            removed_chunks, .aggregation_empty_removed_members()
        ),
        unit_scores = .aggregation_bind_rows(
            unit_score_chunks, .aggregation_empty_unit_scores()
        )
    )
}

.aggregation_audit <- function(
    mat,
    gene_sets,
    groups,
    weights,
    missing = "reject"
) {
    prepared <- .validate_aggregation_inputs(
        mat, gene_sets, groups, weights, missing
    )
    units <- colnames(mat)
    group_order <- unique(unname(groups))
    weight_table <- .aggregation_weight_table(groups, weights, units)
    features <- rownames(mat)
    mappings <- lapply(prepared$members, function(members) {
        unname(match(members, features, nomatch = 0L))
    })
    audited <- lapply(group_order, function(group) {
        .aggregation_audit_group(
            mat, prepared, mappings, groups, weight_table, missing, units, group
        )
    })
    list(
        summary = .aggregation_bind_rows(lapply(audited, `[[`, "summary")),
        weights = weight_table,
        removed_members = .aggregation_bind_rows(
            lapply(audited, `[[`, "removed"),
            .aggregation_empty_removed_members()
        ),
        unit_scores = .aggregation_bind_rows(
            lapply(audited, `[[`, "unit_scores"),
            .aggregation_empty_unit_scores()
        )
    )
}
