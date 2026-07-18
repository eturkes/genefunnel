# Assisted-by: OpenAI Codex.

.sensitivity_result_fields <- function() {
    c(
        "largest_member",
        "largest_absolute_delta",
        "largest_delta",
        "largest_delta_over_sum",
        "median_absolute_delta",
        "effective_size",
        "status"
    )
}

.sensitivity_status_fields <- function() {
    c("semantic", "delta", "normalized")
}

.sensitivity_status_values <- function(field) {
    switch(
        field,
        semantic = c("defined", "too_few_observed"),
        delta = c("ordinary", "unavailable", "not_applicable"),
        normalized = c(
            "ordinary", "zero_total", "unavailable", "not_applicable"
        ),
        stop("Internal sensitivity status field is invalid.", call. = FALSE)
    )
}

.new_sensitivity_result <- function(n_gene_sets, n_columns, dimnames = NULL) {
    new_matrix <- function(value, mode) {
        matrix(
            as.vector(value, mode = mode),
            nrow = n_gene_sets,
            ncol = n_columns,
            dimnames = dimnames
        )
    }
    character_matrix <- function() {
        new_matrix(NA_character_, "character")
    }
    double_matrix <- function() new_matrix(NA_real_, "double")
    list(
        largest_member = character_matrix(),
        largest_absolute_delta = double_matrix(),
        largest_delta = double_matrix(),
        largest_delta_over_sum = double_matrix(),
        median_absolute_delta = double_matrix(),
        effective_size = new_matrix(NA_integer_, "integer"),
        status = list(
            semantic = character_matrix(),
            delta = character_matrix(),
            normalized = character_matrix()
        )
    )
}

.sensitivity_not_applicable <- function(effective_size) {
    list(
        largest_member = NA_character_,
        largest_absolute_delta = NA_real_,
        largest_delta = NA_real_,
        largest_delta_over_sum = NA_real_,
        median_absolute_delta = NA_real_,
        effective_size = as.integer(effective_size),
        semantic = "too_few_observed",
        delta_status = "not_applicable",
        normalized_status = "not_applicable"
    )
}

.sensitivity_exact_summary <- function(exact) {
    magnitudes <- lapply(exact$deltas, `[[`, "magnitude")
    winner <- 1L
    for (index in seq_along(magnitudes)[-1L]) {
        if (.gf_big_compare(magnitudes[[index]], magnitudes[[winner]]) > 0L) {
            winner <- index
        }
    }
    ordered <- .gf_big_order(magnitudes)
    size <- length(magnitudes)
    if (size %% 2L == 1L) {
        median_numerator <- magnitudes[[ordered[[(size + 1L) %/% 2L]]]]
        median_denominator <- exact$denominator
    } else {
        lower <- ordered[[size %/% 2L]]
        upper <- ordered[[size %/% 2L + 1L]]
        median_numerator <- .gf_big_add(
            magnitudes[[lower]],
            magnitudes[[upper]]
        )
        median_denominator <- .gf_big_multiply_small(exact$denominator, 2L)
    }
    list(
        winner = winner,
        largest_numerator = magnitudes[[winner]],
        median_numerator = median_numerator,
        median_denominator = median_denominator
    )
}

.sensitivity_unavailable <- function(effective_size) {
    result <- .sensitivity_not_applicable(effective_size)
    result$semantic <- "defined"
    result$delta_status <- "unavailable"
    result$normalized_status <- "unavailable"
    result
}

.sensitivity_raw_values <- function(exact, summary, size) {
    largest <- .gf_big_ratio_double(
        summary$largest_numerator,
        exact$denominator,
        exact$exponent,
        1L,
        size
    )
    median <- .gf_big_ratio_double(
        summary$median_numerator,
        summary$median_denominator,
        exact$exponent,
        1L,
        size
    )
    if (is.na(largest) || is.na(median)) {
        return(NULL)
    }
    list(largest = largest, median = median)
}

.sensitivity_normalized_value <- function(exact, summary, sign, size) {
    if (.gf_big_is_zero(exact$total)) {
        return(list(value = NA_real_, status = "zero_total"))
    }
    value <- .gf_big_ratio_double(
        summary$largest_numerator,
        .gf_big_multiply(exact$denominator, exact$total),
        0L,
        sign,
        size
    )
    list(
        value = value,
        status = if (is.na(value)) "unavailable" else "ordinary"
    )
}

.sensitivity_from_exact <- function(exact, members) {
    size <- length(members)
    summary <- .sensitivity_exact_summary(exact)
    raw <- .sensitivity_raw_values(exact, summary, size)
    if (is.null(raw)) {
        return(.sensitivity_unavailable(size))
    }
    winner <- summary$winner
    sign <- exact$deltas[[winner]]$sign
    normalized <- .sensitivity_normalized_value(exact, summary, sign, size)
    list(
        largest_member = members[[winner]],
        largest_absolute_delta = raw$largest,
        largest_delta = sign * raw$largest,
        largest_delta_over_sum = normalized$value,
        median_absolute_delta = raw$median,
        effective_size = as.integer(size),
        semantic = "defined",
        delta_status = "ordinary",
        normalized_status = normalized$status
    )
}

.sensitivity_cell <- function(values, members) {
    if (length(values) != length(members)) {
        stop("Internal sensitivity membership is misaligned.", call. = FALSE)
    }
    observed <- !is.na(values)
    values <- unname(values[observed])
    members <- members[observed]
    if (length(values) < 3L) {
        return(.sensitivity_not_applicable(length(values)))
    }
    .sensitivity_from_exact(.gf_exact_deltas_sorted(values), members)
}

.sensitivity_extract_values <- function(mat, indices, column) {
    as.numeric(mat[indices, column, drop = TRUE])
}

.sensitivity_matrix_chunk <- function(mat, gene_indices, gene_members, storage) {
    mat <- switch(
        storage,
        dense = .as_dense_numeric_chunk(mat),
        sparse = .as_sparse_numeric_chunk(mat),
        stop("Internal sensitivity storage kind is invalid.", call. = FALSE)
    )
    result <- .new_sensitivity_result(length(gene_indices), ncol(mat))
    for (gene_set in seq_along(gene_indices)) {
        for (column in seq_len(ncol(mat))) {
            cell <- .sensitivity_cell(
                .sensitivity_extract_values(
                    mat,
                    gene_indices[[gene_set]],
                    column
                ),
                gene_members[[gene_set]]
            )
            result$largest_member[[gene_set, column]] <- cell$largest_member
            result$largest_absolute_delta[[gene_set, column]] <-
                cell$largest_absolute_delta
            result$largest_delta[[gene_set, column]] <- cell$largest_delta
            result$largest_delta_over_sum[[gene_set, column]] <-
                cell$largest_delta_over_sum
            result$median_absolute_delta[[gene_set, column]] <-
                cell$median_absolute_delta
            result$effective_size[[gene_set, column]] <- cell$effective_size
            result$status$semantic[[gene_set, column]] <- cell$semantic
            result$status$delta[[gene_set, column]] <- cell$delta_status
            result$status$normalized[[gene_set, column]] <-
                cell$normalized_status
        }
    }
    result
}

.sensitivity_matrix_task <- function(
    task,
    gene_indices,
    gene_members,
    storage,
    ...
) {
    tryCatch(
        list(
            id = task$id,
            first = task$first,
            last = task$last,
            sensitivity = .sensitivity_matrix_chunk(
                task$mat,
                gene_indices,
                gene_members,
                storage
            )
        ),
        error = function(condition) {
            stop(
                sprintf(
                    "GeneFunnel sensitivity calculation failed in chunk %d (%s): %s",
                    task$id,
                    .format_chunk_context(task),
                    conditionMessage(condition)
                ),
                call. = FALSE
            )
        }
    )
}

.validate_sensitivity_chunk <- function(value, dimensions) {
    valid_schema <- is.list(value) &&
        !is.object(value) &&
        identical(names(value), .sensitivity_result_fields()) &&
        is.list(value$status) &&
        !is.object(value$status) &&
        identical(names(value$status), .sensitivity_status_fields())
    if (!valid_schema) {
        stop("Internal sensitivity schema is invalid.", call. = FALSE)
    }
    for (field in c(
        "largest_absolute_delta", "largest_delta",
        "largest_delta_over_sum", "median_absolute_delta"
    )) {
        .validate_component_matrix(
            value[[field]], "double", dimensions, paste0("`", field, "`")
        )
    }
    .validate_component_matrix(
        value$largest_member, "character", dimensions, "`largest_member`"
    )
    .validate_component_matrix(
        value$effective_size, "integer", dimensions, "`effective_size`"
    )
    for (field in .sensitivity_status_fields()) {
        .validate_component_matrix(
            value$status[[field]], "character", dimensions, field
        )
        if (anyNA(value$status[[field]]) ||
            any(!value$status[[field]] %in% .sensitivity_status_values(field))) {
            stop("Internal sensitivity status is invalid.", call. = FALSE)
        }
    }
    invisible(value)
}

.copy_sensitivity_chunk <- function(result, chunk, n_gene_sets) {
    dimensions <- c(
        as.integer(n_gene_sets),
        chunk$last - chunk$first + 1L
    )
    .validate_sensitivity_chunk(chunk$sensitivity, dimensions)
    columns <- seq.int(chunk$first, chunk$last)
    for (field in .sensitivity_result_fields()[1:6]) {
        result[[field]][, columns] <- chunk$sensitivity[[field]]
    }
    for (field in .sensitivity_status_fields()) {
        result$status[[field]][, columns] <- chunk$sensitivity$status[[field]]
    }
    result
}

.assemble_sensitivity_chunks <- function(
    chunks,
    n_gene_sets,
    n_columns,
    dimnames
) {
    chunks <- .order_component_chunks(
        chunks,
        n_columns,
        calculation = "sensitivity"
    )
    result <- .new_sensitivity_result(n_gene_sets, n_columns, dimnames)
    for (chunk in chunks) {
        result <- .copy_sensitivity_chunk(result, chunk, n_gene_sets)
    }
    result
}

.gene_set_sensitivity <- function(
    mat,
    gene_sets,
    BPPARAM = BiocParallel::bpparam()
) {
    .validate_score_matrix(mat)
    prepared <- .prepare_gene_sets(gene_sets, rownames(mat))
    .validate_bpparam(BPPARAM)
    retained <- prepared$coverage$scoreable
    if (any(!retained)) {
        .warn_unscoreable_sets(prepared$coverage$gene_set[!retained])
    }
    retained_names <- prepared$coverage$gene_set[retained]
    dimnames <- list(retained_names, colnames(mat))
    if (!any(retained)) {
        return(.new_sensitivity_result(0L, ncol(mat), dimnames))
    }
    gene_indices <- prepared$indices[retained]
    features <- rownames(mat)
    gene_members <- lapply(gene_indices, function(indices) features[indices])
    storage <- .matrix_storage(mat)
    ranges <- .column_chunk_ranges(
        ncol(mat),
        BiocParallel::bpnworkers(BPPARAM)
    )
    chunks <- BiocParallel::bpiterate(
        .matrix_chunk_iterator(mat, ranges),
        .sensitivity_matrix_task,
        gene_indices = gene_indices,
        gene_members = gene_members,
        storage = storage,
        BPPARAM = BPPARAM
    )
    .assemble_sensitivity_chunks(
        chunks,
        length(gene_indices),
        ncol(mat),
        dimnames
    )
}
