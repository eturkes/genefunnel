# Assisted-by: OpenAI Codex.

aggregation_cellbench_numbers <- function(registry, key) {
    values <- suppressWarnings(as.numeric(aggregation_registry_vector(
        registry,
        "cellbench",
        key
    )))
    if (anyNA(values) || any(!is.finite(values))) {
        stop("CellBench numeric registry field is invalid: ", key, call. = FALSE)
    }
    values
}

aggregation_cellbench_composition_matrix <- function(registry) {
    labels <- aggregation_registry_vector(
        registry,
        "cellbench",
        "mixture_proportions"
    )
    parsed <- strsplit(labels, "/", fixed = TRUE)
    values <- t(vapply(parsed, function(value) {
        result <- suppressWarnings(as.numeric(value))
        if (length(result) != 3L || anyNA(result) || any(result < 0)) {
            stop("CellBench composition registry is invalid.", call. = FALSE)
        }
        result
    }, numeric(3L)))
    rownames(values) <- labels
    values
}

aggregation_cellbench_validate_contract <- function(registry) {
    keys <- c(
        "curve_gate", "stability_gate", "ratio_zero_rule", "curve_group",
        "stability_unit"
    )
    expected <- c(
        "median_absolute_error<=0.15;quantile90_absolute_error<=0.30",
        "within_platform_split_spearman>=0.75;cross_platform_spearman>=0.60",
        "positive/0=Inf;0/positive=-Inf;0/0=0",
        paste0(
            "platform,composition,mRNA_amount,set;",
            "every_fixed_group_median_and_quantile90_must_pass"
        ),
        paste0(
            "condition_set_median;",
            "192_fixed_composition_amount_set_conditions_per_platform"
        )
    )
    observed <- vapply(keys, function(key) {
        aggregation_registry_value(registry, "cellbench", key)
    }, character(1L))
    if (!identical(unname(observed), expected)) {
        stop("CellBench runner disagrees with B-1.0.3.", call. = FALSE)
    }
    if (!Sys.getlocale("LC_COLLATE") %in% c("C", "POSIX")) {
        stop("CellBench execution requires C collation.", call. = FALSE)
    }
    invisible(registry)
}

aggregation_cellbench_read_counts <- function(path) {
    value <- utils::read.csv(
        path,
        row.names = 1L,
        check.names = FALSE,
        stringsAsFactors = FALSE
    )
    matrix <- as.matrix(value)
    storage.mode(matrix) <- "double"
    valid <- nrow(matrix) >= 128L && ncol(matrix) > 0L &&
        !is.null(rownames(matrix)) && !is.null(colnames(matrix)) &&
        !anyNA(rownames(matrix)) && !anyNA(colnames(matrix)) &&
        all(nzchar(rownames(matrix))) && all(nzchar(colnames(matrix))) &&
        !anyDuplicated(rownames(matrix)) && !anyDuplicated(colnames(matrix)) &&
        !anyNA(matrix) && all(is.finite(matrix)) && all(matrix >= 0) &&
        all(matrix == floor(matrix))
    if (!valid) {
        stop("CellBench count matrix is malformed.", call. = FALSE)
    }
    matrix
}

aggregation_cellbench_match_compositions <- function(proportions, registry) {
    allowed <- aggregation_cellbench_composition_matrix(registry)
    matched <- vapply(seq_len(nrow(proportions)), function(index) {
        candidates <- which(rowSums(
            allowed == matrix(
                proportions[index, ],
                nrow = nrow(allowed),
                ncol = ncol(allowed),
                byrow = TRUE
            )
        ) == ncol(allowed))
        if (length(candidates) == 1L) candidates else NA_integer_
    }, integer(1L))
    if (anyNA(matched)) {
        stop("CellBench composition is outside the registry.", call. = FALSE)
    }
    rownames(allowed)[matched]
}

aggregation_cellbench_assign_split <- function(metadata) {
    metadata$split <- NA_character_
    mixed <- !metadata$pure
    key <- paste(metadata$composition, metadata$mRNA_amount, sep = "::")
    groups <- split(which(mixed), key[mixed], drop = TRUE)
    for (indices in groups) {
        ordered <- indices[order(
            rownames(metadata)[indices],
            method = "radix"
        )]
        metadata$split[ordered] <- ifelse(
            seq_along(ordered) %% 2L == 1L,
            "odd",
            "even"
        )
    }
    metadata
}

aggregation_cellbench_read_metadata <- function(path, counts, registry) {
    value <- utils::read.csv(
        path,
        row.names = 1L,
        check.names = FALSE,
        stringsAsFactors = FALSE
    )
    cell_lines <- aggregation_registry_vector(registry, "cellbench", "cell_lines")
    proportion_fields <- paste0(cell_lines, "_prop")
    required <- c("outliers", proportion_fields, "mRNA_amount")
    valid <- all(required %in% names(value)) &&
        identical(rownames(value), colnames(counts)) &&
        is.logical(value$outliers) && !anyNA(value$outliers) &&
        all(vapply(value[proportion_fields], is.numeric, logical(1L))) &&
        is.numeric(value$mRNA_amount)
    if (!valid) {
        stop("CellBench metadata schema or sample alignment failed.", call. = FALSE)
    }
    proportions <- as.matrix(value[proportion_fields])
    if (anyNA(proportions) || any(!is.finite(proportions)) ||
        any(proportions < 0) || any(rowSums(proportions) <= 0)) {
        stop("CellBench mixture proportions are invalid.", call. = FALSE)
    }
    value$composition <- aggregation_cellbench_match_compositions(
        proportions,
        registry
    )
    value$pure <- rowSums(proportions == 1) == 1L &
        rowSums(proportions == 0) == 2L
    amounts <- aggregation_cellbench_numbers(registry, "mrna_amount")
    if (anyNA(value$mRNA_amount) || any(!value$mRNA_amount %in% amounts)) {
        stop("CellBench mRNA amount is outside the registry.", call. = FALSE)
    }
    value
}

aggregation_cellbench_read_platform <- function(
    count_path,
    metadata_path,
    platform,
    registry
) {
    counts <- aggregation_cellbench_read_counts(count_path)
    metadata <- aggregation_cellbench_read_metadata(
        metadata_path,
        counts,
        registry
    )
    retained <- !metadata$outliers
    counts <- counts[, retained, drop = FALSE]
    metadata <- aggregation_cellbench_assign_split(
        metadata[retained, , drop = FALSE]
    )
    totals <- colSums(counts)
    if (!length(totals) || any(!is.finite(totals)) || any(totals <= 0)) {
        stop("CellBench retained library has zero or invalid total.", call. = FALSE)
    }
    scale <- as.numeric(aggregation_registry_value(
        registry,
        "cellbench",
        "common_scale"
    ))
    cpm <- sweep(counts, 2L, totals, "/") * scale
    list(platform = platform, cpm = cpm, metadata = metadata)
}

aggregation_cellbench_common_genes <- function(first, second) {
    common <- rownames(first$cpm)[rownames(first$cpm) %in% rownames(second$cpm)]
    if (length(common) < 128L || anyDuplicated(common)) {
        stop("CellBench common gene universe is insufficient.", call. = FALSE)
    }
    common
}

aggregation_cellbench_pure_profiles <- function(platform, genes, registry) {
    cell_lines <- aggregation_registry_vector(registry, "cellbench", "cell_lines")
    fields <- paste0(cell_lines, "_prop")
    result <- vapply(seq_along(cell_lines), function(index) {
        selected <- platform$metadata$pure &
            platform$metadata[[fields[[index]]]] == 1
        if (!any(selected)) {
            stop("CellBench pure profile is absent.", call. = FALSE)
        }
        rowMeans(platform$cpm[genes, selected, drop = FALSE])
    }, numeric(length(genes)))
    rownames(result) <- genes
    colnames(result) <- cell_lines
    if (any(!is.finite(result)) || any(result < 0)) {
        stop("CellBench pure profile is invalid.", call. = FALSE)
    }
    result
}

aggregation_cellbench_ratio <- function(first, second) {
    if (length(first) != length(second) || any(first < 0) || any(second < 0)) {
        stop("CellBench ratio inputs are invalid.", call. = FALSE)
    }
    result <- numeric(length(first))
    both_positive <- first > 0 & second > 0
    result[both_positive] <- log2(first[both_positive] / second[both_positive])
    result[first > 0 & second == 0] <- Inf
    result[first == 0 & second > 0] <- -Inf
    result
}

aggregation_cellbench_pair_record <- function(
    profiles,
    first,
    second,
    size
) {
    ratio <- aggregation_cellbench_ratio(profiles[, first], profiles[, second])
    genes <- rownames(profiles)
    first_rank <- order(-ratio, genes, method = "radix")
    second_rank <- order(ratio, genes, method = "radix")
    half <- as.integer(size / 2L)
    first_members <- first_rank[seq_len(half)]
    second_rank <- second_rank[!second_rank %in% first_members]
    second_members <- second_rank[seq_len(half)]
    indices <- c(first_members, second_members)
    list(
        set = paste0("pair_", first, "_", second, "_n", size),
        kind = "pair", first_line = first, second_line = second,
        members = genes[indices],
        side = c(rep(paste0(first, "_over_", second), half),
            rep(paste0(second, "_over_", first), half)),
        metric = c(ratio[first_members], -ratio[second_members])
    )
}

aggregation_cellbench_pair_records <- function(profiles, registry) {
    lines <- aggregation_registry_vector(registry, "cellbench", "cell_lines")
    sizes <- as.integer(aggregation_cellbench_numbers(registry, "set_size"))
    pairs <- utils::combn(lines, 2L, simplify = FALSE)
    records <- list()
    for (pair in pairs) {
        for (size in sizes) {
            records[[length(records) + 1L]] <- aggregation_cellbench_pair_record(
                profiles,
                pair[[1L]],
                pair[[2L]],
                size
            )
        }
    }
    records
}

aggregation_cellbench_complex_records <- function(profiles, registry) {
    sizes <- as.integer(aggregation_cellbench_numbers(registry, "set_size"))
    genes <- rownames(profiles)
    minimum <- apply(profiles, 1L, min)
    ranking <- order(-minimum, genes, method = "radix")
    lapply(sizes, function(size) {
        indices <- ranking[seq_len(size)]
        list(
            set = paste0("complex_n", size), kind = "complex",
            first_line = NA_character_, second_line = NA_character_,
            members = genes[indices], side = rep("minimum_pure_CPM", size),
            metric = minimum[indices]
        )
    })
}

aggregation_cellbench_selection_tables <- function(records) {
    manifest <- do.call(rbind, lapply(records, function(record) {
        data.frame(
            set = record$set, kind = record$kind,
            first_line = record$first_line, second_line = record$second_line,
            declared_size = length(record$members), stringsAsFactors = FALSE,
            check.names = FALSE
        )
    }))
    membership <- do.call(rbind, lapply(records, function(record) {
        data.frame(
            set = record$set, member_rank = seq_along(record$members),
            gene_id = record$members, selection_side = record$side,
            selection_metric = record$metric, stringsAsFactors = FALSE,
            check.names = FALSE
        )
    }))
    rownames(manifest) <- rownames(membership) <- NULL
    list(manifest = manifest, membership = membership)
}

aggregation_cellbench_select_sets <- function(training_profiles, registry) {
    records <- c(
        aggregation_cellbench_pair_records(training_profiles, registry),
        aggregation_cellbench_complex_records(training_profiles, registry)
    )
    sets <- stats::setNames(lapply(records, `[[`, "members"), vapply(
        records,
        `[[`,
        character(1L),
        "set"
    ))
    sizes <- vapply(sets, length, integer(1L))
    if (length(sets) != 12L || anyDuplicated(names(sets)) ||
        any(vapply(sets, anyDuplicated, integer(1L))) ||
        !setequal(unique(sizes), c(8L, 32L, 128L))) {
        stop("CellBench set selection failed.", call. = FALSE)
    }
    c(list(sets = sets), aggregation_cellbench_selection_tables(records))
}

aggregation_cellbench_mixed_compositions <- function(registry) {
    values <- aggregation_cellbench_composition_matrix(registry)
    pure <- rowSums(values == 1) == 1L & rowSums(values == 0) == 2L
    values[!pure, , drop = FALSE]
}

aggregation_cellbench_reference_one <- function(
    platform,
    profiles,
    selection,
    composition,
    registry,
    audit
) {
    weights <- as.numeric(composition) / sum(composition)
    names(weights) <- colnames(profiles)
    label <- rownames(composition)[[1L]]
    groups <- stats::setNames(rep(paste0("reference_", label), ncol(profiles)),
        colnames(profiles))
    value <- audit(
        profiles,
        selection$sets,
        groups,
        weights,
        missing = "reject"
    )$summary
    if (!identical(value$gene_set, names(selection$sets))) {
        stop("CellBench reference audit order failed.", call. = FALSE)
    }
    data.frame(
        platform = platform, composition = label, set = value$gene_set,
        eligible = value$eligible, reason = value$reason,
        aggregate_score = value$aggregate_score,
        weighted_unit_score = value$weighted_unit_score,
        aggregation_gap = value$aggregation_gap,
        reference_R = value$normalized_gap,
        normalized_status = value$normalized_status,
        stringsAsFactors = FALSE, check.names = FALSE
    )
}

aggregation_cellbench_references <- function(
    profiles,
    selection,
    registry,
    audit
) {
    compositions <- aggregation_cellbench_mixed_compositions(registry)
    rows <- list()
    for (platform in names(profiles)) {
        for (index in seq_len(nrow(compositions))) {
            rows[[length(rows) + 1L]] <- aggregation_cellbench_reference_one(
                platform,
                profiles[[platform]],
                selection,
                compositions[index, , drop = FALSE],
                registry,
                audit
            )
        }
    }
    result <- do.call(rbind, rows)
    rownames(result) <- NULL
    result
}

aggregation_cellbench_score_matrix <- function(platform, selection, scorer) {
    mixed <- !platform$metadata$pure
    block <- platform$cpm[, mixed, drop = FALSE]
    value <- scorer(
        block,
        selection$sets,
        BPPARAM = BiocParallel::SerialParam()
    )
    expected <- c(length(selection$sets), sum(mixed))
    if (!identical(dim(value), expected) ||
        !identical(rownames(value), names(selection$sets)) ||
        !identical(colnames(value), colnames(block))) {
        stop("CellBench scorer shape or order failed.", call. = FALSE)
    }
    value
}

aggregation_cellbench_observations <- function(
    platform,
    selection,
    references,
    scorer
) {
    scores <- aggregation_cellbench_score_matrix(platform, selection, scorer)
    grid <- expand.grid(
        set = rownames(scores), sample = colnames(scores),
        KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
    )
    metadata <- platform$metadata[match(grid$sample, rownames(platform$metadata)),
        , drop = FALSE]
    reference <- references[
        match(
            paste(metadata$composition, grid$set, sep = "::"),
            paste(references$composition, references$set, sep = "::")
        ),
        , drop = FALSE
    ]
    set_info <- selection$manifest[match(grid$set, selection$manifest$set), ]
    library_score <- as.vector(scores)
    defined <- is.finite(library_score) & library_score > 0 &
        reference$eligible & is.finite(reference$weighted_unit_score) &
        is.finite(reference$reference_R)
    observed <- rep(NA_real_, nrow(grid))
    observed[defined] <- (
        library_score[defined] - reference$weighted_unit_score[defined]
    ) / library_score[defined]
    data.frame(
        platform = platform$platform, sample = grid$sample,
        composition = metadata$composition, mRNA_amount = metadata$mRNA_amount,
        split = metadata$split, set = grid$set, kind = set_info$kind,
        declared_size = set_info$declared_size, library_score = library_score,
        weighted_pure_unit_score = reference$weighted_unit_score,
        reference_R = reference$reference_R, observed_R = observed,
        absolute_error = abs(observed - reference$reference_R),
        defined = defined, stringsAsFactors = FALSE, check.names = FALSE
    )
}

aggregation_cellbench_fixed_grid <- function(selection, registry) {
    expand.grid(
        platform = c(
            aggregation_registry_value(registry, "cellbench", "training_platform"),
            aggregation_registry_value(registry, "cellbench", "heldout_platform")
        ),
        composition = rownames(aggregation_cellbench_mixed_compositions(registry)),
        mRNA_amount = aggregation_cellbench_numbers(registry, "mrna_amount"),
        set = names(selection$sets), KEEP.OUT.ATTRS = FALSE,
        stringsAsFactors = FALSE
    )
}

aggregation_cellbench_group_key <- function(value) {
    paste(
        value$platform,
        value$composition,
        format(value$mRNA_amount, scientific = FALSE, trim = TRUE),
        value$set,
        sep = "::"
    )
}

aggregation_cellbench_quantile <- function(values, probability) {
    if (!length(values) || any(!is.finite(values))) {
        return(NA_real_)
    }
    unname(stats::quantile(
        values,
        probs = probability,
        names = FALSE,
        type = 8L
    ))
}

aggregation_cellbench_curve_one <- function(value) {
    errors <- value$absolute_error
    complete <- length(errors) > 0L && all(value$defined) &&
        all(is.finite(errors))
    median_error <- if (complete) stats::median(errors) else NA_real_
    quantile_error <- if (complete) {
        aggregation_cellbench_quantile(errors, 0.90)
    } else {
        NA_real_
    }
    data.frame(
        library_count = nrow(value), defined_count = sum(value$defined),
        median_absolute_error = median_error,
        quantile90_absolute_error = quantile_error,
        median_pass = complete && median_error <= 0.15,
        quantile90_pass = complete && quantile_error <= 0.30,
        passed = complete && median_error <= 0.15 && quantile_error <= 0.30,
        stringsAsFactors = FALSE, check.names = FALSE
    )
}

aggregation_cellbench_curve_groups <- function(observations, selection, registry) {
    grid <- aggregation_cellbench_fixed_grid(selection, registry)
    observed_groups <- split(
        observations,
        aggregation_cellbench_group_key(observations),
        drop = TRUE
    )
    keys <- aggregation_cellbench_group_key(grid)
    rows <- lapply(keys, function(key) {
        value <- observed_groups[[key]]
        if (is.null(value)) value <- observations[FALSE, , drop = FALSE]
        aggregation_cellbench_curve_one(value)
    })
    summary <- do.call(rbind, rows)
    rownames(summary) <- NULL
    result <- cbind(grid, summary)
    if (nrow(result) != 384L || anyDuplicated(keys)) {
        stop("CellBench curve grid is invalid.", call. = FALSE)
    }
    result
}

aggregation_cellbench_condition_one <- function(value) {
    summarize <- function(selected) {
        values <- value$observed_R[selected]
        if (length(values) && all(is.finite(values))) median(values) else NA_real_
    }
    complete <- is.finite(value$observed_R)
    odd <- value$split == "odd"
    even <- value$split == "even"
    data.frame(
        library_count = nrow(value), defined_count = sum(complete),
        odd_count = sum(odd), even_count = sum(even),
        complete_median_R = summarize(rep(TRUE, nrow(value))),
        odd_median_R = summarize(odd), even_median_R = summarize(even),
        defined = nrow(value) > 0L && all(complete) && any(odd) && any(even),
        stringsAsFactors = FALSE, check.names = FALSE
    )
}

aggregation_cellbench_condition_medians <- function(
    observations,
    selection,
    registry
) {
    grid <- aggregation_cellbench_fixed_grid(selection, registry)
    observed_groups <- split(
        observations,
        aggregation_cellbench_group_key(observations),
        drop = TRUE
    )
    rows <- lapply(aggregation_cellbench_group_key(grid), function(key) {
        value <- observed_groups[[key]]
        if (is.null(value)) value <- observations[FALSE, , drop = FALSE]
        aggregation_cellbench_condition_one(value)
    })
    result <- cbind(grid, do.call(rbind, rows))
    rownames(result) <- NULL
    result
}

aggregation_cellbench_spearman <- function(first, second) {
    retained <- is.finite(first) & is.finite(second)
    first <- first[retained]
    second <- second[retained]
    if (length(first) < 3L || length(unique(first)) < 2L ||
        length(unique(second)) < 2L) {
        return(NA_real_)
    }
    unname(stats::cor(first, second, method = "spearman"))
}

aggregation_cellbench_endpoint <- function(
    gate,
    endpoint,
    stratum,
    estimate,
    comparison,
    threshold,
    count,
    complete = TRUE
) {
    complete <- isTRUE(complete)
    passed <- complete && is.finite(estimate) && switch(
        comparison,
        `<=` = estimate <= threshold,
        `>=` = estimate >= threshold,
        stop("CellBench endpoint comparison is invalid.", call. = FALSE)
    )
    data.frame(
        gate = gate, endpoint = endpoint, stratum = stratum,
        estimate = estimate, comparison = comparison, threshold = threshold,
        count = as.integer(count), passed = passed, stringsAsFactors = FALSE,
        check.names = FALSE
    )
}

aggregation_cellbench_curve_endpoints <- function(curve_groups) {
    median_complete <- all(is.finite(curve_groups$median_absolute_error))
    quantile_complete <- all(is.finite(curve_groups$quantile90_absolute_error))
    rbind(
        aggregation_cellbench_endpoint(
            "curve", "maximum_group_median_absolute_error", "all_groups",
            if (median_complete) max(curve_groups$median_absolute_error) else NA,
            "<=", 0.15, nrow(curve_groups), all(curve_groups$median_pass)
        ),
        aggregation_cellbench_endpoint(
            "curve", "maximum_group_quantile90_absolute_error", "all_groups",
            if (quantile_complete) max(curve_groups$quantile90_absolute_error) else NA,
            "<=", 0.30, nrow(curve_groups), all(curve_groups$quantile90_pass)
        )
    )
}

aggregation_cellbench_within_endpoint <- function(conditions, platform) {
    selected <- conditions$platform == platform
    value <- conditions[selected, , drop = FALSE]
    aggregation_cellbench_endpoint(
        "stability", "split_spearman", platform,
        aggregation_cellbench_spearman(value$odd_median_R, value$even_median_R),
        ">=", 0.75, nrow(value),
        nrow(value) == 192L && all(value$defined)
    )
}

aggregation_cellbench_cross_endpoint <- function(conditions, registry) {
    first_label <- aggregation_registry_value(
        registry,
        "cellbench",
        "training_platform"
    )
    second_label <- aggregation_registry_value(
        registry,
        "cellbench",
        "heldout_platform"
    )
    first <- conditions[conditions$platform == first_label, , drop = FALSE]
    second <- conditions[conditions$platform == second_label, , drop = FALSE]
    key <- function(value) paste(
        value$composition, value$mRNA_amount, value$set, sep = "::"
    )
    second <- second[match(key(first), key(second)), , drop = FALSE]
    aggregation_cellbench_endpoint(
        "stability", "complete_condition_spearman", "cross_platform",
        aggregation_cellbench_spearman(
            first$complete_median_R,
            second$complete_median_R
        ),
        ">=", 0.60, nrow(first), nrow(first) == 192L &&
            identical(key(first), key(second)) && all(first$defined) &&
            all(second$defined)
    )
}

aggregation_cellbench_endpoints <- function(curve_groups, conditions, registry) {
    platforms <- c(
        aggregation_registry_value(registry, "cellbench", "training_platform"),
        aggregation_registry_value(registry, "cellbench", "heldout_platform")
    )
    result <- rbind(
        aggregation_cellbench_curve_endpoints(curve_groups),
        aggregation_cellbench_within_endpoint(conditions, platforms[[1L]]),
        aggregation_cellbench_within_endpoint(conditions, platforms[[2L]]),
        aggregation_cellbench_cross_endpoint(conditions, registry)
    )
    rownames(result) <- NULL
    result
}

aggregation_cellbench_summary <- function(
    profiles,
    selection,
    observations,
    curve_groups,
    endpoints
) {
    data.frame(
        protocol_version = AGGREGATION_PROTOCOL_VERSION,
        platforms = length(profiles), common_genes = nrow(profiles[[1L]]),
        selected_sets = length(selection$sets),
        mixed_library_set_rows = nrow(observations),
        defined_observations = sum(observations$defined),
        curve_groups = nrow(curve_groups), failed_curve_groups = sum(!curve_groups$passed),
        endpoint_count = nrow(endpoints), failed_endpoints = sum(!endpoints$passed),
        all_cellbench_gates_pass = all(endpoints$passed),
        stringsAsFactors = FALSE, check.names = FALSE
    )
}

aggregation_cellbench_profile_table <- function(profiles) {
    rows <- lapply(names(profiles), function(platform) {
        value <- profiles[[platform]]
        grid <- expand.grid(
            gene_id = rownames(value), cell_line = colnames(value),
            KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
        )
        data.frame(
            platform = platform, grid,
            mean_pure_CPM = as.vector(value), stringsAsFactors = FALSE,
            check.names = FALSE
        )
    })
    result <- do.call(rbind, rows)
    rownames(result) <- NULL
    result
}

aggregation_cellbench_analyze <- function(platforms, registry, audit, scorer) {
    aggregation_cellbench_validate_contract(registry)
    common <- aggregation_cellbench_common_genes(platforms[[1L]], platforms[[2L]])
    profiles <- lapply(platforms, aggregation_cellbench_pure_profiles,
        genes = common, registry = registry)
    selection <- aggregation_cellbench_select_sets(profiles[[1L]], registry)
    references <- aggregation_cellbench_references(
        profiles,
        selection,
        registry,
        audit
    )
    observations <- do.call(rbind, lapply(platforms, function(platform) {
        aggregation_cellbench_observations(
            platform,
            selection,
            references[references$platform == platform$platform, , drop = FALSE],
            scorer
        )
    }))
    rownames(observations) <- NULL
    curves <- aggregation_cellbench_curve_groups(observations, selection, registry)
    conditions <- aggregation_cellbench_condition_medians(
        observations,
        selection,
        registry
    )
    endpoints <- aggregation_cellbench_endpoints(curves, conditions, registry)
    list(
        profiles = profiles, selection = selection, references = references,
        observations = observations, curve_groups = curves,
        conditions = conditions, endpoints = endpoints,
        summary = aggregation_cellbench_summary(
            profiles, selection, observations, curves, endpoints
        )
    )
}

aggregation_cellbench_integer_library <- function(probability, depth = 1e6L) {
    probability <- probability / sum(probability)
    expected <- probability * depth
    value <- floor(expected)
    remainder <- as.integer(depth - sum(value))
    if (remainder > 0L) {
        ranking <- order(-(expected - value), seq_along(value), method = "radix")
        value[ranking[seq_len(remainder)]] <- value[ranking[seq_len(remainder)]] + 1
    }
    as.numeric(value)
}

aggregation_cellbench_smoke_profiles <- function() {
    index <- seq_len(160L)
    value <- cbind(
        H2228 = ifelse(index %% 11L == 0L, 0, 1 + (index %% 23L)^2),
        H1975 = ifelse(index %% 13L == 0L, 0, 1 + ((161L - index) %% 19L)^2),
        HCC827 = ifelse(index %% 17L == 0L, 0, 1 + 50 * abs(sin(index / 7)))
    )
    value <- apply(value, 2L, aggregation_cellbench_integer_library)
    rownames(value) <- sprintf("ENSG%011d", index)
    value
}

aggregation_cellbench_smoke_fixture <- function(registry) {
    profiles <- aggregation_cellbench_smoke_profiles()
    compositions <- aggregation_cellbench_composition_matrix(registry)
    mixed <- aggregation_cellbench_mixed_compositions(registry)
    amounts <- aggregation_cellbench_numbers(registry, "mrna_amount")
    counts <- profiles
    colnames(counts) <- paste0("pure_", colnames(profiles))
    metadata <- data.frame(
        outliers = FALSE, H2228_prop = c(1, 0, 0),
        H1975_prop = c(0, 1, 0), HCC827_prop = c(0, 0, 1),
        mRNA_amount = 30, row.names = colnames(counts), check.names = FALSE
    )
    for (composition in rownames(mixed)) {
        weights <- mixed[composition, ] / sum(mixed[composition, ])
        probability <- drop(profiles %*% weights) / 1e6
        for (amount in amounts) {
            for (replicate in 1:2) {
                sample <- paste("mix", composition, amount, replicate, sep = "_")
                counts <- cbind(counts, aggregation_cellbench_integer_library(probability))
                colnames(counts)[ncol(counts)] <- sample
                row <- data.frame(
                    outliers = FALSE, H2228_prop = mixed[composition, 1L],
                    H1975_prop = mixed[composition, 2L],
                    HCC827_prop = mixed[composition, 3L],
                    mRNA_amount = amount, row.names = sample,
                    check.names = FALSE
                )
                metadata <- rbind(metadata, row)
            }
        }
    }
    list(counts = counts, metadata = metadata)
}

aggregation_cellbench_write_smoke <- function(value, directory, prefix) {
    count_path <- file.path(directory, paste0(prefix, ".count.csv.gz"))
    metadata_path <- file.path(directory, paste0(prefix, ".metadata.csv.gz"))
    count_connection <- gzfile(count_path, open = "wt")
    utils::write.csv(value$counts, count_connection, quote = TRUE)
    close(count_connection)
    metadata_connection <- gzfile(metadata_path, open = "wt")
    utils::write.csv(value$metadata, metadata_connection, quote = TRUE)
    close(metadata_connection)
    c(counts = count_path, metadata = metadata_path)
}

aggregation_cellbench_smoke_platforms <- function(registry, directory) {
    fixture <- aggregation_cellbench_smoke_fixture(registry)
    labels <- c(
        aggregation_registry_value(registry, "cellbench", "training_platform"),
        aggregation_registry_value(registry, "cellbench", "heldout_platform")
    )
    lapply(labels, function(label) {
        paths <- aggregation_cellbench_write_smoke(fixture, directory, label)
        aggregation_cellbench_read_platform(
            paths[["counts"]],
            paths[["metadata"]],
            label,
            registry
        )
    }) |> stats::setNames(labels)
}

aggregation_cellbench_expect_error <- function(expression, pattern) {
    message <- tryCatch(
        {
            force(expression)
            NA_character_
        },
        error = conditionMessage
    )
    if (is.na(message) || !grepl(pattern, message, fixed = TRUE)) {
        stop("CellBench malformed-input smoke did not fail closed.", call. = FALSE)
    }
    invisible(message)
}

aggregation_cellbench_validate_input_smoke <- function(
    fixture,
    directory,
    registry
) {
    negative <- fixture
    negative$counts[1L, 1L] <- -1
    paths <- aggregation_cellbench_write_smoke(negative, directory, "negative")
    aggregation_cellbench_expect_error(
        aggregation_cellbench_read_platform(
            paths[["counts"]], paths[["metadata"]], "negative", registry
        ),
        "count matrix is malformed"
    )
    misaligned <- fixture
    misaligned$metadata <- misaligned$metadata[rev(seq_len(nrow(
        misaligned$metadata
    ))), , drop = FALSE]
    paths <- aggregation_cellbench_write_smoke(misaligned, directory, "misaligned")
    aggregation_cellbench_expect_error(
        aggregation_cellbench_read_platform(
            paths[["counts"]], paths[["metadata"]], "misaligned", registry
        ),
        "metadata schema or sample alignment failed"
    )
    invisible(TRUE)
}

aggregation_validate_cellbench_smoke <- function(registry, audit, scorer) {
    old_locale <- Sys.getlocale("LC_COLLATE")
    on.exit(suppressWarnings(Sys.setlocale("LC_COLLATE", old_locale)), add = TRUE)
    if (is.na(suppressWarnings(Sys.setlocale("LC_COLLATE", "C")))) {
        stop("CellBench smoke cannot set C collation.", call. = FALSE)
    }
    directory <- tempfile("aggregation-cellbench-smoke-")
    dir.create(directory)
    on.exit(unlink(directory, recursive = TRUE, force = TRUE), add = TRUE)
    aggregation_cellbench_validate_input_smoke(
        aggregation_cellbench_smoke_fixture(registry),
        directory,
        registry
    )
    platforms <- aggregation_cellbench_smoke_platforms(registry, directory)
    result <- aggregation_cellbench_analyze(platforms, registry, audit, scorer)
    ratio <- aggregation_cellbench_ratio(c(1, 0, 0), c(0, 1, 0))
    stopifnot(
        identical(ratio, c(Inf, -Inf, 0)),
        nrow(result$selection$manifest) == 12L,
        nrow(result$selection$membership) == 672L,
        nrow(result$observations) == 768L,
        nrow(result$curve_groups) == 384L,
        nrow(result$conditions) == 384L,
        nrow(result$endpoints) == 5L,
        all(result$endpoints$passed),
        isTRUE(result$summary$all_cellbench_gates_pass)
    )
    first <- result$curve_groups[1L, ]
    keep <- aggregation_cellbench_group_key(result$observations) !=
        aggregation_cellbench_group_key(first)
    broken <- aggregation_cellbench_curve_groups(
        result$observations[keep, , drop = FALSE],
        result$selection,
        registry
    )
    stopifnot(!broken$passed[[1L]], broken$library_count[[1L]] == 0L)
    invisible(result)
}
