# Assisted-by: OpenAI Codex.

aggregation_set_seed <- function(seed) {
    valid <- length(seed) == 1L && is.finite(seed) && seed == floor(seed) &&
        seed > 0 && seed <= .Machine$integer.max
    if (!valid) {
        stop("Synthetic seed is invalid.", call. = FALSE)
    }
    RNGkind(
        kind = "Mersenne-Twister",
        normal.kind = "Inversion",
        sample.kind = "Rejection"
    )
    set.seed(as.integer(seed))
}

aggregation_latent_baseline <- function(t, profile, dynamic_range) {
    baseline <- switch(
        profile,
        monotone = dynamic_range^(t - 0.5),
        u_shaped = dynamic_range^(2 * abs(t - 0.5) - 0.5),
        stop("Unknown synthetic baseline profile.", call. = FALSE)
    )
    baseline / mean(baseline)
}

aggregation_archetype_pattern <- function(t, unit_count, archetype) {
    member_count <- length(t)
    unit_index <- seq_len(unit_count)
    common <- sin(2 * pi * t)
    pattern <- switch(
        archetype,
        complex_like = matrix(
            rep(common, times = unit_count),
            nrow = unit_count,
            byrow = TRUE
        ),
        cascade_like = t(vapply(unit_index, function(unit) {
            peak <- exp(-16 * (t - (unit - 0.5) / unit_count)^2)
            peak - mean(peak)
        }, numeric(member_count))),
        regulatory_like = t(vapply(unit_index, function(unit) {
            cos(2 * pi * (t - (unit - 1) / unit_count))
        }, numeric(member_count))),
        mixed_direction = outer(
            (-1)^(unit_index - 1L),
            ifelse(seq_len(member_count) <= member_count / 2L, 1, -1)
        ),
        stop("Unknown synthetic archetype.", call. = FALSE)
    )
    sweep(pattern, 1L, apply(abs(pattern), 1L, max), "/")
}

aggregation_overlap_multiplier <- function(
    member_count,
    unit_count,
    overlap_fraction,
    archetype
) {
    multiplier <- matrix(1, nrow = unit_count, ncol = member_count)
    if (overlap_fraction == 0) {
        return(multiplier)
    }
    affected <- seq_len(member_count / 2L)
    member_sign <- (-1)^(affected - 1L)
    unit_sign <- if (identical(archetype, "complex_like")) {
        rep.int(1, unit_count)
    } else {
        (-1)^(seq_len(unit_count) - 1L)
    }
    multiplier[, affected] <- 1 + 0.5 * overlap_fraction *
        outer(unit_sign, member_sign)
    multiplier
}

aggregation_perturb_profiles <- function(profiles, t, row, archetype) {
    unit_count <- nrow(profiles)
    member_count <- ncol(profiles)
    aggregation_set_seed(row$latent_seed[[1L]])
    unit_normal <- stats::rnorm(unit_count)
    independent <- matrix(
        stats::rnorm(unit_count * member_count),
        nrow = unit_count,
        byrow = TRUE
    )
    if (identical(archetype, "complex_like")) {
        unit_normal[] <- unit_normal[[1L]]
        independent <- matrix(
            rep(independent[1L, ], times = unit_count),
            nrow = unit_count,
            byrow = TRUE
        )
    }
    loading <- 0.5 + t
    correlation <- row$correlation[[1L]]
    log_sd <- row$log_sd[[1L]]
    noise <- sqrt(correlation) * outer(unit_normal, loading) +
        sqrt(1 - correlation) * independent
    correction <- 0.5 * log_sd^2 * matrix(
        correlation * loading^2 + (1 - correlation),
        nrow = unit_count,
        ncol = member_count,
        byrow = TRUE
    )
    profiles * exp(log_sd * noise - correction)
}

aggregation_latent_profiles <- function(row, common_scale) {
    member_count <- as.integer(row$member_count[[1L]])
    unit_count <- as.integer(row$unit_count[[1L]])
    t <- (seq_len(member_count) - 1) / (member_count - 1)
    baseline <- aggregation_latent_baseline(
        t,
        row$baseline_profile[[1L]],
        row$dynamic_range[[1L]]
    )
    archetype <- row$archetype[[1L]]
    pattern <- aggregation_archetype_pattern(t, unit_count, archetype)
    common <- matrix(
        rep(sin(2 * pi * t), times = unit_count),
        nrow = unit_count,
        byrow = TRUE
    )
    mixed <- (1 - row$complementarity[[1L]]) * common +
        row$complementarity[[1L]] * pattern
    profiles <- matrix(baseline, unit_count, member_count, byrow = TRUE) *
        (1 + 0.8 * mixed)
    profiles <- aggregation_perturb_profiles(profiles, t, row, archetype)
    profiles <- profiles * aggregation_overlap_multiplier(
        member_count,
        unit_count,
        row$overlap_fraction[[1L]],
        archetype
    )
    profiles[1L, 1L] <- profiles[1L, 1L] *
        row$outlier_multiplier[[1L]]
    profiles <- sweep(profiles, 1L, rowSums(profiles), "/") * common_scale
    if (any(!is.finite(profiles)) || any(profiles <= 0)) {
        stop(
            "Synthetic latent profiles left the positive domain.",
            call. = FALSE
        )
    }
    profiles
}

aggregation_reference_score <- function(values) {
    member_count <- length(values)
    coefficient <- member_count / (2 * (member_count - 1L))
    sum(values) - coefficient * sum(abs(values - mean(values)))
}

aggregation_reference_audit <- function(profiles, weights) {
    weights <- weights / sum(weights)
    centers <- rowMeans(profiles)
    centered <- sweep(profiles, 1L, centers, "-")
    positive <- colSums(sweep(pmax(centered, 0), 1L, weights, "*"))
    negative <- colSums(sweep(pmax(-centered, 0), 1L, weights, "*"))
    member_count <- ncol(profiles)
    coefficient <- member_count / (2 * (member_count - 1L))
    gap <- 2 * coefficient * sum(pmin(positive, negative))
    aggregate <- colSums(sweep(profiles, 1L, weights, "*"))
    aggregate_score <- aggregation_reference_score(aggregate)
    unit_scores <- apply(profiles, 1L, aggregation_reference_score)
    weighted_unit_score <- sum(weights * unit_scores)
    result <- list(
        aggregate_score = aggregate_score,
        weighted_unit_score = weighted_unit_score,
        aggregation_gap = gap,
        normalized_gap = if (aggregate_score > 0) {
            gap / aggregate_score
        } else {
            NA_real_
        },
        identity_residual = aggregate_score - weighted_unit_score - gap
    )
    tolerance <- 1e-12 * max(1, aggregate_score)
    valid <- all(is.finite(unlist(result))) && result$aggregation_gap >= 0 &&
        result$normalized_gap >= 0 && result$normalized_gap <= 1 &&
        abs(result$identity_residual) <= tolerance
    if (!valid) {
        stop("Synthetic latent reference failed its identity.", call. = FALSE)
    }
    result
}

aggregation_subunit_depths <- function(library_depth, subunit_count) {
    base <- library_depth %/% subunit_count
    remainder <- library_depth %% subunit_count
    base + as.integer(seq_len(subunit_count) <= remainder)
}

aggregation_measure_profiles <- function(profiles, row, common_scale) {
    aggregation_set_seed(row$seed[[1L]])
    unit_count <- nrow(profiles)
    member_count <- ncol(profiles)
    counts <- matrix(0, nrow = unit_count, ncol = member_count)
    depths <- aggregation_subunit_depths(
        as.integer(row$library_depth[[1L]]),
        as.integer(row$subunit_count[[1L]])
    )
    dropout <- row$dropout_fraction[[1L]]
    for (unit in seq_len(unit_count)) {
        probabilities <- profiles[unit, ] / sum(profiles[unit, ])
        for (depth in depths) {
            observed <- as.double(stats::rmultinom(
                1L,
                size = depth,
                prob = probabilities
            ))
            observed[stats::runif(member_count) < dropout] <- 0
            counts[unit, ] <- counts[unit, ] + observed
        }
    }
    totals <- rowSums(counts)
    eligible <- all(totals > 0)
    list(
        eligible = eligible,
        reason = if (eligible) "eligible" else "zero_observed_unit",
        profiles = if (eligible) {
            sweep(counts, 1L, totals, "/") * common_scale
        } else {
            NULL
        },
        raw_total = sum(counts),
        detected_fraction = mean(counts > 0),
        zero_unit_count = as.integer(sum(totals == 0))
    )
}

aggregation_audit_measurements <- function(
    measurements,
    labels,
    weights,
    audit
) {
    eligible <- which(vapply(
        measurements,
        `[[`,
        logical(1L),
        "eligible"
    ))
    if (length(eligible) == 0L) {
        return(NULL)
    }
    member_count <- ncol(measurements[[eligible[[1L]]]]$profiles)
    unit_count <- length(weights)
    blocks <- lapply(eligible, function(index) {
        t(measurements[[index]]$profiles)
    })
    mat <- do.call(cbind, blocks)
    feature_ids <- sprintf("member_%03d", seq_len(member_count))
    unit_ids <- unlist(lapply(labels[eligible], function(label) {
        paste0(label, "_unit_", seq_len(unit_count))
    }), use.names = FALSE)
    dimnames(mat) <- list(feature_ids, unit_ids)
    group_values <- rep(labels[eligible], each = unit_count)
    groups <- stats::setNames(group_values, unit_ids)
    masses <- stats::setNames(rep(weights, times = length(eligible)), unit_ids)
    audit(
        mat,
        list(synthetic = feature_ids),
        groups,
        masses,
        missing = "reject"
    )
}

aggregation_measurement_metrics <- function(
    measurement,
    label,
    audited,
    reference,
    max_weight
) {
    summary <- if (is.null(audited) || !measurement$eligible) {
        NULL
    } else {
        audited$summary[audited$summary$group == label, , drop = FALSE]
    }
    audit_eligible <- !is.null(summary) && nrow(summary) == 1L &&
        isTRUE(summary$eligible[[1L]])
    value <- function(field) {
        if (audit_eligible) summary[[field]][[1L]] else NA_real_
    }
    data.frame(
        eligible = audit_eligible,
        reason = if (measurement$eligible) {
            if (is.null(summary) || nrow(summary) != 1L) {
                "missing_audit_row"
            } else {
                summary$reason[[1L]]
            }
        } else {
            measurement$reason
        },
        zero_unit_count = measurement$zero_unit_count,
        raw_total = measurement$raw_total,
        detected_fraction = measurement$detected_fraction,
        max_weight = max_weight,
        latent_aggregate_score = reference$aggregate_score,
        latent_weighted_unit_score = reference$weighted_unit_score,
        latent_aggregation_gap = reference$aggregation_gap,
        latent_normalized_gap = reference$normalized_gap,
        observed_aggregate_score = value("aggregate_score"),
        observed_weighted_unit_score = value("weighted_unit_score"),
        observed_aggregation_gap = value("aggregation_gap"),
        observed_normalized_gap = value("normalized_gap"),
        identity_residual = value("identity_residual"),
        stringsAsFactors = FALSE,
        check.names = FALSE
    )
}

aggregation_simulate_pair <- function(pair, audit, common_scale) {
    invariant_fields <- setdiff(
        names(pair),
        c("scenario_id", "seed", "measurement_replicate")
    )
    invariant <- vapply(invariant_fields, function(field) {
        identical(pair[[field]][[1L]], pair[[field]][[2L]])
    }, logical(1L))
    if (nrow(pair) != 2L ||
        !identical(pair$measurement_replicate, c("A", "B")) ||
        !all(invariant)) {
        stop("Synthetic A/B pair is malformed.", call. = FALSE)
    }
    latent <- aggregation_latent_profiles(
        pair[1L, , drop = FALSE],
        common_scale
    )
    weights <- aggregation_weights(
        pair$unit_count[[1L]],
        pair$weight_profile[[1L]]
    )
    reference <- aggregation_reference_audit(latent, weights)
    measurements <- lapply(seq_len(2L), function(index) {
        aggregation_measure_profiles(
            latent,
            pair[index, , drop = FALSE],
            common_scale
        )
    })
    labels <- pair$measurement_replicate
    audited <- aggregation_audit_measurements(
        measurements,
        labels,
        weights,
        audit
    )
    metrics <- do.call(rbind, lapply(seq_len(2L), function(index) {
        aggregation_measurement_metrics(
            measurements[[index]],
            labels[[index]],
            audited,
            reference,
            max(weights)
        )
    }))
    rownames(metrics) <- NULL
    cbind(pair, metrics)
}

aggregation_smoke_design <- function(design) {
    latent <- design[design$measurement_replicate == "A", , drop = FALSE]
    queries <- list(
        complex = latent$archetype == "complex_like" &
            latent$outlier_multiplier == 1 & latent$dropout_fraction == 0,
        cascade = latent$archetype == "cascade_like" &
            latent$subunit_count == 16L & latent$unit_count == 4L,
        regulatory = latent$archetype == "regulatory_like" &
            latent$dropout_fraction == 0.5 & latent$library_depth == 100000L,
        mixed = latent$archetype == "mixed_direction" &
            latent$weight_profile == "dominant" &
            latent$outlier_multiplier == 20
    )
    selected <- vapply(queries, function(query) {
        match(TRUE, query, nomatch = 0L)
    }, integer(1L))
    if (any(selected == 0L)) {
        stop("Synthetic smoke strata are absent.", call. = FALSE)
    }
    ids <- latent$latent_id[selected]
    rows <- unlist(lapply(ids, function(id) {
        which(design$latent_id == id)
    }), use.names = FALSE)
    design[rows, , drop = FALSE]
}

aggregation_validate_smoke <- function(observed, design) {
    if (!identical(observed$scenario_id, design$scenario_id) ||
        nrow(observed) != 8L || !all(table(observed$latent_id) == 2L)) {
        stop("Synthetic smoke output identity is invalid.", call. = FALSE)
    }
    pairs <- split(observed, observed$latent_id, drop = TRUE)
    paired_target <- vapply(pairs, function(pair) {
        identical(pair$latent_normalized_gap[[1L]],
            pair$latent_normalized_gap[[2L]])
    }, logical(1L))
    raw_limit <- observed$library_depth * observed$unit_count
    valid_metrics <- observed$raw_total >= 0 & observed$raw_total <= raw_limit &
        observed$detected_fraction >= 0 & observed$detected_fraction <= 1 &
        (!observed$eligible |
            (is.finite(observed$observed_normalized_gap) &
                observed$observed_normalized_gap >= 0 &
                observed$observed_normalized_gap <= 1))
    complex <- observed$archetype == "complex_like"
    if (!all(paired_target) || !all(valid_metrics) ||
        !all(observed$latent_normalized_gap[complex] == 0) ||
        !all(observed$seed[observed$measurement_replicate == "B"] -
            observed$seed[observed$measurement_replicate == "A"] == 1000000L)) {
        stop("Synthetic smoke invariants failed.", call. = FALSE)
    }
    invisible(observed)
}

aggregation_chunk_path <- function(checkpoint_dir, chunk_id) {
    file.path(checkpoint_dir, sprintf("chunk-%04d.rds", chunk_id))
}

aggregation_simulate_chunk <- function(
    chunk_id,
    latent_indices,
    design,
    audit,
    common_scale,
    checkpoint_dir,
    fingerprint
) {
    rows <- as.vector(rbind(2L * latent_indices - 1L, 2L * latent_indices))
    expected_ids <- design$scenario_id[rows]
    path <- if (is.null(checkpoint_dir)) NULL else {
        aggregation_chunk_path(checkpoint_dir, chunk_id)
    }
    if (!is.null(path) && file.exists(path)) {
        cached <- readRDS(path)
        if (!identical(attr(cached, "fingerprint"), fingerprint) ||
            !identical(cached$scenario_id, expected_ids)) {
            stop("Synthetic checkpoint identity mismatch.", call. = FALSE)
        }
        return(cached)
    }
    pairs <- lapply(latent_indices, function(index) {
        pair_rows <- c(2L * index - 1L, 2L * index)
        aggregation_simulate_pair(
            design[pair_rows, , drop = FALSE],
            audit,
            common_scale
        )
    })
    observed <- do.call(rbind, pairs)
    rownames(observed) <- NULL
    attr(observed, "fingerprint") <- fingerprint
    if (!is.null(path)) {
        temporary <- paste0(path, ".tmp-", Sys.getpid())
        saveRDS(observed, temporary, version = 3L)
        if (!file.rename(temporary, path)) {
            unlink(temporary, force = TRUE)
            stop("Cannot publish a synthetic checkpoint.", call. = FALSE)
        }
    }
    observed
}

aggregation_execution_shape <- function(design, workers, chunk_size) {
    paired_rows <- nrow(design) %% 2L == 0L
    expected_replicates <- if (paired_rows) {
        rep(c("A", "B"), nrow(design) / 2L)
    } else {
        character()
    }
    valid <- paired_rows && workers >= 1L && chunk_size >= 1L &&
        identical(design$measurement_replicate, expected_replicates)
    if (!valid) {
        stop(
            "Synthetic execution settings or pair order are invalid.",
            call. = FALSE
        )
    }
    as.integer(nrow(design) / 2L)
}

aggregation_execute_chunks <- function(chunk_sequence, run_chunk, workers) {
    if (workers == 1L || .Platform$OS.type == "windows") {
        return(lapply(chunk_sequence, run_chunk))
    }
    parallel::mclapply(
        chunk_sequence,
        run_chunk,
        mc.cores = workers,
        mc.preschedule = TRUE,
        mc.set.seed = FALSE
    )
}

aggregation_simulate_design <- function(
    design,
    audit,
    common_scale,
    workers = 1L,
    chunk_size = 128L,
    checkpoint_dir = NULL,
    fingerprint = AGGREGATION_PROTOCOL_SHA256
) {
    workers <- as.integer(workers)
    chunk_size <- as.integer(chunk_size)
    latent_count <- aggregation_execution_shape(design, workers, chunk_size)
    if (!is.null(checkpoint_dir)) {
        dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)
    }
    chunk_ids <- split(
        seq_len(latent_count),
        ceiling(seq_len(latent_count) / chunk_size)
    )
    run_chunk <- function(chunk_id) {
        aggregation_simulate_chunk(
            chunk_id,
            chunk_ids[[chunk_id]],
            design,
            audit,
            common_scale,
            checkpoint_dir,
            fingerprint
        )
    }
    chunk_sequence <- seq_along(chunk_ids)
    observed <- aggregation_execute_chunks(
        chunk_sequence,
        run_chunk,
        workers
    )
    if (any(vapply(observed, inherits, logical(1L), "try-error"))) {
        stop("One or more synthetic chunks failed.", call. = FALSE)
    }
    result <- do.call(rbind, observed)
    attr(result, "fingerprint") <- NULL
    rownames(result) <- NULL
    if (!identical(result$scenario_id, design$scenario_id)) {
        stop("Synthetic chunks returned out of order.", call. = FALSE)
    }
    result
}
