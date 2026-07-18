# Assisted-by: OpenAI Codex.

sensitivity_controlled_set_seed <- function(seed) {
    valid <- length(seed) == 1L && is.finite(seed) && seed == floor(seed) &&
        seed > 0 && seed <= .Machine$integer.max
    if (!valid) {
        stop("Sensitivity controlled seed is invalid.", call. = FALSE)
    }
    RNGkind("Mersenne-Twister", "Inversion", "Rejection")
    set.seed(as.integer(seed))
}

sensitivity_controlled_seed <- function(registry, key, scenario_id) {
    base <- as.integer(sensitivity_registry_value(
        registry, "controlled", key
    ))
    seed <- base + as.integer(scenario_id) - 1L
    if (is.na(seed) || seed > .Machine$integer.max) {
        stop("Sensitivity controlled scenario seed overflowed.", call. = FALSE)
    }
    seed
}

sensitivity_controlled_profile <- function(t, archetype, dynamic_range) {
    switch(
        archetype,
        monotone = dynamic_range^(t - 0.5),
        u_shaped = dynamic_range^(2 * abs(t - 0.5) - 0.5),
        wave = dynamic_range^(sin(2 * pi * t) / 2),
        dominant = ifelse(
            seq_along(t) == 1L,
            sqrt(dynamic_range),
            dynamic_range^(-0.5)
        ),
        stop("Sensitivity controlled archetype is invalid.", call. = FALSE)
    )
}

sensitivity_controlled_latent <- function(row, registry) {
    size <- as.integer(row$member_count[[1L]])
    t <- (seq_len(size) - 1) / (size - 1)
    range <- as.numeric(row$dynamic_range[[1L]])
    sigma <- as.numeric(row$log_sd[[1L]])
    sensitivity_controlled_set_seed(sensitivity_controlled_seed(
        registry, "latent_seed_base", row$scenario_id[[1L]]
    ))
    profile <- sensitivity_controlled_profile(
        t, row$archetype[[1L]], range
    ) * exp(sigma * stats::rnorm(size) - sigma^2 / 2)
    shares <- profile / sum(profile)
    capture <- 2^sin(4 * pi * t)
    depth <- as.numeric(row$expected_total[[1L]])
    means <- depth * shares * capture / sum(shares * capture)
    phi <- as.numeric(row$dispersion[[1L]])
    zero_probability <- if (phi == 0) exp(-means) else (1 + phi * means)^(-1 / phi)
    detection <- (1 - as.numeric(row$dropout_fraction[[1L]])) *
        (1 - zero_probability)
    valid <- all(is.finite(c(profile, shares, means, detection))) &&
        all(profile > 0) && all(shares > 0) && all(means > 0) &&
        all(detection >= 0 & detection <= 1) &&
        abs(sum(shares) - 1) <= 16 * .Machine$double.eps &&
        abs(sum(means) - depth) <= 64 * .Machine$double.eps * depth
    if (!valid) stop("Sensitivity controlled latent state is invalid.", call. = FALSE)
    list(
        members = sprintf("member_%03d", seq_len(size)),
        shares = shares,
        means = means,
        detection = detection
    )
}

sensitivity_controlled_measurement <- function(row, latent, registry, key) {
    sensitivity_controlled_set_seed(sensitivity_controlled_seed(
        registry, key, row$scenario_id[[1L]]
    ))
    size <- length(latent$members)
    phi <- as.numeric(row$dispersion[[1L]])
    counts <- if (phi == 0) {
        stats::rpois(size, latent$means)
    } else {
        stats::rnbinom(size, size = 1 / phi, mu = latent$means)
    }
    dropout <- stats::rbinom(
        size, size = 1L, prob = as.numeric(row$dropout_fraction[[1L]])
    )
    counts[dropout == 1] <- 0
    if (any(!is.finite(counts)) || any(counts < 0) || any(counts != floor(counts))) {
        stop("Sensitivity controlled measurement is invalid.", call. = FALSE)
    }
    as.numeric(counts)
}

sensitivity_controlled_removed_count <- function(size, fraction) {
    as.integer(max(1, min(size - 3L, floor(fraction * size + 0.5))))
}

sensitivity_controlled_mask_weights <- function(mechanism, latent) {
    weights <- switch(
        mechanism,
        uniform = rep.int(1, length(latent$shares)),
        low_abundance = 1 / latent$shares,
        high_abundance = latent$shares,
        low_detection = 1 - latent$detection + 1 / 16,
        high_detection = latent$detection + 1 / 16,
        stop("Sensitivity controlled mask mechanism is invalid.", call. = FALSE)
    )
    if (any(!is.finite(weights)) || any(weights <= 0)) {
        stop("Sensitivity controlled mask weights are invalid.", call. = FALSE)
    }
    weights / max(weights)
}

sensitivity_controlled_mask_seed <- function(row, mask, registry) {
    base <- as.integer(sensitivity_registry_value(
        registry, "controlled", "mask_seed_base"
    ))
    seed <- base + 1000L * (as.integer(row$scenario_id[[1L]]) - 1L) +
        100L * (mask$fraction_index[[1L]] - 1L) +
        10L * (mask$mechanism_index[[1L]] - 1L) +
        mask$mask_repeat[[1L]] - 1L
    if (is.na(seed) || seed > .Machine$integer.max) {
        stop("Sensitivity controlled mask seed overflowed.", call. = FALSE)
    }
    seed
}

sensitivity_controlled_masks <- function(row, latent, registry) {
    design <- sensitivity_mask_design(registry)
    selected <- vector("list", nrow(design))
    seeds <- integer(nrow(design))
    removed <- integer(nrow(design))
    for (index in seq_len(nrow(design))) {
        mask <- design[index, , drop = FALSE]
        seeds[[index]] <- sensitivity_controlled_mask_seed(row, mask, registry)
        removed[[index]] <- sensitivity_controlled_removed_count(
            length(latent$members), mask$removed_fraction[[1L]]
        )
        sensitivity_controlled_set_seed(seeds[[index]])
        selected[[index]] <- sort(sample.int(
            length(latent$members), removed[[index]], replace = FALSE,
            prob = sensitivity_controlled_mask_weights(
                mask$mask_mechanism[[1L]], latent
            )
        ))
    }
    design$mask_seed <- seeds
    design$removed_count <- removed
    design$removed_members <- vapply(selected, function(indices) {
        paste(latent$members[indices], collapse = ";")
    }, character(1L))
    list(design = design, selected = selected)
}

sensitivity_controlled_reference_score <- function(values) {
    values <- values[!is.na(values)]
    size <- length(values)
    sum(values) - size / (2 * (size - 1L)) * sum(abs(values - mean(values)))
}

sensitivity_controlled_scores <- function(
    first, second, masks, members, scorer, backend
) {
    partials <- vapply(masks$selected, function(indices) {
        value <- first
        value[indices] <- NA_real_
        value
    }, numeric(length(first)))
    mat <- cbind(first, second, partials)
    dimnames(mat) <- list(
        members,
        c("full_A", "full_B", sprintf("mask_%02d", seq_len(ncol(partials))))
    )
    observed <- scorer(mat, list(controlled = members), backend)
    scores <- drop(observed[1L, , drop = TRUE])
    reference <- apply(mat, 2L, sensitivity_controlled_reference_score)
    tolerance <- 1e-12 * pmax(1, abs(reference))
    valid <- identical(dim(observed), c(1L, ncol(mat))) &&
        all(is.finite(scores)) && all(abs(scores - reference) <= tolerance)
    if (!valid) stop("Sensitivity controlled scores are invalid.", call. = FALSE)
    list(scores = scores, partials = partials)
}

sensitivity_controlled_sensitivity_facts <- function(cell, observed_sum) {
    zero <- observed_sum == 0
    valid <- identical(cell$semantic, "defined") &&
        identical(cell$delta_status, "ordinary") &&
        identical(cell$normalized_status, if (zero) "zero_total" else "ordinary") &&
        is.finite(cell$largest_absolute_delta) && is.finite(cell$largest_delta) &&
        is.finite(cell$median_absolute_delta) &&
        (zero || is.finite(cell$largest_delta_over_sum))
    if (!valid) {
        stop("Sensitivity controlled diagnostic status is invalid.", call. = FALSE)
    }
    data.frame(
        largest_member = cell$largest_member,
        largest_absolute_delta = cell$largest_absolute_delta,
        largest_delta = cell$largest_delta,
        median_absolute_delta = cell$median_absolute_delta,
        sensitivity_delta_status = cell$delta_status,
        sensitivity_normalized_status = cell$normalized_status,
        largest_absolute_delta_over_sum = if (zero) 0 else {
            abs(cell$largest_delta_over_sum)
        },
        largest_delta_over_sum = if (zero) 0 else cell$largest_delta_over_sum,
        median_absolute_delta_over_sum = if (zero) 0 else {
            cell$median_absolute_delta / observed_sum
        },
        stringsAsFactors = FALSE,
        check.names = FALSE
    )
}

sensitivity_controlled_scenario_row <- function(row) {
    data.frame(
        scenario_id = as.integer(row$scenario_id[[1L]]),
        fold = as.integer(row$fold[[1L]]),
        member_count = as.integer(row$member_count[[1L]]),
        archetype = row$archetype[[1L]],
        dynamic_range = as.numeric(row$dynamic_range[[1L]]),
        log_sd = as.numeric(row$log_sd[[1L]]),
        expected_total = as.numeric(row$expected_total[[1L]]),
        dispersion = as.numeric(row$dispersion[[1L]]),
        dropout_fraction = as.numeric(row$dropout_fraction[[1L]]),
        profile_replicate = as.integer(row$profile_replicate[[1L]]),
        stringsAsFactors = FALSE,
        check.names = FALSE
    )
}

sensitivity_controlled_baseline <- function(
    score, observed_sum, declared_size, coverage, observed_fraction, effective
) {
    score_ratio <- numeric(length(score))
    positive <- observed_sum > 0
    score_ratio[positive] <- score[positive] / observed_sum[positive]
    data.frame(
        declared_size = as.integer(declared_size),
        declared_coverage = coverage,
        observed_fraction = observed_fraction,
        observed_sum = observed_sum,
        effective_size = as.integer(effective),
        score = score,
        log2_declared_size = log2(declared_size),
        log1p_observed_sum = log1p(observed_sum),
        score_over_sum = score_ratio,
        zero_total = as.integer(!positive),
        stringsAsFactors = FALSE,
        check.names = FALSE
    )
}

sensitivity_controlled_feature_rows <- function(
    row, first, masks, scores, sensitivity
) {
    mask_count <- nrow(masks$design)
    index <- rep(seq_len(mask_count), each = 2L)
    mode <- rep(c("global_absence", "sample_missing"), times = mask_count)
    partial_score <- scores$scores[-c(1L, 2L)]
    partial_sum <- colSums(scores$partials, na.rm = TRUE)
    full_sum <- sum(first)
    target <- if (full_sum == 0) rep.int(0, mask_count) else
        abs(partial_score - scores$scores[[1L]]) / full_sum
    size <- length(first)
    effective <- size - masks$design$removed_count
    if (!identical(
        as.integer(vapply(sensitivity, `[[`, integer(1L), "effective_size")),
        as.integer(effective)
    )) stop("Sensitivity controlled feature support is invalid.", call. = FALSE)
    baseline <- sensitivity_controlled_baseline(
        partial_score[index], partial_sum[index], size,
        ifelse(mode == "global_absence", effective[index] / size, 1),
        ifelse(mode == "global_absence", 1, effective[index] / size),
        effective[index]
    )
    sensitivity_rows <- do.call(rbind, Map(
        sensitivity_controlled_sensitivity_facts, sensitivity, partial_sum
    ))
    scenario <- sensitivity_controlled_scenario_row(row)
    scenario <- scenario[rep.int(1L, length(index)), , drop = FALSE]
    mask_rows <- masks$design[index, , drop = FALSE]
    result <- cbind(
        scenario, mask_rows, absence_mode = mode,
        full_score_A = scores$scores[[1L]],
        partial_score_A = partial_score[index], full_sum_A = full_sum,
        target = target[index], baseline,
        absence_sample_missing = as.integer(mode == "sample_missing"),
        sensitivity_rows[index, , drop = FALSE]
    )
    rownames(result) <- NULL
    result
}

sensitivity_controlled_technical_row <- function(
    row, first, second, scores, cell
) {
    first_sum <- sum(first)
    second_sum <- sum(second)
    denominator <- first_sum + second_sum
    target <- if (denominator == 0) 0 else {
        2 * abs(scores$scores[[1L]] - scores$scores[[2L]]) / denominator
    }
    if (!identical(cell$effective_size, as.integer(length(first)))) {
        stop("Sensitivity controlled technical support is invalid.", call. = FALSE)
    }
    cbind(
        sensitivity_controlled_scenario_row(row),
        full_score_A = scores$scores[[1L]], full_score_B = scores$scores[[2L]],
        full_sum_A = first_sum, full_sum_B = second_sum,
        detected_fraction_A = mean(first > 0),
        detected_fraction_B = mean(second > 0), target = target,
        sensitivity_controlled_baseline(
            scores$scores[[1L]], first_sum, length(first), 1, 1, length(first)
        ),
        sensitivity_controlled_sensitivity_facts(cell, first_sum)
    )
}

sensitivity_controlled_validate_encoding <- function(
    first, masks, members, scorer, sensitivity, backend
) {
    checked <- unique(c(1L, length(masks$selected)))
    for (index in checked) {
        removed <- masks$selected[[index]]
        missing <- first
        missing[removed] <- NA_real_
        missing_mat <- matrix(missing, dimnames = list(members, "sample"))
        global_mat <- matrix(
            first[-removed], dimnames = list(members[-removed], "sample")
        )
        missing_score <- scorer(
            missing_mat, list(controlled = members), backend
        )[[1L]]
        global_score <- scorer(
            global_mat, list(controlled = members), backend
        )[[1L]]
        missing_cell <- sensitivity(missing, members)
        global_cell <- sensitivity(first[-removed], members[-removed])
        if (!identical(missing_score, global_score) ||
            !identical(missing_cell, global_cell)) {
            stop("Sensitivity absence encodings diverged.", call. = FALSE)
        }
    }
    invisible(TRUE)
}

sensitivity_controlled_observe_scenario <- function(
    row, registry, scorer, sensitivity, backend, verify_encoding = FALSE
) {
    if (nrow(row) != 1L) {
        stop("Sensitivity controlled scenario row is invalid.", call. = FALSE)
    }
    latent <- sensitivity_controlled_latent(row, registry)
    first <- sensitivity_controlled_measurement(
        row, latent, registry, "measurement_A_seed_base"
    )
    second <- sensitivity_controlled_measurement(
        row, latent, registry, "measurement_B_seed_base"
    )
    masks <- sensitivity_controlled_masks(row, latent, registry)
    scores <- sensitivity_controlled_scores(
        first, second, masks, latent$members, scorer, backend
    )
    cells <- lapply(seq_len(ncol(scores$partials)), function(index) {
        sensitivity(scores$partials[, index], latent$members)
    })
    full_cell <- sensitivity(first, latent$members)
    if (verify_encoding) sensitivity_controlled_validate_encoding(
        first, masks, latent$members, scorer, sensitivity, backend
    )
    list(
        feature = sensitivity_controlled_feature_rows(
            row, first, masks, scores, cells
        ),
        technical = sensitivity_controlled_technical_row(
            row, first, second, scores, full_cell
        )
    )
}

sensitivity_controlled_bind_observations <- function(observed) {
    result <- lapply(c("feature", "technical"), function(kind) {
        value <- do.call(rbind, lapply(observed, `[[`, kind))
        rownames(value) <- NULL
        value
    })
    names(result) <- c("feature", "technical")
    result
}

sensitivity_controlled_feature_fields <- function() {
    c(
        "scenario_id", "fold", "member_count", "archetype", "dynamic_range",
        "log_sd", "expected_total", "dispersion", "dropout_fraction",
        "profile_replicate", "removed_fraction", "mask_mechanism",
        "mask_repeat", "fraction_index", "mechanism_index", "mask_seed",
        "removed_count", "removed_members", "absence_mode", "full_score_A",
        "partial_score_A", "full_sum_A", "target", "declared_size",
        "declared_coverage", "observed_fraction", "observed_sum",
        "effective_size", "score", "log2_declared_size",
        "log1p_observed_sum", "score_over_sum", "zero_total",
        "absence_sample_missing", "largest_member", "largest_absolute_delta",
        "largest_delta", "median_absolute_delta", "sensitivity_delta_status",
        "sensitivity_normalized_status", "largest_absolute_delta_over_sum",
        "largest_delta_over_sum", "median_absolute_delta_over_sum"
    )
}

sensitivity_controlled_technical_fields <- function() {
    c(
        "scenario_id", "fold", "member_count", "archetype", "dynamic_range",
        "log_sd", "expected_total", "dispersion", "dropout_fraction",
        "profile_replicate", "full_score_A", "full_score_B", "full_sum_A",
        "full_sum_B", "detected_fraction_A", "detected_fraction_B", "target",
        "declared_size", "declared_coverage", "observed_fraction",
        "observed_sum", "effective_size", "score", "log2_declared_size",
        "log1p_observed_sum", "score_over_sum", "zero_total", "largest_member",
        "largest_absolute_delta", "largest_delta", "median_absolute_delta",
        "sensitivity_delta_status", "sensitivity_normalized_status",
        "largest_absolute_delta_over_sum", "largest_delta_over_sum",
        "median_absolute_delta_over_sum"
    )
}

sensitivity_controlled_predictors <- function(kind, augmented = FALSE) {
    baseline <- c(
        "log2_declared_size", "declared_coverage", "observed_fraction",
        "log1p_observed_sum", "effective_size", "score_over_sum", "zero_total"
    )
    if (kind == "feature") baseline <- c(baseline, "absence_sample_missing")
    if (!kind %in% c("feature", "technical")) {
        stop("Sensitivity controlled model kind is invalid.", call. = FALSE)
    }
    if (augmented) c(
        baseline, "largest_absolute_delta_over_sum", "largest_delta_over_sum",
        "median_absolute_delta_over_sum"
    ) else baseline
}

sensitivity_controlled_validate_schema <- function(value) {
    valid <- is.list(value) && identical(names(value), c("feature", "technical")) &&
        is.data.frame(value$feature) && is.data.frame(value$technical) &&
        identical(names(value$feature), sensitivity_controlled_feature_fields()) &&
        identical(names(value$technical), sensitivity_controlled_technical_fields()) &&
        !anyDuplicated(names(value$feature)) && !anyDuplicated(names(value$technical))
    if (!valid) {
        stop("Sensitivity controlled observation schema is invalid.", call. = FALSE)
    }
    invisible(value)
}

sensitivity_controlled_expected_scenarios <- function(design) {
    rows <- lapply(seq_len(nrow(design)), function(index) {
        sensitivity_controlled_scenario_row(design[index, , drop = FALSE])
    })
    result <- do.call(rbind, rows)
    rownames(result) <- NULL
    result
}

sensitivity_controlled_validate_scenario_identity <- function(value, design) {
    expected <- sensitivity_controlled_expected_scenarios(design)
    feature_index <- rep(seq_len(nrow(expected)), each = 60L)
    valid <- all(vapply(names(expected), function(field) {
        identical(value$feature[[field]], expected[[field]][feature_index]) &&
            identical(value$technical[[field]], expected[[field]])
    }, logical(1L)))
    if (!valid) {
        stop("Sensitivity controlled scenario identity is invalid.", call. = FALSE)
    }
    invisible(value)
}

sensitivity_controlled_expected_mask_seeds <- function(design, masks, registry) {
    unlist(lapply(seq_len(nrow(design)), function(scenario) {
        vapply(seq_len(nrow(masks)), function(mask) {
            sensitivity_controlled_mask_seed(
                design[scenario, , drop = FALSE],
                masks[mask, , drop = FALSE],
                registry
            )
        }, integer(1L))
    }), use.names = FALSE)
}

sensitivity_controlled_validate_mask_identity <- function(value, design, registry) {
    masks <- sensitivity_mask_design(registry)
    mask_index <- rep(seq_len(nrow(masks)), each = 2L)
    tiled_index <- rep(mask_index, times = nrow(design))
    fields <- c(
        "removed_fraction", "mask_mechanism", "mask_repeat",
        "fraction_index", "mechanism_index"
    )
    valid <- all(vapply(fields, function(field) {
        identical(value$feature[[field]], masks[[field]][tiled_index])
    }, logical(1L)))
    modes <- rep(c("global_absence", "sample_missing"), nrow(masks) * nrow(design))
    seeds <- rep(sensitivity_controlled_expected_mask_seeds(
        design, masks, registry
    ), each = 2L)
    sizes <- rep(as.integer(design$member_count), each = 60L)
    removed <- vapply(seq_len(nrow(value$feature)), function(index) {
        sensitivity_controlled_removed_count(
            sizes[[index]], value$feature$removed_fraction[[index]]
        )
    }, integer(1L))
    selected_count <- lengths(strsplit(
        value$feature$removed_members, ";", fixed = TRUE
    ))
    valid <- valid && identical(value$feature$absence_mode, modes) &&
        identical(value$feature$mask_seed, seeds) &&
        identical(value$feature$removed_count, removed) &&
        identical(selected_count, removed)
    if (!valid) stop("Sensitivity controlled mask identity is invalid.", call. = FALSE)
    invisible(value)
}

sensitivity_controlled_validate_baseline <- function(frame, kind) {
    predictors <- sensitivity_controlled_predictors(kind, augmented = TRUE)
    finite <- c("target", "declared_size", "observed_sum", "score", predictors)
    finite <- unique(finite)
    positive <- frame$observed_sum > 0
    expected_ratio <- numeric(nrow(frame))
    expected_ratio[positive] <- frame$score[positive] / frame$observed_sum[positive]
    valid <- all(vapply(finite, function(field) {
        is.numeric(frame[[field]]) && all(is.finite(frame[[field]]))
    }, logical(1L))) && all(frame$declared_size >= frame$effective_size) &&
        all(frame$declared_coverage > 0 & frame$declared_coverage <= 1) &&
        all(frame$observed_fraction > 0 & frame$observed_fraction <= 1) &&
        identical(frame$log2_declared_size, log2(frame$declared_size)) &&
        identical(frame$log1p_observed_sum, log1p(frame$observed_sum)) &&
        identical(frame$score_over_sum, expected_ratio) &&
        identical(frame$zero_total, as.integer(!positive))
    if (!valid) {
        stop("Sensitivity controlled baseline facts are invalid.", call. = FALSE)
    }
    invisible(frame)
}

sensitivity_controlled_validate_diagnostics <- function(frame) {
    zero <- frame$zero_total == 1L
    expected_status <- ifelse(zero, "zero_total", "ordinary")
    positive <- !zero
    signed_ratio <- frame$largest_delta[positive] /
        frame$observed_sum[positive]
    median_ratio <- frame$median_absolute_delta[positive] /
        frame$observed_sum[positive]
    close <- function(first, second) isTRUE(all.equal(
        first, second, tolerance = 1e-12, check.attributes = FALSE
    ))
    valid <- !anyNA(frame$largest_member) &&
        all(frame$sensitivity_delta_status == "ordinary") &&
        identical(frame$sensitivity_normalized_status, expected_status) &&
        identical(frame$largest_absolute_delta, abs(frame$largest_delta)) &&
        identical(
            frame$largest_absolute_delta_over_sum,
            abs(frame$largest_delta_over_sum)
        ) && all(frame$median_absolute_delta >= 0) &&
        all(frame$largest_absolute_delta_over_sum[zero] == 0) &&
        all(frame$largest_delta_over_sum[zero] == 0) &&
        all(frame$median_absolute_delta_over_sum[zero] == 0) &&
        close(frame$largest_delta_over_sum[positive], signed_ratio) &&
        close(frame$median_absolute_delta_over_sum[positive], median_ratio)
    if (!valid) {
        stop("Sensitivity controlled diagnostic facts are invalid.", call. = FALSE)
    }
    invisible(frame)
}

sensitivity_controlled_validate_targets <- function(value) {
    feature <- value$feature
    expected_feature <- numeric(nrow(feature))
    positive <- feature$full_sum_A > 0
    expected_feature[positive] <- abs(
        feature$partial_score_A[positive] - feature$full_score_A[positive]
    ) / feature$full_sum_A[positive]
    technical <- value$technical
    denominator <- technical$full_sum_A + technical$full_sum_B
    expected_technical <- numeric(nrow(technical))
    positive <- denominator > 0
    expected_technical[positive] <- 2 * abs(
        technical$full_score_A[positive] - technical$full_score_B[positive]
    ) / denominator[positive]
    valid <- identical(feature$target, expected_feature) &&
        identical(technical$target, expected_technical) &&
        identical(feature$score, feature$partial_score_A) &&
        identical(technical$score, technical$full_score_A) &&
        identical(technical$observed_sum, technical$full_sum_A)
    if (!valid) stop("Sensitivity controlled targets are invalid.", call. = FALSE)
    invisible(value)
}

sensitivity_controlled_validate_encoding_pairs <- function(feature) {
    pairs <- matrix(seq_len(nrow(feature)), nrow = 2L)
    pair_fields <- c(
        "partial_score_A", "target", "observed_sum", "effective_size",
        "largest_member", "largest_absolute_delta", "largest_delta",
        "median_absolute_delta", "sensitivity_delta_status",
        "sensitivity_normalized_status", "largest_absolute_delta_over_sum",
        "largest_delta_over_sum", "median_absolute_delta_over_sum"
    )
    pair_equal <- all(vapply(pair_fields, function(field) {
        identical(feature[[field]][pairs[1L, ]], feature[[field]][pairs[2L, ]])
    }, logical(1L)))
    coverage <- feature$declared_coverage[pairs[1L, ]] < 1 &
        feature$declared_coverage[pairs[2L, ]] == 1 &
        feature$observed_fraction[pairs[1L, ]] == 1 &
        feature$observed_fraction[pairs[2L, ]] < 1 &
        feature$absence_sample_missing[pairs[1L, ]] == 0L &
        feature$absence_sample_missing[pairs[2L, ]] == 1L
    if (!pair_equal || !all(coverage)) {
        stop("Sensitivity controlled encoding facts are invalid.", call. = FALSE)
    }
    invisible(feature)
}

sensitivity_controlled_smoke_design <- function(design) {
    queries <- list(
        small = design$member_count == "8" & design$archetype == "monotone" &
            design$dispersion == "0" & design$dropout_fraction == "0" &
            design$profile_replicate == "1",
        missing = design$member_count == "32" & design$archetype == "u_shaped" &
            design$dynamic_range == "64" & design$dispersion == "0.1" &
            design$dropout_fraction == "0.2" & design$profile_replicate == "2",
        wave = design$member_count == "128" & design$archetype == "wave" &
            design$expected_total == "10000" & design$dispersion == "0" &
            design$profile_replicate == "9",
        dominant = design$member_count == "128" & design$archetype == "dominant" &
            design$log_sd == "0.5" & design$dispersion == "0.1" &
            design$dropout_fraction == "0.2" & design$profile_replicate == "10"
    )
    selected <- vapply(queries, function(query) match(TRUE, query, 0L), integer(1L))
    if (any(selected == 0L) || anyDuplicated(selected)) {
        stop("Sensitivity controlled smoke strata are absent.", call. = FALSE)
    }
    design[selected, , drop = FALSE]
}

sensitivity_controlled_validate_observations <- function(value, design, registry) {
    sensitivity_controlled_validate_schema(value)
    sensitivity_controlled_validate_scenario_identity(value, design)
    sensitivity_controlled_validate_mask_identity(value, design, registry)
    sensitivity_controlled_validate_baseline(value$feature, "feature")
    sensitivity_controlled_validate_baseline(value$technical, "technical")
    sensitivity_controlled_validate_diagnostics(value$feature)
    sensitivity_controlled_validate_diagnostics(value$technical)
    sensitivity_controlled_validate_targets(value)
    sensitivity_controlled_validate_encoding_pairs(value$feature)
    invisible(value)
}

sensitivity_controlled_digest <- function(value) {
    path <- tempfile("sensitivity-controlled-", fileext = ".rds")
    on.exit(unlink(path, force = TRUE), add = TRUE)
    saveRDS(value, path, compress = FALSE, version = 3L)
    unname(tools::md5sum(path))
}
