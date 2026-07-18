# Assisted-by: OpenAI Codex.

aggregation_quantile <- function(values, probability) {
    values <- values[is.finite(values)]
    if (length(values) == 0L) {
        return(NA_real_)
    }
    unname(stats::quantile(
        values,
        probs = probability,
        names = FALSE,
        type = 8L
    ))
}

aggregation_spearman <- function(first, second) {
    retained <- is.finite(first) & is.finite(second)
    first <- first[retained]
    second <- second[retained]
    if (length(first) < 3L || length(unique(first)) < 2L ||
        length(unique(second)) < 2L) {
        return(NA_real_)
    }
    unname(stats::cor(first, second, method = "spearman"))
}

aggregation_observation_pairs <- function(observations) {
    first <- observations[
        observations$measurement_replicate == "A",
        ,
        drop = FALSE
    ]
    second_pool <- observations[
        observations$measurement_replicate == "B",
        ,
        drop = FALSE
    ]
    second <- second_pool[
        match(first$latent_id, second_pool$latent_id),
        ,
        drop = FALSE
    ]
    valid <- nrow(first) * 2L == nrow(observations) &&
        !anyDuplicated(first$latent_id) &&
        !anyDuplicated(second_pool$latent_id) &&
        identical(first$latent_id, second$latent_id) &&
        identical(first$latent_seed, second$latent_seed) &&
        identical(first$fold, second$fold) &&
        identical(first$latent_normalized_gap, second$latent_normalized_gap)
    if (!valid) {
        stop(
            "Synthetic observations do not form exact A/B pairs.",
            call. = FALSE
        )
    }
    list(first = first, second = second)
}

aggregation_model_frame <- function(observations, registry) {
    paired <- aggregation_observation_pairs(observations)
    factor_fields <- aggregation_registry_vector(
        registry,
        "synthetic",
        "design_factor_fields"
    )
    prepare <- function(value) {
        value$target <- value$latent_normalized_gap
        value$aggregate_score <- value$observed_aggregate_score
        value$aggregation_gap <- value$observed_aggregation_gap
        value$normalized_gap <- value$observed_normalized_gap
        for (field in factor_fields) {
            levels <- aggregation_registry_vector(
                registry,
                "synthetic",
                field
            )
            value[[field]] <- factor(as.character(value[[field]]), levels)
        }
        value
    }
    first <- prepare(paired$first)
    second <- prepare(paired$second)
    factor_invalid <- vapply(factor_fields, function(field) {
        anyNA(first[[field]]) || anyNA(second[[field]])
    }, logical(1L))
    if (any(factor_invalid)) {
        stop("Synthetic model factor encoding failed.", call. = FALSE)
    }
    numeric_fields <- c(
        "target", "aggregate_score", "raw_total", "detected_fraction",
        "max_weight", "aggregation_gap", "normalized_gap"
    )
    complete <- first$eligible & second$eligible & Reduce(`&`, lapply(
        numeric_fields,
        function(field) is.finite(first[[field]]) & is.finite(second[[field]])
    ))
    if (!any(complete)) {
        stop("No complete synthetic A/B model pairs remain.", call. = FALSE)
    }
    list(
        first = first[complete, , drop = FALSE],
        second = second[complete, , drop = FALSE],
        factor_fields = factor_fields,
        excluded_pairs = sum(!complete)
    )
}

aggregation_model_matrices <- function(model_frame) {
    baseline_terms <- c(
        "aggregate_score", "raw_total", "detected_fraction", "max_weight",
        model_frame$factor_fields
    )
    baseline_formula <- stats::reformulate(baseline_terms)
    augmented_formula <- stats::update.formula(
        baseline_formula,
        ~ . + aggregation_gap + normalized_gap
    )
    combined <- rbind(model_frame$first, model_frame$second)
    contrasts <- stats::setNames(lapply(
        model_frame$factor_fields,
        function(field) stats::contr.treatment(levels(combined[[field]]))
    ), model_frame$factor_fields)
    baseline <- stats::model.matrix(
        baseline_formula,
        combined,
        contrasts.arg = contrasts
    )
    augmented <- stats::model.matrix(
        augmented_formula,
        combined,
        contrasts.arg = contrasts
    )
    first_rows <- seq_len(nrow(model_frame$first))
    second_rows <- nrow(model_frame$first) + first_rows
    list(
        baseline_first = baseline[first_rows, , drop = FALSE],
        baseline_second = baseline[second_rows, , drop = FALSE],
        augmented_first = augmented[first_rows, , drop = FALSE],
        augmented_second = augmented[second_rows, , drop = FALSE]
    )
}

aggregation_qr_model <- function(matrix, response) {
    fitted <- stats::lm.fit(matrix, response)
    coefficients <- fitted$coefficients
    aliased <- is.na(coefficients)
    coefficients[aliased] <- 0
    if (any(!is.finite(coefficients))) {
        stop(
            "Synthetic QR model produced non-finite coefficients.",
            call. = FALSE
        )
    }
    list(
        coefficients = coefficients,
        aliased = aliased,
        rank = fitted$rank
    )
}

aggregation_model_coefficients <- function(model, fold, label) {
    data.frame(
        fold = as.integer(fold),
        model = label,
        term = names(model$coefficients),
        estimate = unname(model$coefficients),
        aliased = unname(model$aliased),
        rank = as.integer(model$rank),
        stringsAsFactors = FALSE,
        check.names = FALSE
    )
}

aggregation_fit_fold_model <- function(
    training_matrix,
    testing_matrix,
    training_response,
    testing_response
) {
    model <- aggregation_qr_model(training_matrix, training_response)
    prediction <- drop(testing_matrix %*% model$coefficients)
    list(
        model = model,
        prediction = prediction,
        squared_error = (prediction - testing_response)^2,
        rmse = sqrt(mean((prediction - testing_response)^2))
    )
}

aggregation_cross_validation_fold <- function(
    fold,
    model_frame,
    matrices
) {
    training <- model_frame$first$fold != fold
    testing <- model_frame$second$fold == fold
    target <- model_frame$second$target[testing]
    baseline <- aggregation_fit_fold_model(
        matrices$baseline_first[training, , drop = FALSE],
        matrices$baseline_second[testing, , drop = FALSE],
        model_frame$first$target[training],
        target
    )
    augmented <- aggregation_fit_fold_model(
        matrices$augmented_first[training, , drop = FALSE],
        matrices$augmented_second[testing, , drop = FALSE],
        model_frame$first$target[training],
        target
    )
    reduction <- (baseline$rmse - augmented$rmse) / baseline$rmse
    list(
        fold = data.frame(
            fold = fold, training_pairs = sum(training),
            testing_pairs = sum(testing), baseline_rmse = baseline$rmse,
            augmented_rmse = augmented$rmse, rmse_reduction = reduction,
            stringsAsFactors = FALSE, check.names = FALSE
        ),
        predictions = data.frame(
            latent_id = model_frame$second$latent_id[testing], fold = fold,
            target = target, baseline_prediction = baseline$prediction,
            augmented_prediction = augmented$prediction,
            baseline_squared_error = baseline$squared_error,
            augmented_squared_error = augmented$squared_error,
            stringsAsFactors = FALSE, check.names = FALSE
        ),
        coefficients = rbind(
            aggregation_model_coefficients(baseline$model, fold, "baseline"),
            aggregation_model_coefficients(augmented$model, fold, "augmented")
        )
    )
}

aggregation_cross_validate <- function(observations, registry) {
    model_frame <- aggregation_model_frame(observations, registry)
    matrices <- aggregation_model_matrices(model_frame)
    fold_count <- as.integer(aggregation_registry_value(
        registry,
        "synthetic",
        "cv_folds"
    ))
    results <- lapply(seq_len(fold_count), function(fold) {
        aggregation_cross_validation_fold(fold, model_frame, matrices)
    })
    folds <- do.call(rbind, lapply(results, `[[`, "fold"))
    predictions <- do.call(rbind, lapply(results, `[[`, "predictions"))
    coefficients <- do.call(rbind, lapply(results, `[[`, "coefficients"))
    rownames(folds) <- rownames(predictions) <- rownames(coefficients) <- NULL
    required <- c(
        folds$baseline_rmse, folds$augmented_rmse, folds$rmse_reduction,
        predictions$baseline_prediction, predictions$augmented_prediction
    )
    if (any(!is.finite(required))) {
        stop(
            "Synthetic cross-validation produced non-finite output.",
            call. = FALSE
        )
    }
    list(
        folds = folds,
        predictions = predictions,
        coefficients = coefficients,
        excluded_pairs = model_frame$excluded_pairs
    )
}

aggregation_bootstrap_once <- function(predictions, fold_rows) {
    reductions <- vapply(fold_rows, function(rows) {
        sampled <- rows[sample.int(
            length(rows),
            length(rows),
            replace = TRUE
        )]
        baseline <- sqrt(mean(predictions$baseline_squared_error[sampled]))
        augmented <- sqrt(mean(predictions$augmented_squared_error[sampled]))
        (baseline - augmented) / baseline
    }, numeric(1L))
    stats::median(reductions)
}

aggregation_bootstrap_reductions <- function(
    predictions, registry, repeats = NULL
) {
    locked_repeats <- as.integer(aggregation_registry_value(
        registry,
        "synthetic",
        "bootstrap_repeats"
    ))
    if (is.null(repeats)) {
        repeats <- locked_repeats
    }
    repeats <- as.integer(repeats)
    if (repeats < 1L || repeats > locked_repeats) {
        stop("Synthetic bootstrap repeat count is invalid.", call. = FALSE)
    }
    base_seed <- as.integer(aggregation_registry_value(
        registry,
        "general",
        "simulation_seed"
    ))
    offset <- as.integer(aggregation_registry_value(
        registry,
        "synthetic",
        "bootstrap_seed_offset"
    ))
    aggregation_set_seed(base_seed + offset)
    fold_rows <- split(seq_len(nrow(predictions)), predictions$fold)
    estimates <- numeric(repeats)
    for (replicate in seq_len(repeats)) {
        estimates[[replicate]] <- aggregation_bootstrap_once(
            predictions,
            fold_rows
        )
    }
    data.frame(
        bootstrap_replicate = seq_len(repeats),
        median_fold_rmse_reduction = estimates,
        stringsAsFactors = FALSE,
        check.names = FALSE
    )
}

aggregation_endpoint_row <- function(
    gate,
    endpoint,
    stratum,
    estimate,
    comparison,
    threshold,
    count
) {
    passed <- is.finite(estimate) && switch(
        comparison,
        `<=` = estimate <= threshold,
        `>=` = estimate >= threshold,
        stop("Unknown synthetic endpoint comparison.", call. = FALSE)
    )
    data.frame(
        gate = gate,
        endpoint = endpoint,
        stratum = stratum,
        estimate = estimate,
        comparison = comparison,
        threshold = threshold,
        count = as.integer(count),
        passed = passed,
        stringsAsFactors = FALSE,
        check.names = FALSE
    )
}

aggregation_validate_endpoint_contract <- function(registry) {
    keys <- c(
        "curve_gate", "stability_gate", "null_gate", "incremental_gate",
        "quantile_type", "bootstrap_tail_probability"
    )
    expected <- c(
        "median_absolute_error<=0.10;quantile90_absolute_error<=0.25",
        "overall_spearman>=0.80;each_nonnull_archetype_spearman>=0.65",
        paste0(
            "complex_like_no_outlier_dropout0_median_R<=0.05;",
            "complex_like_no_outlier_dropout0_quantile95_R<=0.20"
        ),
        "median_fold_RMSE_reduction>=0.10;bootstrap95_lower>=0.05",
        "8",
        "0.025"
    )
    observed <- vapply(keys, function(key) {
        aggregation_registry_value(registry, "synthetic", key)
    }, character(1L))
    if (!identical(unname(observed), expected)) {
        stop(
            "Synthetic endpoint contract disagrees with the runner.",
            call. = FALSE
        )
    }
    invisible(registry)
}

aggregation_curve_endpoints <- function(observations) {
    curve <- observations$archetype != "complex_like" & observations$eligible &
        is.finite(observations$observed_normalized_gap) &
        is.finite(observations$latent_normalized_gap)
    curve_error <- abs(
        observations$observed_normalized_gap[curve] -
            observations$latent_normalized_gap[curve]
    )
    list(
        aggregation_endpoint_row(
            "curve", "median_absolute_error", "non_complex",
            median(curve_error), "<=", 0.10, length(curve_error)
        ),
        aggregation_endpoint_row(
            "curve", "quantile90_absolute_error", "non_complex",
            aggregation_quantile(curve_error, 0.90), "<=", 0.25,
            length(curve_error)
        )
    )
}

aggregation_stability_estimate <- function(paired, stable, archetype = NULL) {
    selected <- stable
    if (!is.null(archetype)) {
        selected <- selected & paired$first$archetype == archetype
    }
    c(
        estimate = aggregation_spearman(
            paired$first$observed_normalized_gap[selected],
            paired$second$observed_normalized_gap[selected]
        ),
        count = sum(selected)
    )
}

aggregation_stability_endpoints <- function(observations) {
    paired <- aggregation_observation_pairs(observations)
    stable <- paired$first$eligible & paired$second$eligible &
        is.finite(paired$first$observed_normalized_gap) &
        is.finite(paired$second$observed_normalized_gap)
    strata <- c("overall", "cascade_like", "regulatory_like", "mixed_direction")
    estimates <- lapply(strata, function(stratum) {
        aggregation_stability_estimate(
            paired,
            stable,
            if (stratum == "overall") NULL else stratum
        )
    })
    Map(function(stratum, estimate) {
        aggregation_endpoint_row(
            "stability", "spearman", stratum, estimate[["estimate"]],
            ">=", if (stratum == "overall") 0.80 else 0.65,
            estimate[["count"]]
        )
    }, strata, estimates)
}

aggregation_null_endpoints <- function(observations) {
    null <- observations$archetype == "complex_like" &
        observations$outlier_multiplier == 1 &
        observations$dropout_fraction == 0 & observations$eligible &
        is.finite(observations$observed_normalized_gap)
    values <- observations$observed_normalized_gap[null]
    list(
        aggregation_endpoint_row(
            "null", "median_normalized_gap", "complex_no_outlier_dropout0",
            median(values), "<=", 0.05, length(values)
        ),
        aggregation_endpoint_row(
            "null", "quantile95_normalized_gap",
            "complex_no_outlier_dropout0",
            aggregation_quantile(values, 0.95), "<=", 0.20, length(values)
        )
    )
}

aggregation_incremental_endpoints <- function(cross_validation, bootstrap) {
    folds <- cross_validation$folds$rmse_reduction
    bootstrap_values <- bootstrap$median_fold_rmse_reduction
    list(
        aggregation_endpoint_row(
            "incremental", "median_fold_rmse_reduction", "paired_A_to_B",
            median(folds), ">=", 0.10, length(folds)
        ),
        aggregation_endpoint_row(
            "incremental", "bootstrap95_lower", "paired_A_to_B",
            aggregation_quantile(bootstrap_values, 0.025), ">=", 0.05,
            length(bootstrap_values)
        )
    )
}

aggregation_synthetic_endpoints <- function(
    observations,
    cross_validation,
    bootstrap,
    registry
) {
    aggregation_validate_endpoint_contract(registry)
    rows <- c(
        aggregation_curve_endpoints(observations),
        aggregation_stability_endpoints(observations),
        aggregation_null_endpoints(observations),
        aggregation_incremental_endpoints(cross_validation, bootstrap)
    )
    endpoints <- do.call(rbind, rows)
    rownames(endpoints) <- NULL
    endpoints
}

aggregation_stratified_errors <- function(observations, registry) {
    selected <- observations$archetype != "complex_like" &
        observations$eligible &
        is.finite(observations$observed_normalized_gap) &
        is.finite(observations$latent_normalized_gap)
    value <- observations[selected, , drop = FALSE]
    value$absolute_error <- abs(
        value$observed_normalized_gap - value$latent_normalized_gap
    )
    fields <- aggregation_registry_vector(
        registry,
        "synthetic",
        "design_factor_fields"
    )
    rows <- list()
    for (field in fields) {
        levels <- aggregation_registry_vector(registry, "synthetic", field)
        for (level in levels) {
            retained <- as.character(value[[field]]) == level
            errors <- value$absolute_error[retained]
            rows[[length(rows) + 1L]] <- data.frame(
                factor = field,
                level = level,
                count = length(errors),
                median_absolute_error = if (length(errors)) {
                    median(errors)
                } else {
                    NA_real_
                },
                quantile90_absolute_error = aggregation_quantile(errors, 0.90),
                median_observed_gap = if (length(errors)) {
                    median(value$observed_normalized_gap[retained])
                } else {
                    NA_real_
                },
                median_latent_gap = if (length(errors)) {
                    median(value$latent_normalized_gap[retained])
                } else {
                    NA_real_
                },
                stringsAsFactors = FALSE,
                check.names = FALSE
            )
        }
    }
    result <- do.call(rbind, rows)
    rownames(result) <- NULL
    result
}

aggregation_synthetic_summary <- function(
    observations,
    endpoints,
    cross_validation,
    bootstrap
) {
    gate_pass <- tapply(endpoints$passed, endpoints$gate, all)
    required <- c("curve", "stability", "null", "incremental")
    data.frame(
        protocol_version = AGGREGATION_PROTOCOL_VERSION,
        measurement_rows = nrow(observations),
        latent_scenarios = length(unique(observations$latent_id)),
        eligible_measurements = sum(observations$eligible),
        ineligible_measurements = sum(!observations$eligible),
        model_excluded_pairs = cross_validation$excluded_pairs,
        bootstrap_repeats = nrow(bootstrap),
        curve_pass = unname(gate_pass[["curve"]]),
        stability_pass = unname(gate_pass[["stability"]]),
        null_pass = unname(gate_pass[["null"]]),
        incremental_pass = unname(gate_pass[["incremental"]]),
        all_synthetic_gates_pass = all(gate_pass[required]),
        stringsAsFactors = FALSE,
        check.names = FALSE
    )
}

aggregation_summary_smoke_observations <- function(design) {
    latent_integer <- as.integer(sub("^B-S0*", "", design$latent_id))
    replicate_sign <- ifelse(design$measurement_replicate == "A", -1, 1)
    target <- 0.08 + 0.18 * design$complementarity +
        0.04 * (design$archetype == "mixed_direction") +
        0.01 * sin(latent_integer)
    observed <- pmin(0.99, pmax(
        0,
        target + replicate_sign * 0.03 * cos(1.7 * latent_integer)
    ))
    design$eligible <- TRUE
    design$reason <- "eligible"
    design$zero_unit_count <- 0L
    design$raw_total <- design$library_depth * design$unit_count *
        (1 - 0.2 * design$dropout_fraction)
    design$detected_fraction <- pmin(
        1,
        0.4 + log10(design$library_depth) / 10
    )
    design$max_weight <- ifelse(
        design$weight_profile == "dominant",
        0.8,
        ifelse(
            design$weight_profile == "moderate",
            0.75,
            1 / design$unit_count
        )
    )
    design$latent_aggregate_score <- 100 + design$member_count
    design$latent_weighted_unit_score <- 90
    design$latent_aggregation_gap <- target * design$latent_aggregate_score
    design$latent_normalized_gap <- target
    design$observed_aggregate_score <- 100 + design$member_count +
        2 * cos(latent_integer)
    design$observed_weighted_unit_score <- 90
    design$observed_aggregation_gap <- observed *
        design$observed_aggregate_score
    design$observed_normalized_gap <- observed
    design$identity_residual <- 0
    design
}

aggregation_validate_summary_smoke <- function(design, registry) {
    observations <- aggregation_summary_smoke_observations(design)
    cross_validation <- aggregation_cross_validate(observations, registry)
    bootstrap <- aggregation_bootstrap_reductions(
        cross_validation$predictions,
        registry,
        repeats = 8L
    )
    endpoints <- aggregation_synthetic_endpoints(
        observations,
        cross_validation,
        bootstrap,
        registry
    )
    strata <- aggregation_stratified_errors(observations, registry)
    summary <- aggregation_synthetic_summary(
        observations,
        endpoints,
        cross_validation,
        bootstrap
    )
    valid <- nrow(cross_validation$folds) == 10L &&
        nrow(cross_validation$predictions) == 62208L &&
        nrow(bootstrap) == 8L && nrow(endpoints) == 10L &&
        nrow(strata) == 35L && summary$measurement_rows == 124416L
    if (!valid) {
        stop("Synthetic summary smoke dimensions failed.", call. = FALSE)
    }
    invisible(summary)
}

aggregation_synthetic_report_views <- function(endpoints, strata) {
    endpoint_view <- endpoints
    endpoint_view$estimate <- benchmark_format_number(endpoint_view$estimate)
    endpoint_view$threshold <- benchmark_format_number(endpoint_view$threshold)
    endpoint_view$passed <- ifelse(endpoint_view$passed, "yes", "NO")
    strata_view <- strata
    numeric_fields <- names(strata_view)[vapply(
        strata_view,
        is.numeric,
        logical(1L)
    )]
    strata_view[numeric_fields] <- lapply(
        strata_view[numeric_fields],
        benchmark_format_number
    )
    list(endpoints = endpoint_view, strata = strata_view)
}

aggregation_write_synthetic_report <- function(
    path,
    endpoints,
    summary,
    strata,
    metadata,
    workers
) {
    views <- aggregation_synthetic_report_views(endpoints, strata)
    outcome <- if (summary$all_synthetic_gates_pass) {
        "PASS - every frozen synthetic co-primary gate passed."
    } else {
        "FAIL - at least one frozen synthetic co-primary gate failed."
    }
    lines <- c(
        "# GeneFunnel aggregation synthetic validation",
        "",
        paste0("Protocol `", AGGREGATION_PROTOCOL_VERSION, "`. ", outcome),
        "",
        paste0(
            "This controlled result tests recovery under the frozen ",
            "generator; it does not validate a biological claim or alter ",
            "the exact theorem."
        ),
        "",
        "## Co-primary endpoints",
        "",
        benchmark_markdown_table(views$endpoints),
        "",
        "## Factor-stratified descriptive errors",
        "",
        benchmark_markdown_table(views$strata),
        "",
        "## Environment",
        "",
        benchmark_markdown_table(metadata),
        "",
        "## Reproduce",
        "",
        "```sh",
        paste0(
            "R_LIBS_USER=\"$PWD/.agent/R-library\" Rscript --vanilla ",
            "benchmark/run-aggregation-synthetic.R --workers=", workers
        ),
        "```",
        ""
    )
    writeLines(lines, path)
}
