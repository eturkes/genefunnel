# Assisted-by: OpenAI Codex.

sensitivity_controlled_binary_predictors <- function(kind) {
    switch(
        kind,
        feature = c("zero_total", "absence_sample_missing"),
        technical = "zero_total",
        stop("Sensitivity controlled model kind is invalid.", call. = FALSE)
    )
}

sensitivity_controlled_scale_schema <- function(frame, predictors, training, kind) {
    values <- as.matrix(frame[, predictors, drop = FALSE])
    storage.mode(values) <- "double"
    means <- colMeans(values[training, , drop = FALSE])
    sds <- apply(values[training, , drop = FALSE], 2L, stats::sd)
    binary <- predictors %in% sensitivity_controlled_binary_predictors(kind)
    keep <- sds != 0
    if (any(!is.finite(values)) || any(!is.finite(means)) ||
        any(!is.finite(sds))) {
        stop("Sensitivity controlled scaling input is non-finite.", call. = FALSE)
    }
    continuous <- which(!binary & keep)
    values[, continuous] <- sweep(
        sweep(values[, continuous, drop = FALSE], 2L, means[continuous], "-"),
        2L, sds[continuous], "/"
    )
    matrix <- cbind(`(Intercept)` = 1, values[, keep, drop = FALSE])
    if (any(!is.finite(matrix))) {
        stop("Sensitivity controlled scaled matrix is non-finite.", call. = FALSE)
    }
    list(
        matrix = matrix,
        scaling = data.frame(
            term = predictors, training_mean = unname(means),
            training_sd = unname(sds), binary = binary, kept = keep,
            stringsAsFactors = FALSE, check.names = FALSE
        )
    )
}

sensitivity_controlled_fit_model <- function(
    frame, predictors, training, testing, kind, label, fold
) {
    schema <- sensitivity_controlled_scale_schema(
        frame, predictors, training, kind
    )
    fit <- stats::lm.fit(
        schema$matrix[training, , drop = FALSE],
        frame$target[training],
        tol = 1e-7
    )
    coefficients <- fit$coefficients
    aliased <- is.na(coefficients)
    coefficients[aliased] <- 0
    prediction <- drop(
        schema$matrix[testing, , drop = FALSE] %*% coefficients
    )
    if (any(!is.finite(coefficients)) || any(!is.finite(prediction))) {
        stop("Sensitivity controlled QR output is non-finite.", call. = FALSE)
    }
    coefficient_rows <- data.frame(
        target = kind, fold = as.integer(fold), model = label,
        term = names(coefficients), estimate = unname(coefficients),
        aliased = unname(aliased), rank = as.integer(fit$rank),
        stringsAsFactors = FALSE, check.names = FALSE
    )
    scaling <- cbind(
        target = kind, fold = as.integer(fold), model = label, schema$scaling
    )
    list(prediction = prediction, coefficients = coefficient_rows, scaling = scaling)
}

sensitivity_controlled_strata_fields <- function(kind) {
    scenario <- c(
        "member_count", "archetype", "dynamic_range", "log_sd",
        "expected_total", "dispersion", "dropout_fraction", "profile_replicate"
    )
    if (kind == "feature") c(
        scenario, "removed_fraction", "mask_mechanism", "mask_repeat",
        "absence_mode"
    ) else if (kind == "technical") scenario else
        stop("Sensitivity controlled strata kind is invalid.", call. = FALSE)
}

sensitivity_controlled_prediction_identity <- function(frame, testing, kind) {
    fields <- c("scenario_id", "fold", sensitivity_controlled_strata_fields(kind))
    missing <- setdiff(fields, names(frame))
    if (length(missing)) {
        stop("Sensitivity controlled prediction identity is incomplete.", call. = FALSE)
    }
    frame[testing, fields, drop = FALSE]
}

sensitivity_controlled_fold_result <- function(fold, frame, kind) {
    training <- frame$fold != fold
    testing <- frame$fold == fold
    baseline <- sensitivity_controlled_fit_model(
        frame, sensitivity_controlled_predictors(kind), training, testing,
        kind, "baseline", fold
    )
    augmented <- sensitivity_controlled_fit_model(
        frame, sensitivity_controlled_predictors(kind, TRUE), training, testing,
        kind, "augmented", fold
    )
    target <- frame$target[testing]
    baseline_error <- (baseline$prediction - target)^2
    augmented_error <- (augmented$prediction - target)^2
    baseline_rmse <- sqrt(mean(baseline_error))
    augmented_rmse <- sqrt(mean(augmented_error))
    if (!is.finite(baseline_rmse) || baseline_rmse == 0 ||
        !is.finite(augmented_rmse)) {
        stop("Sensitivity controlled fold RMSE is invalid.", call. = FALSE)
    }
    list(
        fold = data.frame(
            target = kind, fold = as.integer(fold), training_rows = sum(training),
            testing_rows = sum(testing), baseline_rmse = baseline_rmse,
            augmented_rmse = augmented_rmse,
            rmse_reduction = (baseline_rmse - augmented_rmse) / baseline_rmse,
            stringsAsFactors = FALSE, check.names = FALSE
        ),
        predictions = cbind(
            row_id = which(testing),
            sensitivity_controlled_prediction_identity(frame, testing, kind),
            target = target, baseline_prediction = baseline$prediction,
            augmented_prediction = augmented$prediction,
            baseline_squared_error = baseline_error,
            augmented_squared_error = augmented_error
        ),
        coefficients = rbind(baseline$coefficients, augmented$coefficients),
        scaling = rbind(baseline$scaling, augmented$scaling)
    )
}

sensitivity_controlled_validate_model_frame <- function(frame, kind) {
    predictors <- sensitivity_controlled_predictors(kind, TRUE)
    required <- c("scenario_id", "fold", "target", predictors)
    missing <- setdiff(required, names(frame))
    finite <- if (length(missing)) FALSE else all(vapply(
        c("target", predictors),
        function(field) is.numeric(frame[[field]]) && all(is.finite(frame[[field]])),
        logical(1L)
    ))
    cluster_folds <- split(frame$fold, frame$scenario_id, drop = TRUE)
    fixed_cluster <- all(vapply(cluster_folds, function(value) {
        length(unique(value)) == 1L
    }, logical(1L)))
    valid <- !length(missing) && nrow(frame) > 0L && finite && fixed_cluster &&
        identical(sort(unique(frame$fold)), seq_len(10L)) &&
        all(frame$scenario_id == as.integer(frame$scenario_id))
    if (!valid) {
        stop("Sensitivity controlled model frame is invalid.", call. = FALSE)
    }
    invisible(frame)
}

sensitivity_controlled_cross_validate <- function(frame, kind) {
    sensitivity_controlled_validate_model_frame(frame, kind)
    results <- lapply(seq_len(10L), function(fold) {
        sensitivity_controlled_fold_result(fold, frame, kind)
    })
    bind <- function(field) {
        value <- do.call(rbind, lapply(results, `[[`, field))
        rownames(value) <- NULL
        value
    }
    folds <- bind("fold")
    predictions <- bind("predictions")
    predictions <- predictions[order(predictions$row_id), , drop = FALSE]
    rownames(predictions) <- NULL
    required <- c(
        folds$baseline_rmse, folds$augmented_rmse, folds$rmse_reduction,
        predictions$baseline_prediction, predictions$augmented_prediction
    )
    if (any(!is.finite(required)) ||
        !identical(predictions$row_id, seq_len(nrow(frame)))) {
        stop("Sensitivity controlled cross-validation is invalid.", call. = FALSE)
    }
    list(
        folds = folds, predictions = predictions,
        coefficients = bind("coefficients"), scaling = bind("scaling")
    )
}

sensitivity_controlled_cluster_table <- function(predictions, fold) {
    selected <- predictions$fold == fold
    frame <- predictions[selected, , drop = FALSE]
    scenario_ids <- unique(frame$scenario_id)
    if (!identical(scenario_ids, sort(scenario_ids))) {
        stop("Sensitivity bootstrap cluster order is invalid.", call. = FALSE)
    }
    groups <- split(
        seq_len(nrow(frame)),
        factor(frame$scenario_id, levels = scenario_ids)
    )
    rows <- Map(function(scenario_id, retained) {
        data.frame(
            scenario_id = scenario_id, rows = length(retained),
            baseline_sse = sum(frame$baseline_squared_error[retained]),
            augmented_sse = sum(frame$augmented_squared_error[retained]),
            stringsAsFactors = FALSE, check.names = FALSE
        )
    }, scenario_ids, groups)
    result <- do.call(rbind, rows)
    rownames(result) <- NULL
    result
}

sensitivity_controlled_bootstrap_once <- function(clusters) {
    reductions <- vapply(clusters, function(frame) {
        sampled <- sample.int(nrow(frame), nrow(frame), replace = TRUE)
        count <- sum(frame$rows[sampled])
        baseline <- sqrt(sum(frame$baseline_sse[sampled]) / count)
        augmented <- sqrt(sum(frame$augmented_sse[sampled]) / count)
        if (!is.finite(baseline) || baseline == 0 || !is.finite(augmented)) {
            stop("Sensitivity bootstrap RMSE is invalid.", call. = FALSE)
        }
        (baseline - augmented) / baseline
    }, numeric(1L))
    stats::median(reductions)
}

sensitivity_controlled_bootstrap <- function(
    cross_validation, registry, kind, repeats = NULL
) {
    locked <- as.integer(sensitivity_registry_value(
        registry, "controlled", "bootstrap_repeats"
    ))
    if (is.null(repeats)) repeats <- locked
    repeats <- as.integer(repeats)
    if (repeats < 1L || repeats > locked) {
        stop("Sensitivity bootstrap repeat count is invalid.", call. = FALSE)
    }
    key <- if (kind == "feature") {
        "feature_bootstrap_seed"
    } else if (kind == "technical") {
        "technical_bootstrap_seed"
    } else stop("Sensitivity bootstrap kind is invalid.", call. = FALSE)
    sensitivity_controlled_set_seed(as.integer(sensitivity_registry_value(
        registry, "controlled", key
    )))
    clusters <- lapply(seq_len(10L), function(fold) {
        sensitivity_controlled_cluster_table(cross_validation$predictions, fold)
    })
    estimates <- numeric(repeats)
    for (replicate in seq_len(repeats)) {
        estimates[[replicate]] <- sensitivity_controlled_bootstrap_once(clusters)
    }
    data.frame(
        target = kind, bootstrap_replicate = seq_len(repeats),
        median_fold_rmse_reduction = estimates,
        stringsAsFactors = FALSE, check.names = FALSE
    )
}

sensitivity_controlled_quantile <- function(values, probability) {
    if (length(values) == 0L || any(!is.finite(values))) {
        stop("Sensitivity controlled quantile input is invalid.", call. = FALSE)
    }
    unname(stats::quantile(
        values, probs = probability, names = FALSE, type = 8L
    ))
}

sensitivity_controlled_validate_endpoint_contract <- function(registry) {
    keys <- c(
        "point_estimator", "bootstrap_tail_probability",
        "feature_gate", "technical_gate"
    )
    expected <- c(
        "median_of_10_heldout_fold_RMSE_reductions", "0.05",
        "median_RMSE_reduction>=0.10;bootstrap95_lower>=0.05",
        "median_RMSE_reduction>=0.10;bootstrap95_lower>=0.05"
    )
    observed <- vapply(keys, function(key) {
        sensitivity_registry_value(registry, "controlled", key)
    }, character(1L))
    quantile_type <- sensitivity_registry_value(registry, "general", "quantile_type")
    if (!identical(unname(observed), expected) || quantile_type != "8") {
        stop("Sensitivity endpoint contract disagrees with runner.", call. = FALSE)
    }
    invisible(registry)
}

sensitivity_controlled_endpoint_rows <- function(
    cross_validation, bootstrap, kind
) {
    point <- stats::median(cross_validation$folds$rmse_reduction)
    lower <- sensitivity_controlled_quantile(
        bootstrap$median_fold_rmse_reduction, 0.05
    )
    data.frame(
        target = kind,
        endpoint = c("median_fold_rmse_reduction", "bootstrap95_lower"),
        estimate = c(point, lower), comparison = c(">=", ">="),
        threshold = c(0.10, 0.05),
        count = c(nrow(cross_validation$folds), nrow(bootstrap)),
        passed = c(point >= 0.10, lower >= 0.05),
        stringsAsFactors = FALSE, check.names = FALSE
    )
}

sensitivity_controlled_endpoints <- function(
    feature_cv, feature_bootstrap, technical_cv, technical_bootstrap, registry
) {
    sensitivity_controlled_validate_endpoint_contract(registry)
    result <- rbind(
        sensitivity_controlled_endpoint_rows(
            feature_cv, feature_bootstrap, "feature"
        ),
        sensitivity_controlled_endpoint_rows(
            technical_cv, technical_bootstrap, "technical"
        )
    )
    rownames(result) <- NULL
    result
}

sensitivity_controlled_strata_levels <- function(registry, kind) {
    fields <- sensitivity_scenario_keys()
    levels <- lapply(fields, sensitivity_registry_levels, registry = registry)
    names(levels) <- fields
    if (kind == "feature") {
        levels$removed_fraction <- sensitivity_registry_levels(
            registry, "removed_fraction"
        )
        levels$mask_mechanism <- sensitivity_registry_levels(
            registry, "mask_mechanism"
        )
        levels$mask_repeat <- as.character(seq_len(as.integer(
            sensitivity_registry_value(registry, "controlled", "mask_repeats")
        )))
        levels$absence_mode <- sensitivity_registry_levels(registry, "absence_mode")
    }
    levels
}

sensitivity_controlled_stratum_row <- function(
    predictions, kind, field, level
) {
    selected <- as.character(predictions[[field]]) == level
    baseline <- sqrt(mean(predictions$baseline_squared_error[selected]))
    augmented <- sqrt(mean(predictions$augmented_squared_error[selected]))
    reduction <- if (is.finite(baseline) && baseline > 0) {
        (baseline - augmented) / baseline
    } else NA_real_
    data.frame(
        target = kind, factor = field, level = level,
        rows = sum(selected), scenarios = length(unique(
            predictions$scenario_id[selected]
        )), baseline_rmse = baseline, augmented_rmse = augmented,
        rmse_reduction = reduction,
        median_target = stats::median(predictions$target[selected]),
        stringsAsFactors = FALSE, check.names = FALSE
    )
}

sensitivity_controlled_strata <- function(cross_validation, registry, kind) {
    levels <- sensitivity_controlled_strata_levels(registry, kind)
    rows <- unlist(lapply(names(levels), function(field) {
        lapply(levels[[field]], function(level) {
            sensitivity_controlled_stratum_row(
                cross_validation$predictions, kind, field, level
            )
        })
    }), recursive = FALSE)
    result <- do.call(rbind, rows)
    rownames(result) <- NULL
    if (any(result$rows == 0L) || any(!is.finite(result$median_target))) {
        stop("Sensitivity controlled strata are incomplete.", call. = FALSE)
    }
    result
}

sensitivity_controlled_curve_row <- function(feature, rows, candidate) {
    data.frame(
        parent_protocol = SENSITIVITY_PROTOCOL_VERSION,
        execution_protocol = SENSITIVITY_CONTROLLED_VERSION,
        source_git_head = candidate,
        removed_fraction = feature$removed_fraction[rows[[1L]]],
        mask_mechanism = feature$mask_mechanism[rows[[1L]]],
        absence_mode = feature$absence_mode[rows[[1L]]],
        rows = length(rows), median_feature_loss = stats::median(feature$target[rows]),
        quantile90_feature_loss = sensitivity_controlled_quantile(
            feature$target[rows], 0.9
        ),
        median_largest_absolute_delta_over_sum = stats::median(
            feature$largest_absolute_delta_over_sum[rows]
        ),
        median_absolute_delta_over_sum = stats::median(
            feature$median_absolute_delta_over_sum[rows]
        ),
        median_declared_coverage = stats::median(feature$declared_coverage[rows]),
        median_observed_fraction = stats::median(feature$observed_fraction[rows]),
        zero_total_rows = sum(feature$zero_total[rows] == 1L),
        study_composition_dependent = TRUE, rescue_endpoint = FALSE,
        stringsAsFactors = FALSE, check.names = FALSE
    )
}

sensitivity_controlled_curves <- function(feature, candidate) {
    fields <- c("removed_fraction", "mask_mechanism", "absence_mode")
    if (any(!fields %in% names(feature)) ||
        !grepl("^[0-9a-f]{40}$", candidate)) {
        stop("Sensitivity controlled curve input is invalid.", call. = FALSE)
    }
    key <- interaction(feature[fields], drop = TRUE, lex.order = TRUE)
    groups <- split(seq_len(nrow(feature)), key)
    result <- do.call(rbind, lapply(groups, function(rows) {
        sensitivity_controlled_curve_row(feature, rows, candidate)
    }))
    mechanism <- c(
        "uniform", "low_abundance", "high_abundance", "low_detection",
        "high_detection"
    )
    order <- order(
        match(result$removed_fraction, c(0.125, 0.25, 0.5)),
        match(result$mask_mechanism, mechanism),
        match(result$absence_mode, c("global_absence", "sample_missing"))
    )
    result <- result[order, , drop = FALSE]
    rownames(result) <- NULL
    result
}

sensitivity_controlled_summary <- function(
    feature_rows, technical_rows, endpoints, feature_bootstrap, technical_bootstrap
) {
    feature_pass <- all(endpoints$passed[endpoints$target == "feature"])
    technical_pass <- all(endpoints$passed[endpoints$target == "technical"])
    data.frame(
        parent_protocol = SENSITIVITY_PROTOCOL_VERSION,
        execution_protocol = SENSITIVITY_CONTROLLED_VERSION,
        feature_rows = nrow(feature_rows), technical_rows = nrow(technical_rows),
        latent_scenarios = length(unique(technical_rows$scenario_id)),
        feature_bootstrap_repeats = nrow(feature_bootstrap),
        technical_bootstrap_repeats = nrow(technical_bootstrap),
        feature_gate_pass = feature_pass,
        technical_gate_pass = technical_pass,
        controlled_gate_pass = feature_pass && technical_pass,
        public_api_permitted = FALSE,
        stringsAsFactors = FALSE, check.names = FALSE
    )
}

sensitivity_controlled_analyze <- function(frames, registry, repeats = NULL) {
    feature_cv <- sensitivity_controlled_cross_validate(frames$feature, "feature")
    technical_cv <- sensitivity_controlled_cross_validate(
        frames$technical, "technical"
    )
    feature_bootstrap <- sensitivity_controlled_bootstrap(
        feature_cv, registry, "feature", repeats
    )
    technical_bootstrap <- sensitivity_controlled_bootstrap(
        technical_cv, registry, "technical", repeats
    )
    endpoints <- sensitivity_controlled_endpoints(
        feature_cv, feature_bootstrap, technical_cv, technical_bootstrap, registry
    )
    list(
        feature_cv = feature_cv, technical_cv = technical_cv,
        bootstrap = rbind(feature_bootstrap, technical_bootstrap),
        endpoints = endpoints,
        strata = rbind(
            sensitivity_controlled_strata(feature_cv, registry, "feature"),
            sensitivity_controlled_strata(technical_cv, registry, "technical")
        ),
        summary = sensitivity_controlled_summary(
            frames$feature, frames$technical, endpoints,
            feature_bootstrap, technical_bootstrap
        )
    )
}

sensitivity_controlled_model_smoke <- function(registry) {
    design <- sensitivity_controlled_design(registry)
    scenarios <- sensitivity_controlled_expected_scenarios(design)
    phase <- scenarios$scenario_id
    base <- data.frame(
        scenarios,
        log2_declared_size = log2(scenarios$member_count),
        declared_coverage = 0.65 + 0.3 * ((phase %% 7L) / 6),
        observed_fraction = 0.6 + 0.35 * ((phase %% 11L) / 10),
        log1p_observed_sum = 3 + (phase %% 13L) / 4,
        effective_size = pmax(3, scenarios$member_count - phase %% 3L),
        score_over_sum = 0.35 + 0.5 * ((phase %% 17L) / 16),
        zero_total = 0L,
        largest_absolute_delta_over_sum = 0.01 + abs(sin(phase / 7)) / 8,
        largest_delta_over_sum = sin(phase / 11) / 9,
        median_absolute_delta_over_sum = 0.005 + abs(cos(phase / 13)) / 12,
        stringsAsFactors = FALSE, check.names = FALSE
    )
    signal <- 0.9 * base$largest_absolute_delta_over_sum -
        0.7 * base$largest_delta_over_sum +
        0.6 * base$median_absolute_delta_over_sum
    base$target <- 0.05 + 0.02 * base$declared_coverage + signal +
        0.003 * sin(phase * 1.7)
    technical <- base
    feature <- base[rep(seq_len(nrow(base)), each = 6L), , drop = FALSE]
    within <- rep(seq_len(6L), times = nrow(base))
    feature$removed_fraction <- rep(c(0.125, 0.25, 0.5), length.out = nrow(feature))
    feature$mask_mechanism <- rep(
        c("uniform", "low_abundance", "high_abundance", "low_detection",
            "high_detection"), length.out = nrow(feature)
    )
    feature$mask_repeat <- 1L + (within %% 2L)
    feature$absence_mode <- ifelse(
        within %% 2L == 0L, "sample_missing", "global_absence"
    )
    feature$absence_sample_missing <- as.integer(
        feature$absence_mode == "sample_missing"
    )
    feature$target <- feature$target + 0.002 * cos(within * 2.3)
    rownames(feature) <- rownames(technical) <- NULL
    list(feature = feature, technical = technical)
}

sensitivity_controlled_validate_model_adversaries <- function(
    frames, feature_cv, registry
) {
    changed <- frames$feature
    held_out <- changed$fold == 1L
    changed$target[held_out] <- changed$target[held_out] + 1
    changed_cv <- sensitivity_controlled_cross_validate(changed, "feature")
    original_fold <- feature_cv$predictions$fold == 1L
    changed_fold <- changed_cv$predictions$fold == 1L
    no_leakage <- identical(
        feature_cv$predictions$baseline_prediction[original_fold],
        changed_cv$predictions$baseline_prediction[changed_fold]
    ) && identical(
        feature_cv$predictions$augmented_prediction[original_fold],
        changed_cv$predictions$augmented_prediction[changed_fold]
    )
    aliased <- frames$feature
    aliased$observed_fraction <- aliased$declared_coverage
    aliased_cv <- sensitivity_controlled_cross_validate(aliased, "feature")
    zero_sd <- feature_cv$scaling$term == "zero_total"
    malformed <- frames$feature
    malformed$score_over_sum[[1L]] <- NA_real_
    rejected <- inherits(try(
        sensitivity_controlled_cross_validate(malformed, "feature"),
        silent = TRUE
    ), "try-error")
    first <- sensitivity_controlled_bootstrap(
        feature_cv, registry, "feature", repeats = 8L
    )
    invisible(stats::runif(9L))
    second <- sensitivity_controlled_bootstrap(
        feature_cv, registry, "feature", repeats = 8L
    )
    valid <- no_leakage && any(aliased_cv$coefficients$aliased) &&
        all(!feature_cv$scaling$kept[zero_sd]) && rejected && identical(first, second)
    if (!valid) stop("Sensitivity controlled model adversary failed.", call. = FALSE)
    invisible(TRUE)
}

sensitivity_controlled_validate_summary_smoke <- function(registry) {
    frames <- sensitivity_controlled_model_smoke(registry)
    analysis <- sensitivity_controlled_analyze(frames, registry, repeats = 8L)
    sensitivity_controlled_validate_model_adversaries(
        frames, analysis$feature_cv, registry
    )
    valid <- nrow(analysis$feature_cv$folds) == 10L &&
        nrow(analysis$technical_cv$folds) == 10L &&
        nrow(analysis$bootstrap) == 16L && nrow(analysis$endpoints) == 4L &&
        all(is.finite(analysis$endpoints$estimate)) &&
        all(analysis$endpoints$passed) && nrow(analysis$strata) == 68L &&
        analysis$summary$feature_rows == 34560L &&
        analysis$summary$technical_rows == 5760L
    if (!valid) stop("Sensitivity controlled summary smoke failed.", call. = FALSE)
    invisible(analysis)
}
