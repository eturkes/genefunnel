# Assisted-by: OpenAI Codex.

catalogue_one_sided_upper <- function(ratios) {
    logs <- log(ratios)
    estimate <- exp(mean(logs))
    standard_error <- stats::sd(logs) / sqrt(length(logs))
    upper <- if (is.finite(standard_error)) {
        exp(
            mean(logs) +
                stats::qt(0.95, df = length(logs) - 1L) * standard_error
        )
    } else {
        estimate
    }
    c(estimate = estimate, upper = upper)
}

catalogue_pair_results <- function(results, mode, batches) {
    list_results <- results[results$method == "list", , drop = FALSE]
    compiled_results <- results[results$method == "compiled", , drop = FALSE]
    keys <- c("scenario_id", "repeat_id")
    paired <- merge(
        list_results,
        compiled_results,
        by = keys,
        suffixes = c("_list", "_compiled"),
        sort = FALSE
    )
    if (nrow(paired) * 2L != nrow(results)) {
        stop("Catalogue benchmark pairs are incomplete.", call. = FALSE)
    }
    primary <- paste0("cumulative_", batches, "_sec")
    paired$primary_ratio <-
        paired[[paste0(primary, "_compiled")]] /
        paired[[paste0(primary, "_list")]]
    for (batch in seq_len(batches)) {
        field <- paste0("cumulative_", batch, "_sec")
        paired[[paste0("cumulative_", batch, "_ratio")]] <-
            paired[[paste0(field, "_compiled")]] /
            paired[[paste0(field, "_list")]]
    }
    paired$digest_equal <-
        paired$output_sha256_list == paired$output_sha256_compiled
    paired$output_facts_equal <-
        paired$output_bytes_list == paired$output_bytes_compiled &
        paired$output_na_cells_list == paired$output_na_cells_compiled &
        paired$output_sum_list == paired$output_sum_compiled
    environment_fields <- c(
        "R_version", "genefunnel_version", "Matrix_version", "Matrix_path",
        "BiocParallel_version", "BiocParallel_path", "Rcpp_version",
        "Rcpp_path", "RcppArmadillo_version", "RcppArmadillo_path",
        "R_platform"
    )
    paired$environment_equal <- Reduce(`&`, lapply(
        environment_fields,
        function(field) {
            paired[[paste0(field, "_list")]] ==
                paired[[paste0(field, "_compiled")]]
        }
    ))
    paired
}

catalogue_scenario_summary <- function(paired, mode, batches) {
    groups <- split(paired, paired$scenario_id, drop = TRUE)
    summary <- do.call(rbind, lapply(groups, function(group) {
        interval <- catalogue_one_sided_upper(group$primary_ratio)
        cumulative_medians <- vapply(seq_len(5L), function(batch) {
            field <- paste0("cumulative_", batch, "_ratio")
            if (field %in% names(group)) stats::median(group[[field]]) else NA_real_
        }, numeric(1))
        list_alloc <- stats::median(group$scoring_alloc_bytes_list)
        compiled_alloc <- stats::median(group$scoring_alloc_bytes_compiled)
        allocation_limit <- list_alloc + max(65536, 0.01 * list_alloc)
        list_process_max <- stats::median(group$process_max_rss_kib_list)
        compiled_process_max <- stats::median(
            group$process_max_rss_kib_compiled
        )
        list_rss <- stats::median(group$scoring_peak_increment_kib_list)
        compiled_rss <- stats::median(
            group$scoring_peak_increment_kib_compiled
        )
        rss_limit <- list_rss + max(4096, 0.10 * list_rss)
        compile_rss_increment <- stats::median(
            group$compile_peak_increment_kib_compiled
        )
        compile_rss_limit <- 16384 +
            512 * group$canonical_memberships_compiled[[1L]] / 1024
        object_bpm <- stats::median(group$object_bytes_per_membership_compiled)
        wire_bpm <- stats::median(group$wire_bytes_per_membership_compiled)
        required <- c(
            interval,
            list_alloc,
            compiled_alloc,
            list_process_max,
            compiled_process_max,
            list_rss, compiled_rss,
            compile_rss_increment,
            object_bpm,
            wire_bpm
        )
        monitor_quality <- isTRUE(
            all(group$scoring_memory_samples_list >= 2L) &&
                all(group$scoring_memory_samples_compiled >= 2L) &&
                all(group$compile_memory_samples_compiled >= 2L) &&
                all(group$scoring_memory_scope_list != "unavailable") &&
                all(group$scoring_memory_scope_compiled != "unavailable") &&
                all(group$compile_memory_scope_compiled != "unavailable")
        )
        gate <- identical(mode, "gate")
        data.frame(
            scenario_id = group$scenario_id[[1L]],
            storage = group$storage_list[[1L]],
            overlap = group$overlap_list[[1L]],
            repeats = nrow(group),
            primary_ratio_geomean = interval[["estimate"]],
            primary_ratio_upper_95 = interval[["upper"]],
            cumulative_1_ratio_median = cumulative_medians[[1L]],
            cumulative_2_ratio_median = cumulative_medians[[2L]],
            cumulative_3_ratio_median = cumulative_medians[[3L]],
            cumulative_4_ratio_median = cumulative_medians[[4L]],
            cumulative_5_ratio_median = cumulative_medians[[5L]],
            list_scoring_alloc_median_bytes = list_alloc,
            compiled_scoring_alloc_median_bytes = compiled_alloc,
            scoring_allocation_limit_bytes = allocation_limit,
            list_whole_process_max_median_kib = list_process_max,
            compiled_whole_process_max_median_kib = compiled_process_max,
            list_scoring_peak_increment_median_kib = list_rss,
            compiled_scoring_peak_increment_median_kib = compiled_rss,
            scoring_peak_increment_limit_kib = rss_limit,
            compile_peak_increment_median_kib = compile_rss_increment,
            compile_peak_increment_limit_kib = compile_rss_limit,
            object_bytes_per_membership = object_bpm,
            wire_bytes_per_membership = wire_bpm,
            digest_pass = isTRUE(all(group$digest_equal)),
            output_facts_pass = isTRUE(all(group$output_facts_equal)),
            environment_pass = isTRUE(all(group$environment_equal)),
            fingerprint_stability_pass =
                all(!is.na(group$content_sha256_compiled)) &&
                all(!is.na(group$wire_sha256_compiled)) &&
                length(unique(group$content_sha256_compiled)) == 1L &&
                length(unique(group$wire_sha256_compiled)) == 1L,
            metric_availability_pass = !gate ||
                (all(is.finite(required)) && monitor_quality),
            timing_pass = !gate ||
                (is.finite(interval[["upper"]]) && interval[["upper"]] <= 0.85),
            amortization_pass = !gate ||
                (is.finite(cumulative_medians[[batches]]) &&
                    cumulative_medians[[batches]] <= 1),
            allocation_pass = !gate ||
                (is.finite(compiled_alloc) && compiled_alloc <= allocation_limit),
            rss_pass = !gate ||
                (is.finite(compiled_rss) && compiled_rss <= rss_limit),
            compile_rss_pass = !gate ||
                (is.finite(compile_rss_increment) &&
                    compile_rss_increment <= compile_rss_limit),
            object_size_pass = !gate ||
                (is.finite(object_bpm) && object_bpm <= 256),
            wire_size_pass = !gate ||
                (is.finite(wire_bpm) && wire_bpm <= 192),
            stringsAsFactors = FALSE,
            check.names = FALSE
        )
    }))
    rownames(summary) <- NULL
    summary
}

catalogue_summary_pass <- function(summary) {
    pass_fields <- grep("_pass$", names(summary), value = TRUE)
    isTRUE(all(vapply(summary[pass_fields], function(value) {
        all(value)
    }, logical(1))))
}

catalogue_write_report <- function(
    path,
    summary,
    metadata,
    correctness,
    options
) {
    pass_fields <- grep("_pass$", names(summary), value = TRUE)
    all_gates <- apply(
        summary[pass_fields],
        1L,
        function(value) isTRUE(all(value))
    )
    view <- data.frame(
        Scenario = summary$scenario_id,
        Storage = summary$storage,
        Overlap = summary$overlap,
        Repeats = summary$repeats,
        `Primary ratio` = benchmark_format_number(summary$primary_ratio_geomean),
        `One-sided 95% upper` = benchmark_format_number(
            summary$primary_ratio_upper_95
        ),
        `Object B/member` = benchmark_format_number(
            summary$object_bytes_per_membership
        ),
        `Wire B/member` = benchmark_format_number(
            summary$wire_bytes_per_membership
        ),
        `Identity pass` = summary$digest_pass & summary$output_facts_pass,
        `Fingerprint pass` = summary$fingerprint_stability_pass,
        `Environment pass` = summary$environment_pass,
        `Timing pass` = summary$timing_pass,
        `Amortization pass` = summary$amortization_pass,
        `Resource pass` = summary$allocation_pass & summary$rss_pass &
            summary$compile_rss_pass & summary$object_size_pass &
            summary$wire_size_pass & summary$metric_availability_pass,
        `All gates` = all_gates,
        check.names = FALSE
    )
    scope <- if (identical(options$mode, "smoke")) {
        "Smoke mode validates orchestration and applies no performance threshold."
    } else {
        "Gate mode applies every frozen correctness, timing, and resource threshold."
    }
    lines <- c(
        "# GeneFunnel compiled-catalogue comparison",
        "",
        paste0(
            "Protocol `", CATALOGUE_PROTOCOL_VERSION, "`; mode `", options$mode,
            "`; list `", options$baseline_id, "`; compiled `",
            options$candidate_id, "`."
        ),
        "",
        paste0(
            "Synthetic computational evidence only; no biological claim. ",
            scope
        ),
        "",
        benchmark_markdown_table(view),
        "",
        "## Correctness",
        "",
        benchmark_markdown_table(correctness),
        "",
        "## Environment",
        "",
        benchmark_markdown_table(metadata),
        ""
    )
    writeLines(lines, path)
}
