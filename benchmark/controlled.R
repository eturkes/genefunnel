# Assisted-by: OpenAI Codex.

benchmark_controlled_scenarios <- c(
    "equal_values",
    "equal_sums_deviation",
    "all_zero",
    "one_nonzero",
    "sample_missingness",
    "partial_coverage",
    "sparse_single_cell",
    "sample_independence",
    "gene_set_independence"
)

benchmark_controlled_format <- function(value) {
    formatted <- vapply(value, function(element) {
        if (is.nan(element)) {
            return("NaN")
        }
        if (is.na(element)) {
            return("NA")
        }
        formatC(element, digits = 17L, format = "fg", flag = "#")
    }, character(1))
    paste(formatted, collapse = ", ")
}

benchmark_controlled_results <- function() {
    backend <- BiocParallel::SerialParam(progressbar = FALSE)
    rows <- list()

    record <- function(
        scenario_id,
        assertion_id,
        metric,
        expected,
        observed,
        tolerance,
        passed
    ) {
        rows[[length(rows) + 1L]] <<- data.frame(
            scenario_id = scenario_id,
            assertion_id = assertion_id,
            metric = metric,
            expected = expected,
            observed = observed,
            tolerance = as.double(tolerance),
            passed = isTRUE(passed),
            stringsAsFactors = FALSE,
            check.names = FALSE
        )
    }

    record_numeric <- function(
        scenario_id,
        assertion_id,
        metric,
        observed,
        expected,
        tolerance = NULL
    ) {
        scale <- max(1, abs(expected), na.rm = TRUE)
        if (!is.finite(scale)) {
            scale <- 1
        }
        if (is.null(tolerance)) {
            tolerance <- 64 * .Machine$double.eps * scale
        }
        passed <- if (length(observed) != length(expected)) {
            FALSE
        } else {
            missing_match <- all(is.na(observed) == is.na(expected))
            nan_match <- all(is.nan(observed) == is.nan(expected))
            finite <- !is.na(expected)
            missing_match && nan_match &&
                all(abs(observed[finite] - expected[finite]) <= tolerance)
        }
        record(
            scenario_id,
            assertion_id,
            metric,
            benchmark_controlled_format(expected),
            benchmark_controlled_format(observed),
            tolerance,
            passed
        )
    }

    record_true <- function(scenario_id, assertion_id, metric, observed, detail) {
        record(
            scenario_id,
            assertion_id,
            metric,
            detail,
            if (isTRUE(observed)) "TRUE" else "FALSE",
            0,
            observed
        )
    }

    record_matrix <- function(scenario_id, assertion_id, observed, expected) {
        same_dimensions <- identical(dim(observed), dim(expected))
        same_names <- identical(dimnames(observed), dimnames(expected))
        same_missingness <- same_dimensions &&
            identical(is.na(observed), is.na(expected)) &&
            identical(is.nan(observed), is.nan(expected))
        finite <- if (same_dimensions && same_missingness) !is.na(expected) else logical()
        differences <- if (length(finite) > 0L) {
            abs(observed[finite] - expected[finite])
        } else {
            Inf
        }
        maximum_difference <- if (length(differences) == 0L) 0 else max(differences)
        scale <- if (same_dimensions && any(finite)) {
            max(1, abs(expected[finite]))
        } else {
            1
        }
        tolerance <- 64 * .Machine$double.eps * scale
        passed <- same_dimensions && same_names && same_missingness &&
            is.finite(maximum_difference) && maximum_difference <= tolerance
        record(
            scenario_id,
            assertion_id,
            "matrix equality",
            "equal dimensions, dimnames, missingness, and finite values",
            paste0(
                "max_abs_difference=", benchmark_controlled_format(maximum_difference),
                "; dimnames_equal=", same_names
            ),
            tolerance,
            passed
        )
    }

    score <- function(mat, gene_sets) {
        genefunnel::genefunnel(mat, gene_sets, BPPARAM = backend)
    }

    equal_mat <- matrix(
        c(4, 4, 4),
        ncol = 1L,
        dimnames = list(c("A", "B", "C"), "sample")
    )
    equal_score <- score(equal_mat, list(equal = c("A", "B", "C")))[1L, 1L]
    record_numeric(
        "equal_values", "score_equals_sum", "score", equal_score, 12
    )

    deviation_mat <- cbind(
        low_deviation = c(A = 2, B = 2, C = 2),
        high_deviation = c(A = 6, B = 0, C = 0)
    )
    deviation_scores <- score(
        deviation_mat,
        list(activity = c("A", "B", "C"))
    )[1L, ]
    record_numeric(
        "equal_sums_deviation",
        "input_sums_match",
        "difference between observed sums",
        diff(colSums(deviation_mat)),
        0
    )
    record_numeric(
        "equal_sums_deviation",
        "low_deviation_score",
        "score",
        deviation_scores[["low_deviation"]],
        6
    )
    record_numeric(
        "equal_sums_deviation",
        "high_deviation_score",
        "score",
        deviation_scores[["high_deviation"]],
        0
    )

    zero_mat <- matrix(
        0,
        nrow = 3L,
        dimnames = list(c("A", "B", "C"), "sample")
    )
    zero_score <- score(zero_mat, list(zero = c("A", "B", "C")))[1L, 1L]
    record_numeric("all_zero", "measured_zeros_score_zero", "score", zero_score, 0)

    one_mat <- matrix(
        c(4, 0, 0),
        ncol = 1L,
        dimnames = list(c("A", "B", "C"), "sample")
    )
    one_score <- score(one_mat, list(one = c("A", "B", "C")))[1L, 1L]
    record_numeric("one_nonzero", "one_nonzero_scores_zero", "score", one_score, 0)

    missing_mat <- rbind(
        A = c(complete = 1, one_missing = 1, insufficient = 4),
        B = c(complete = 2, one_missing = 2, insufficient = NA_real_),
        C = c(complete = 3, one_missing = NA_real_, insufficient = NaN)
    )
    missing_scores <- score(
        missing_mat,
        list(missingness = c("A", "B", "C"))
    )[1L, ]
    record_numeric(
        "sample_missingness",
        "sample_specific_effective_sizes",
        "scores for effective sizes 3, 2, and 1",
        missing_scores,
        c(4.5, 2, NA_real_)
    )

    partial_mat <- matrix(
        c(1, 2),
        ncol = 1L,
        dimnames = list(c("A", "B"), "sample")
    )
    partial_sets <- list(partial = c("A", "B", "C"))
    partial_coverage <- genefunnel::gene_set_coverage(
        partial_sets,
        rownames(partial_mat)
    )
    partial_score <- score(partial_mat, partial_sets)[1L, 1L]
    record_numeric(
        "partial_coverage",
        "matched_members",
        "matched size",
        partial_coverage$matched_size,
        2
    )
    record_numeric(
        "partial_coverage",
        "coverage_fraction",
        "coverage",
        partial_coverage$coverage,
        2 / 3
    )
    record_numeric(
        "partial_coverage",
        "matched_pair_score",
        "score",
        partial_score,
        2
    )

    sparse_mat <- benchmark_sparse_matrix(
        n_features = 2000L,
        n_samples = 80L,
        density = 0.03,
        stored_missing_fraction = 0.02,
        seed = 26071L
    )
    sparse_sets <- benchmark_gene_sets(
        rownames(sparse_mat),
        n_sets = 100L,
        set_size = 16L,
        overlap = "low",
        seed = 26072L
    )
    sparse_scores <- score(sparse_mat, sparse_sets)
    dense_scores <- score(as.matrix(sparse_mat), sparse_sets)
    record_matrix(
        "sparse_single_cell",
        "dense_sparse_equivalence",
        sparse_scores,
        dense_scores
    )
    record_true(
        "sparse_single_cell",
        "sparse_input_is_compact",
        "input object size",
        as.numeric(object.size(sparse_mat)) < 8 * nrow(sparse_mat) * ncol(sparse_mat),
        "sparse object is smaller than its logical dense double payload"
    )

    independence_mat <- rbind(
        A = c(retained = 1, removed = 4),
        B = c(retained = 2, removed = 3),
        C = c(retained = 3, removed = 2),
        D = c(retained = 7, removed = 1)
    )
    target_set <- list(target = c("A", "B", "C"))
    baseline_score <- score(independence_mat, target_set)[1L, "retained"]
    sample_variant <- cbind(
        added = c(A = 9, B = 0, C = 5, D = 2),
        retained = independence_mat[, "retained"]
    )
    variant_score <- score(sample_variant, target_set)[1L, "retained"]
    record_numeric(
        "sample_independence",
        "retained_sample_unchanged",
        "retained sample score",
        variant_score,
        baseline_score
    )

    baseline_sets <- list(
        target = c("A", "B", "C"),
        existing_overlap = c("A", "D")
    )
    baseline_target <- score(
        independence_mat,
        baseline_sets
    )["target", , drop = FALSE]
    expanded_sets <- list(
        existing_overlap = c("C", "D"),
        target = c("A", "B", "C"),
        added_other = c("B", "D")
    )
    expanded_target <- score(
        independence_mat,
        expanded_sets
    )["target", , drop = FALSE]
    record_matrix(
        "gene_set_independence",
        "target_set_unchanged",
        expanded_target,
        baseline_target
    )

    results <- do.call(rbind, rows)
    rownames(results) <- NULL
    missing_scenarios <- setdiff(benchmark_controlled_scenarios, results$scenario_id)
    if (length(missing_scenarios) > 0L) {
        stop(
            "Controlled benchmark omitted scenarios: ",
            paste(missing_scenarios, collapse = ", "),
            call. = FALSE
        )
    }
    results
}
