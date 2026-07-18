# Assisted-by: OpenAI Codex.

aggregation_named <- function(values, units) {
    stats::setNames(values, units)
}

aggregation_summary_order <- function(summary) {
    ordered <- summary[
        order(summary$group, summary$gene_set),
        ,
        drop = FALSE
    ]
    rownames(ordered) <- NULL
    ordered
}

test_that("aggregation audit returns the locked relational schema", {
    mat <- cbind(
        u1 = c(4, 0, 2),
        u2 = c(0, 4, 2),
        u3 = c(-1, Inf, NA_real_),
        u4 = c(3, 1, 2)
    )
    rownames(mat) <- c("A", "B", "C")
    groups <- aggregation_named(c("g1", "g1", "g1", "g2"), colnames(mat))
    weights <- aggregation_named(c(0.5, 0.5, 0, 2), colnames(mat))
    gene_sets <- list(
        pair = c("A", "B", "A", "absent"),
        flat = c("A", "B", "C")
    )

    result <- .aggregation_audit(mat, gene_sets, groups, weights)

    expect_named(
        result,
        c("summary", "weights", "removed_members", "unit_scores"),
        ignore.order = FALSE
    )
    expect_true(all(vapply(result, is.data.frame, logical(1L))))
    expect_named(
        result$summary,
        c(
            "group", "gene_set", "eligible", "reason", "aggregate_score",
            "weighted_unit_score", "aggregation_gap", "normalized_gap",
            "normalized_status", "identity_residual", "declared_size",
            "matched_size", "retained_size", "active_unit_count",
            "excluded_unit_count", "removed_member_count"
        ),
        ignore.order = FALSE
    )
    expect_named(
        result$weights,
        c("group", "unit", "input_weight", "effective_weight", "active"),
        ignore.order = FALSE
    )
    expect_named(
        result$removed_members,
        c("group", "gene_set", "member", "unit", "reason"),
        ignore.order = FALSE
    )
    expect_named(
        result$unit_scores,
        c("group", "gene_set", "unit", "effective_weight", "score"),
        ignore.order = FALSE
    )
    expect_identical(result$summary$group, c("g1", "g1", "g2", "g2"))
    expect_identical(result$summary$gene_set, rep(c("pair", "flat"), 2L))
    expect_true(all(result$summary$eligible))
    expect_true(all(result$summary$reason == "eligible"))
    expect_identical(
        result$summary$normalized_status,
        rep("defined", 4L)
    )

    expect_identical(result$summary$aggregate_score, c(4, 6, 2, 4.5))
    expect_identical(result$summary$weighted_unit_score, c(0, 3, 2, 4.5))
    expect_identical(result$summary$aggregation_gap, c(4, 3, 0, 0))
    expect_identical(result$summary$normalized_gap, c(1, 0.5, 0, 0))
    expect_true(all(result$summary$identity_residual == 0))
    expect_identical(result$summary$declared_size, c(3L, 3L, 3L, 3L))
    expect_identical(result$summary$matched_size, c(2L, 3L, 2L, 3L))
    expect_identical(result$summary$retained_size, c(2L, 3L, 2L, 3L))
    expect_identical(result$summary$active_unit_count, c(2L, 2L, 1L, 1L))
    expect_identical(result$summary$excluded_unit_count, c(1L, 1L, 0L, 0L))
    expect_identical(result$summary$removed_member_count, c(1L, 0L, 1L, 0L))

    expect_identical(result$weights$unit, colnames(mat))
    expect_identical(result$weights$input_weight, c(0.5, 0.5, 0, 2))
    expect_identical(result$weights$effective_weight, c(0.5, 0.5, 0, 1))
    expect_identical(result$weights$active, c(TRUE, TRUE, FALSE, TRUE))
    expect_identical(
        result$removed_members,
        data.frame(
            group = c("g1", "g2"),
            gene_set = c("pair", "pair"),
            member = c("absent", "absent"),
            unit = c(NA_character_, NA_character_),
            reason = c("unmatched", "unmatched"),
            stringsAsFactors = FALSE
        )
    )
    expect_identical(result$unit_scores$score, c(0, 0, 3, 3, 2, 4.5))
})

test_that("missingness policies reject or recompute one common support", {
    mat <- rbind(
        A = c(u1 = 4, u2 = 0),
        B = c(u1 = 0, u2 = 4),
        C = c(u1 = NA_real_, u2 = NA_real_)
    )
    groups <- aggregation_named(c("g", "g"), colnames(mat))
    weights <- aggregation_named(c(1, 1), colnames(mat))
    gene_sets <- list(triple = c("A", "B", "C"), pair = c("A", "C"))

    rejected <- .aggregation_audit(
        mat,
        gene_sets,
        groups,
        weights,
        missing = "reject"
    )
    expect_false(any(rejected$summary$eligible))
    expect_true(all(rejected$summary$reason == "active_missing_values"))
    expect_identical(rejected$summary$retained_size, c(3L, 2L))
    expect_true(all(is.na(rejected$summary$aggregation_gap)))
    expect_identical(
        rejected$removed_members$member,
        rep("C", 4L)
    )
    expect_identical(rejected$removed_members$unit, rep(c("u1", "u2"), 2L))
    expect_true(all(rejected$removed_members$reason == "missing"))
    expect_identical(nrow(rejected$unit_scores), 0L)

    intersected <- .aggregation_audit(
        mat,
        gene_sets,
        groups,
        weights,
        missing = "intersection"
    )
    expect_identical(intersected$summary$eligible, c(TRUE, FALSE))
    expect_identical(
        intersected$summary$reason,
        c("eligible", "too_few_common_members")
    )
    expect_identical(intersected$summary$retained_size, c(2L, 1L))
    expect_identical(intersected$summary$aggregate_score, c(4, NA_real_))
    expect_identical(intersected$summary$weighted_unit_score, c(0, NA_real_))
    expect_identical(intersected$summary$aggregation_gap, c(4, NA_real_))
    expect_identical(intersected$summary$normalized_gap, c(1, NA_real_))
    expect_identical(intersected$summary$removed_member_count, c(1L, 1L))

    ordered <- .aggregation_audit(
        mat,
        list(ordered = c("C", "absent", "A", "B")),
        groups,
        weights,
        missing = "intersection"
    )
    expect_identical(
        ordered$removed_members,
        data.frame(
            group = rep("g", 3L),
            gene_set = rep("ordered", 3L),
            member = c("C", "C", "absent"),
            unit = c("u1", "u2", NA_character_),
            reason = c("missing", "missing", "unmatched"),
            stringsAsFactors = FALSE
        )
    )
    expect_identical(ordered$summary$removed_member_count, 2L)
})

test_that("all-zero groups ignore every unit value and report no estimand", {
    mat <- rbind(
        A = c(u1 = -1, u2 = Inf),
        B = c(u1 = NA_real_, u2 = NaN)
    )
    groups <- aggregation_named(c("g", "g"), colnames(mat))
    weights <- aggregation_named(c(0, 0), colnames(mat))

    result <- .aggregation_audit(
        mat,
        list(pair = c("A", "B")),
        groups,
        weights
    )
    expect_identical(result$summary$eligible, FALSE)
    expect_identical(result$summary$reason, "no_positive_weight")
    expect_identical(result$summary$retained_size, 2L)
    expect_identical(result$summary$active_unit_count, 0L)
    expect_identical(result$summary$excluded_unit_count, 2L)
    expect_true(all(is.na(result$summary[c(
        "aggregate_score", "weighted_unit_score", "aggregation_gap",
        "normalized_gap", "identity_residual"
    )])))
    expect_true(all(is.na(result$weights$effective_weight)))
    expect_false(any(result$weights$active))
})

test_that("insufficient matched support is reported without scoring", {
    mat <- matrix(
        c(NA_real_, 2),
        nrow = 1L,
        dimnames = list("A", c("u1", "u2"))
    )
    groups <- aggregation_named(c("g", "g"), colnames(mat))
    weights <- aggregation_named(c(1, 1), colnames(mat))
    result <- .aggregation_audit(
        mat,
        list(empty = character(), one = c("A", "absent")),
        groups,
        weights
    )

    expect_identical(
        result$summary$reason,
        rep("too_few_matched_members", 2L)
    )
    expect_identical(result$summary$declared_size, c(0L, 2L))
    expect_identical(result$summary$matched_size, c(0L, 1L))
    expect_identical(result$summary$retained_size, c(0L, 1L))
    expect_identical(result$summary$removed_member_count, c(0L, 1L))
    expect_identical(result$removed_members$member, "absent")
    expect_true(is.na(result$removed_members$unit))
    expect_identical(nrow(result$unit_scores), 0L)
})

test_that(
    "aggregation inputs fail before ambiguous alignment or active values",
    {
    mat <- rbind(A = c(u1 = 1, u2 = 2), B = c(u1 = 2, u2 = 1))
    sets <- list(pair = c("A", "B"))
    groups <- aggregation_named(c("g", "g"), colnames(mat))
    weights <- aggregation_named(c(1, 1), colnames(mat))

    expect_error(
        .aggregation_audit(mat, sets, unname(groups), weights),
        "named exactly"
    )
    expect_error(
        .aggregation_audit(mat, sets, groups, unname(weights)),
        "named exactly"
    )
    expect_error(
        .aggregation_audit(mat, sets, groups, weights, missing = "omit"),
        "missing"
    )
    expect_error(
        .aggregation_audit(
            mat, sets, structure(groups, class = "spoof"), weights
        ),
        "groups"
    )
    expect_error(
        .aggregation_audit(
            mat, sets, groups, structure(weights, class = "spoof")
        ),
        "weights"
    )
    negative <- weights
    negative[[1L]] <- -1
    expect_error(
        .aggregation_audit(mat, sets, groups, negative),
        "non-negative"
    )
    duplicate_units <- mat
    colnames(duplicate_units) <- c("u1", "u1")
    expect_error(
        .aggregation_audit(duplicate_units, sets, groups, weights),
        "duplicated"
    )

    active_infinite <- mat
    active_infinite[[1L]] <- Inf
    expect_error(
        .aggregation_audit(active_infinite, sets, groups, weights),
        "infinite"
    )
    inactive_infinite <- active_infinite
    inactive_weights <- aggregation_named(c(0, 1), colnames(mat))
    expect_no_error(.aggregation_audit(
        inactive_infinite,
        sets,
        groups,
        inactive_weights
    ))
    }
)

test_that("dense and sparse audits agree with the independent oracle", {
    set.seed(20260719)
    mat <- matrix(
        stats::rexp(9L * 7L, rate = 0.3),
        nrow = 9L,
        dimnames = list(paste0("f", 1:9), paste0("u", 1:7))
    )
    mat[sample.int(length(mat), 12L)] <- 0
    groups <- aggregation_named(
        c("g1", "g1", "g1", "g1", "g2", "g2", "g2"),
        colnames(mat)
    )
    weights <- aggregation_named(c(1, 2, 0, 4, 3, 1, 2), colnames(mat))
    gene_sets <- list(
        first = c("f1", "f3", "f5", "f7"),
        second = c("f2", "f4", "f6", "f8", "f9")
    )

    dense <- .aggregation_audit(mat, gene_sets, groups, weights)
    sparse <- .aggregation_audit(
        Matrix::Matrix(mat, sparse = TRUE),
        gene_sets,
        groups,
        weights
    )
    expect_identical(sparse, dense)

    for (row in seq_len(nrow(dense$summary))) {
        group <- dense$summary$group[[row]]
        set <- dense$summary$gene_set[[row]]
        units <- groups == group & weights > 0
        effective <- weights[units] / sum(weights[units])
        expected <- reference_aggregation_gap(
            t(mat[gene_sets[[set]], units, drop = FALSE]),
            unname(effective)
        )
        scale <- max(1, expected$aggregate_score)
        expect_equal(
            dense$summary$aggregate_score[[row]],
            expected$aggregate_score,
            tolerance = 2e-14 * scale
        )
        expect_equal(
            dense$summary$weighted_unit_score[[row]],
            expected$weighted_unit_score,
            tolerance = 2e-14 * scale
        )
        expect_equal(
            dense$summary$aggregation_gap[[row]],
            expected$formula_gap,
            tolerance = 2e-14 * scale
        )
        expect_equal(
            dense$summary$normalized_gap[[row]],
            expected$normalized_gap,
            tolerance = 2e-14
        )
    }

    unit_order <- c(3L, 1L, 4L, 2L, 7L, 5L, 6L)
    feature_order <- c(9L, 1L, 7L, 3L, 5L, 2L, 8L, 4L, 6L)
    permuted <- .aggregation_audit(
        mat[feature_order, unit_order],
        gene_sets,
        groups[unit_order],
        weights[unit_order]
    )
    expect_equal(
        aggregation_summary_order(permuted$summary),
        aggregation_summary_order(dense$summary),
        tolerance = 2e-14
    )
})

test_that("audit equality uses coordinate signs rather than proportionality", {
    mat <- cbind(
        p1 = c(0, 1, 3, 2),
        p2 = c(0, 2, 6, 4),
        n1 = c(4, 3, 1, 0),
        n2 = c(6, 3, 2, 1)
    )
    rownames(mat) <- LETTERS[1:4]
    groups <- aggregation_named(
        c("proportional", "proportional", "signs", "signs"),
        colnames(mat)
    )
    weights <- aggregation_named(c(2, 5, 1, 3), colnames(mat))

    result <- .aggregation_audit(
        mat,
        list(profile = rownames(mat)),
        groups,
        weights
    )
    expect_true(all(result$summary$eligible))
    expect_equal(result$summary$aggregation_gap, c(0, 0), tolerance = 1e-15)
    expect_equal(result$summary$identity_residual, c(0, 0), tolerance = 1e-14)
})

test_that("randomized group audits match independent weighted references", {
    set.seed(20260720)

    for (fixture in seq_len(128L)) {
        feature_count <- sample(2:16, 1L)
        unit_count <- sample(2:8, 1L)
        group_count <- sample(seq_len(min(3L, unit_count)), 1L)
        features <- paste0("f", seq_len(feature_count))
        units <- paste0("u", seq_len(unit_count))
        mat <- matrix(
            stats::rexp(feature_count * unit_count, rate = 0.4),
            nrow = feature_count,
            dimnames = list(features, units)
        )
        mat[sample.int(length(mat), sample.int(length(mat), 1L) - 1L)] <- 0
        raw_groups <- sample(rep(
            paste0("g", seq_len(group_count)),
            length.out = unit_count
        ))
        groups <- aggregation_named(raw_groups, units)
        raw_weights <- stats::rexp(unit_count)
        raw_weights[stats::runif(unit_count) < 0.2] <- 0
        for (group in unique(raw_groups)) {
            selected <- which(raw_groups == group)
            if (!any(raw_weights[selected] > 0)) {
                raw_weights[selected[[1L]]] <- 1
            }
        }
        weights <- aggregation_named(raw_weights, units)
        set_count <- sample(1:4, 1L)
        gene_sets <- lapply(seq_len(set_count), function(set) {
            size <- sample.int(feature_count - 1L, 1L) + 1L
            members <- sample(features, size)
            if (set %% 2L == 0L) c(members, paste0("absent_", set)) else members
        })
        names(gene_sets) <- paste0("set", seq_len(set_count))

        observed <- .aggregation_audit(mat, gene_sets, groups, weights)
        expect_true(
            all(observed$summary$eligible),
            info = paste("fixture", fixture)
        )
        for (row in seq_len(nrow(observed$summary))) {
            group <- observed$summary$group[[row]]
            set <- observed$summary$gene_set[[row]]
            active <- groups == group & weights > 0
            effective <- weights[active] / sum(weights[active])
            members <- intersect(unique(gene_sets[[set]]), features)
            expected <- reference_aggregation_gap(
                t(mat[members, active, drop = FALSE]),
                unname(effective)
            )
            scale <- max(1, expected$aggregate_score)
            expect_equal(
                observed$summary$aggregate_score[[row]],
                expected$aggregate_score,
                tolerance = 4e-13 * scale,
                info = paste("aggregate fixture", fixture, row)
            )
            expect_equal(
                observed$summary$weighted_unit_score[[row]],
                expected$weighted_unit_score,
                tolerance = 4e-13 * scale,
                info = paste("weighted fixture", fixture, row)
            )
            expect_equal(
                observed$summary$aggregation_gap[[row]],
                expected$formula_gap,
                tolerance = 4e-13 * scale,
                info = paste("gap fixture", fixture, row)
            )
        }
    }
})

test_that("audit preserves mass, scale, physical-sum, and order identities", {
    mat <- rbind(
        A = c(u1 = 0, u2 = 7, u3 = 2, u4 = 4),
        B = c(u1 = 3, u2 = 1, u3 = 8, u4 = 0),
        C = c(u1 = 9, u2 = 2, u3 = 1, u4 = 5),
        D = c(u1 = 4, u2 = 6, u3 = 0, u4 = 3),
        E = c(u1 = 2, u2 = 0, u3 = 7, u4 = 1)
    )
    groups <- aggregation_named(
        c("g1", "g1", "g2", "g2"),
        colnames(mat)
    )
    weights <- aggregation_named(c(2, 5, 3, 1), colnames(mat))
    gene_sets <- list(
        alpha = c("A", "C", "E"),
        beta = c("B", "D", "E", "A")
    )
    baseline <- .aggregation_audit(mat, gene_sets, groups, weights)

    scaled <- .aggregation_audit(13 * mat, gene_sets, groups, weights)
    scale_columns <- c(
        "aggregate_score", "weighted_unit_score", "aggregation_gap",
        "identity_residual"
    )
    expect_equal(
        scaled$summary[scale_columns],
        13 * baseline$summary[scale_columns],
        tolerance = 2e-13
    )
    expect_equal(
        scaled$summary$normalized_gap,
        baseline$summary$normalized_gap,
        tolerance = 2e-14
    )

    rescaled_weights <- weights
    rescaled_weights[groups == "g1"] <- 8 * rescaled_weights[groups == "g1"]
    rescaled_weights[groups == "g2"] <- 16 * rescaled_weights[groups == "g2"]
    mass_invariant <- .aggregation_audit(
        mat,
        gene_sets,
        groups,
        rescaled_weights
    )
    expect_identical(
        mass_invariant$weights$effective_weight,
        baseline$weights$effective_weight
    )
    expect_equal(mass_invariant$summary, baseline$summary, tolerance = 2e-14)

    for (row in seq_len(nrow(baseline$summary))) {
        group <- baseline$summary$group[[row]]
        set <- baseline$summary$gene_set[[row]]
        selected <- groups == group & weights > 0
        masses <- unname(weights[selected])
        physical_profile <- rowSums(sweep(
            mat[gene_sets[[set]], selected, drop = FALSE],
            2L,
            masses,
            "*"
        ))
        physical_score <- reference_score(physical_profile)
        unit_rows <- baseline$unit_scores$group == group &
            baseline$unit_scores$gene_set == set
        unit_scores <- baseline$unit_scores$score[unit_rows]
        total_mass <- sum(masses)
        expect_equal(
            physical_score,
            total_mass * baseline$summary$aggregate_score[[row]],
            tolerance = 2e-13
        )
        expect_equal(
            physical_score - sum(masses * unit_scores),
            total_mass * baseline$summary$aggregation_gap[[row]],
            tolerance = 2e-13
        )
    }

    unit_order <- c(3L, 4L, 2L, 1L)
    feature_order <- c(5L, 3L, 1L, 4L, 2L)
    permuted_sets <- rev(lapply(gene_sets, rev))
    permuted <- .aggregation_audit(
        mat[feature_order, unit_order],
        permuted_sets,
        groups[unit_order],
        weights[unit_order]
    )
    expect_identical(permuted$summary$group, rep(c("g2", "g1"), each = 2L))
    expect_identical(
        permuted$summary$gene_set,
        rep(c("beta", "alpha"), 2L)
    )
    expect_equal(
        aggregation_summary_order(permuted$summary),
        aggregation_summary_order(baseline$summary),
        tolerance = 2e-13
    )
})

test_that("aggregation numerics retain extremes or fail closed", {
    maximum <- .Machine$double.xmax
    mat <- rbind(
        A = c(u1 = maximum, u2 = 0),
        B = c(u1 = 0, u2 = maximum)
    )
    groups <- aggregation_named(c("g", "g"), colnames(mat))
    weights <- aggregation_named(c(1, 1), colnames(mat))
    extreme <- .aggregation_audit(
        mat,
        list(pair = c("A", "B")),
        groups,
        weights
    )
    expect_true(extreme$summary$eligible)
    expect_identical(extreme$summary$aggregate_score, maximum)
    expect_identical(extreme$summary$aggregation_gap, maximum)
    expect_identical(extreme$summary$normalized_gap, 1)

    huge_weights <- aggregation_named(c(maximum, maximum), colnames(mat))
    huge <- .aggregation_audit(
        mat,
        list(pair = c("A", "B")),
        groups,
        huge_weights
    )
    expect_identical(huge$weights$effective_weight, c(0.5, 0.5))
    expect_identical(huge$summary$aggregation_gap, maximum)

    zero <- .aggregation_audit(
        matrix(
            0,
            nrow = 2L,
            ncol = 2L,
            dimnames = list(c("A", "B"), c("u1", "u2"))
        ),
        list(pair = c("A", "B")),
        groups,
        weights
    )
    expect_true(zero$summary$eligible)
    expect_identical(zero$summary$aggregate_score, 0)
    expect_identical(zero$summary$weighted_unit_score, 0)
    expect_identical(zero$summary$aggregation_gap, 0)
    expect_true(is.na(zero$summary$normalized_gap))
    expect_identical(zero$summary$normalized_status, "zero_aggregate")

    smallest <- 2^-1074
    underflow <- .aggregation_audit(
        rbind(A = c(u1 = smallest, u2 = 0), B = c(u1 = 0, u2 = smallest)),
        list(pair = c("A", "B")),
        groups,
        weights
    )
    expect_false(underflow$summary$eligible)
    expect_identical(underflow$summary$reason, "numerically_unavailable")
    expect_true(is.na(underflow$summary$aggregation_gap))

    expect_error(
        .aggregation_audit(
            mat,
            list(pair = c("A", "B")),
            groups,
            aggregation_named(c(smallest, maximum), colnames(mat))
        ),
        "underflow"
    )
})
