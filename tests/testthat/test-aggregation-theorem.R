# Assisted-by: OpenAI Codex.

test_that("weighted aggregation gap satisfies both independent formulas", {
    set.seed(20260718)

    for (index in seq_len(256L)) {
        unit_count <- sample(2:12, 1L)
        member_count <- sample(2:32, 1L)
        units <- matrix(
            stats::rexp(unit_count * member_count, rate = 0.2),
            nrow = unit_count
        )
        zero_count <- sample.int(length(units), 1L) - 1L
        if (zero_count > 0L) {
            units[sample.int(length(units), zero_count)] <- 0
        }
        weights <- stats::rexp(unit_count)
        if (unit_count > 2L) {
            zero_weights <- sample.int(unit_count, sample(0:2, 1L))
            weights[zero_weights] <- 0
        }
        weights <- weights / sum(weights)
        result <- reference_aggregation_gap(units, weights)
        scale <- max(
            1,
            result$aggregate_score,
            result$weighted_unit_score,
            result$formula_gap
        )
        margin <- 2e-12 * scale

        expect_lte(
            abs(result$score_gap - result$formula_gap),
            margin,
            label = paste("score formula fixture", index)
        )
        expect_lte(
            abs(result$formula_gap - result$cancellation_gap),
            margin,
            label = paste("cancellation formula fixture", index)
        )
        expect_gte(
            result$formula_gap,
            -margin,
            label = paste("non-negative fixture", index)
        )
        expect_lte(
            result$formula_gap,
            result$aggregate_score + margin,
            label = paste("aggregate bound fixture", index)
        )
        if (result$aggregate_score > 0) {
            expect_gte(result$normalized_gap, -2e-12)
            expect_lte(result$normalized_gap, 1 + 2e-12)
        }
    }
})

test_that("equality is no opposing coordinate, not profile proportionality", {
    proportional <- rbind(
        c(0, 1, 3, 2),
        c(0, 2, 6, 4),
        c(0, 0, 0, 0)
    )
    proportional_result <- reference_aggregation_gap(
        proportional,
        c(0.2, 0.5, 0.3)
    )
    expect_equal(proportional_result$formula_gap, 0, tolerance = 1e-15)
    expect_true(proportional_result$no_opposing_deviations)

    nonproportional <- rbind(c(4, 3, 1, 0), c(6, 3, 2, 1))
    nonproportional_result <- reference_aggregation_gap(
        nonproportional,
        c(0.25, 0.75)
    )
    expect_false(isTRUE(all.equal(
        nonproportional[1L, ] / nonproportional[1L, 1L],
        nonproportional[2L, ] / nonproportional[2L, 1L]
    )))
    expect_equal(nonproportional_result$formula_gap, 0, tolerance = 1e-15)
    expect_true(nonproportional_result$no_opposing_deviations)

    complementary <- reference_aggregation_gap(
        rbind(c(4, 0), c(0, 4)),
        c(0.5, 0.5)
    )
    expect_identical(complementary$aggregate, c(2, 2))
    expect_equal(complementary$weighted_unit_score, 0, tolerance = 1e-15)
    expect_equal(complementary$formula_gap, 4, tolerance = 1e-15)
    expect_equal(complementary$normalized_gap, 1, tolerance = 1e-15)
    expect_false(complementary$no_opposing_deviations)
})

test_that("zero weights and invariances preserve aggregation facts", {
    units <- rbind(
        c(0, 1, 4, 8, 2),
        c(7, 3, 0, 1, 5),
        c(NA_real_, NaN, Inf, -1, NA_real_),
        c(2, 6, 1, 0, 3)
    )
    weights <- c(0.2, 0.5, 0, 0.3)
    expected <- reference_aggregation_gap(units, weights)
    without_zero <- reference_aggregation_gap(
        units[c(1L, 2L, 4L), , drop = FALSE],
        weights[c(1L, 2L, 4L)]
    )
    comparable <- setdiff(names(expected), "active")
    expect_equal(
        expected[comparable],
        without_zero[comparable],
        tolerance = 2e-15
    )
    expect_identical(expected$active, c(TRUE, TRUE, FALSE, TRUE))

    scaled <- reference_aggregation_gap(13 * units, weights)
    expect_equal(
        scaled$formula_gap,
        13 * expected$formula_gap,
        tolerance = 2e-14
    )
    expect_equal(
        scaled$normalized_gap,
        expected$normalized_gap,
        tolerance = 2e-15
    )

    unit_order <- c(4L, 1L, 3L, 2L)
    member_order <- c(5L, 2L, 4L, 1L, 3L)
    permuted <- reference_aggregation_gap(
        units[unit_order, member_order],
        weights[unit_order]
    )
    expect_equal(permuted$formula_gap, expected$formula_gap, tolerance = 2e-14)
    expect_equal(
        permuted$normalized_gap,
        expected$normalized_gap,
        tolerance = 2e-15
    )

    active_units <- units[c(1L, 2L, 4L), , drop = FALSE]
    masses <- c(2, 5, 3)
    mass_result <- reference_aggregation_gap(active_units, masses / sum(masses))
    physical_sum <- colSums(sweep(active_units, 1L, masses, "*"))
    physical_gap <- reference_score(physical_sum) - sum(
        masses * apply(active_units, 1L, reference_score)
    )
    expect_equal(
        physical_gap,
        sum(masses) * mass_result$formula_gap,
        tolerance = 2e-13
    )
})

test_that("zero aggregate score makes normalized gap undefined", {
    result <- reference_aggregation_gap(
        matrix(0, nrow = 4L, ncol = 6L),
        c(0.1, 0.2, 0.3, 0.4)
    )
    expect_identical(result$aggregate_score, 0)
    expect_identical(result$weighted_unit_score, 0)
    expect_identical(result$formula_gap, 0)
    expect_true(is.na(result$normalized_gap))
    expect_true(result$no_opposing_deviations)
})
