# Assisted-by: OpenAI Codex.

test_that("component semantics distinguish missingness and zero total", {
    none_observed <- reference_components(c(NA_real_, NaN))
    expect_identical(none_observed$effective_size, 0L)
    expect_identical(none_observed$observed_sum, 0)
    expect_identical(none_observed$observed_fraction, 0)
    expect_true(all(is.na(unlist(none_observed[c(
        "score", "penalty", "balance"
    )]))))

    one_observed <- reference_components(c(4, NA_real_))
    expect_identical(one_observed$effective_size, 1L)
    expect_identical(one_observed$observed_sum, 4)
    expect_identical(one_observed$observed_fraction, 0.5)
    expect_true(all(is.na(unlist(one_observed[c(
        "score", "penalty", "balance"
    )]))))

    zero_total <- reference_components(c(0, 0, NA_real_))
    expect_identical(zero_total$effective_size, 2L)
    expect_identical(zero_total$observed_sum, 0)
    expect_equal(zero_total$observed_fraction, 2 / 3)
    expect_identical(zero_total$score, 0)
    expect_identical(zero_total$penalty, 0)
    expect_true(is.na(zero_total$balance))

    partial <- reference_components(c(1, 2, NA_real_))
    expect_identical(partial$effective_size, 2L)
    expect_identical(partial$observed_sum, 3)
    expect_equal(partial$observed_fraction, 2 / 3)
    expect_identical(partial$score, 2)
    expect_identical(partial$penalty, 1)
    expect_equal(partial$balance, 2 / 3)
})

test_that("total variation gives the GeneFunnel factorization", {
    set.seed(20260717)
    fixtures <- lapply(seq_len(96L), function(index) {
        size <- sample(2:24, 1L)
        values <- runif(size, min = 0, max = 100)
        zero_count <- sample.int(size, 1L) - 1L
        if (zero_count > 0L) {
            values[sample.int(size, zero_count)] <- 0
        }
        values
    })

    for (index in seq_along(fixtures)) {
        values <- fixtures[[index]]
        components <- reference_components(values)
        shares <- values / sum(values)
        uniform <- 1 / length(values)
        total_variation <- 0.5 * sum(abs(shares - uniform))
        pietra <- mean(abs(values - mean(values))) / (2 * mean(values))
        above_uniform_mass <- sum(shares[shares > uniform] - uniform)
        overlap <- sum(pmin(shares, uniform))
        bulla_evenness <- (overlap - uniform) / (1 - uniform)

        expect_equal(
            reference_score(values),
            components$score,
            tolerance = 2e-13,
            info = paste("factorization fixture", index)
        )
        expect_equal(
            components$penalty,
            components$observed_sum - components$score,
            tolerance = 2e-13,
            info = paste("penalty fixture", index)
        )
        expect_equal(pietra, total_variation, tolerance = 2e-15)
        expect_equal(above_uniform_mass, total_variation, tolerance = 2e-15)
        expect_equal(bulla_evenness, components$balance, tolerance = 2e-15)
        expect_equal(
            components$balance,
            1 - length(values) / (length(values) - 1L) * pietra,
            tolerance = 2e-15
        )
        expect_gte(components$balance, -2e-15)
        expect_lte(components$balance, 1 + 2e-15)
        expect_gte(components$score, -2e-12)
        expect_lte(
            components$score,
            components$observed_sum + 2e-12
        )
    }
})

test_that("components are symmetric and positively homogeneous", {
    values <- c(0, 0.25, 2, 7, 19)
    expected <- reference_components(values)

    for (index in seq_len(32L)) {
        permuted <- reference_components(sample(values))
        expect_equal(permuted, expected, tolerance = 2e-15)
    }

    for (scale in c(2^-20, 0.25, 3, 2^20)) {
        scaled <- reference_components(scale * values)
        expect_equal(scaled$observed_sum, scale * expected$observed_sum)
        expect_equal(scaled$score, scale * expected$score, tolerance = 2e-14)
        expect_equal(
            scaled$penalty,
            scale * expected$penalty,
            tolerance = 2e-14
        )
        expect_equal(scaled$balance, expected$balance, tolerance = 2e-15)
        expect_identical(scaled$effective_size, expected$effective_size)
        expect_identical(scaled$observed_fraction, expected$observed_fraction)
    }

    collapsed <- reference_components(0 * values)
    expect_identical(collapsed$score, 0)
    expect_identical(collapsed$penalty, 0)
    expect_true(is.na(collapsed$balance))
})

test_that("the GeneFunnel score is concave on common support", {
    set.seed(17202607)

    for (index in seq_len(128L)) {
        size <- sample(2:32, 1L)
        left <- rexp(size)
        right <- rexp(size)
        weight <- runif(1L)
        mixture <- weight * left + (1 - weight) * right

        mixture_score <- reference_score(mixture)
        weighted_score <- weight * reference_score(left) +
            (1 - weight) * reference_score(right)
        margin <- 5e-14 * max(1, mixture_score, weighted_score)
        expect_true(
            mixture_score + margin >= weighted_score,
            info = paste("concavity fixture", index)
        )
    }
})

test_that("committed cases expose information hidden by the scalar score", {
    equal_pair <- reference_components(c(2, 2))
    unequal_pair <- reference_components(c(6, 2))
    equal_triple <- reference_components(rep(4 / 3, 3L))
    one_positive <- reference_components(c(4, 0))

    expect_equal(
        vapply(
            list(equal_pair, unequal_pair, equal_triple),
            `[[`,
            numeric(1L),
            "score"
        ),
        rep(4, 3L),
        tolerance = 2e-15
    )
    expect_identical(equal_pair$observed_sum, 4)
    expect_identical(equal_pair$balance, 1)
    expect_identical(equal_pair$effective_size, 2L)
    expect_identical(unequal_pair$observed_sum, 8)
    expect_identical(unequal_pair$balance, 0.5)
    expect_identical(unequal_pair$effective_size, 2L)
    expect_equal(equal_triple$observed_sum, 4, tolerance = 2e-15)
    expect_identical(equal_triple$balance, 1)
    expect_identical(equal_triple$effective_size, 3L)
    expect_identical(one_positive$score, 0)
    expect_identical(one_positive$penalty, 4)
    expect_identical(one_positive$balance, 0)
})
