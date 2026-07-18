# Assisted-by: OpenAI Codex.

test_that("scaled component oracle agrees in the ordinary safe region", {
    fixtures <- list(
        c(2, 2),
        c(6, 2),
        rep(4 / 3, 3L),
        c(1, 2, NA_real_),
        c(0, 0.25, 2, 7, 19)
    )

    for (values in fixtures) {
        ordinary <- reference_components(values)
        scaled <- reference_scaled_components(values)

        expect_identical(scaled$semantic, "scoreable")
        expect_identical(scaled$conditioning, "safe")
        expect_equal(
            scaled_ref_as_double(scaled$observed_sum),
            ordinary$observed_sum,
            tolerance = 2e-15
        )
        expect_equal(
            scaled_ref_as_double(scaled$score),
            ordinary$score,
            tolerance = 2e-15
        )
        expect_equal(
            scaled_ref_as_double(scaled$penalty),
            ordinary$penalty,
            tolerance = 2e-15
        )
        expect_equal(
            scaled_ref_as_double(scaled$balance),
            ordinary$balance,
            tolerance = 2e-15
        )

        reconstructed <- scaled_ref_subtract(
            scaled$observed_sum,
            scaled$penalty
        )
        factored <- scaled_ref_multiply(
            scaled$observed_sum,
            scaled$balance
        )
        expect_lt(scaled_ref_relative_error(reconstructed, scaled$score), 2e-29)
        expect_lt(scaled_ref_relative_error(factored, scaled$score), 2e-29)
    }
})

test_that("scaled component oracle preserves overflowed diagnostics", {
    maximum <- .Machine$double.xmax

    sum_overflow <- reference_scaled_components(c(maximum, maximum / 2))
    expect_identical(sum_overflow$semantic, "scoreable")
    expect_identical(sum_overflow$conditioning, "ill_conditioned")
    expect_true(is.na(scaled_ref_as_double(sum_overflow$observed_sum)))
    expect_identical(scaled_ref_as_double(sum_overflow$score), maximum)
    expect_identical(scaled_ref_as_double(sum_overflow$penalty), maximum / 2)
    expect_equal(scaled_ref_as_double(sum_overflow$balance), 2 / 3)
    expect_equal(
        scaled_ref_pair(sum_overflow$observed_sum)[["exponent"]],
        1025
    )

    penalty_overflow <- reference_scaled_components(c(
        maximum,
        maximum,
        0,
        0
    ))
    expect_true(is.na(scaled_ref_as_double(penalty_overflow$observed_sum)))
    expect_true(is.na(scaled_ref_as_double(penalty_overflow$penalty)))
    expect_true(is.finite(scaled_ref_as_double(penalty_overflow$score)))
    expect_equal(scaled_ref_as_double(penalty_overflow$balance), 1 / 3)
    expect_identical(penalty_overflow$conditioning, "ill_conditioned")
})

test_that("scaled oracle exposes cancellation and balance underflow", {
    large <- 2^1023
    small <- 2^-1022
    components <- reference_scaled_components(c(large, small, 0))

    expect_identical(scaled_ref_as_double(components$observed_sum), large)
    expect_identical(scaled_ref_as_double(components$penalty), large)
    expect_equal(
        scaled_ref_as_double(components$score),
        1.5 * small,
        tolerance = 0
    )
    expect_true(is.na(scaled_ref_as_double(components$balance)))
    expect_equal(
        scaled_ref_pair(components$balance),
        c(mantissa = 0.75, exponent = -2044),
        tolerance = 2e-15
    )
    expect_identical(components$conditioning, "ill_conditioned")

    binary64_difference <- scaled_ref_as_double(components$observed_sum) -
        scaled_ref_as_double(components$penalty)
    expect_identical(binary64_difference, 0)
    expect_gt(scaled_ref_as_double(components$score), 0)
})

test_that("scaled oracle retains semantic edge states", {
    too_few <- reference_scaled_components(c(-0, NA_real_, NaN))
    expect_identical(too_few$semantic, "too_few_observed")
    expect_identical(too_few$effective_size, 1L)
    expect_identical(scaled_ref_as_double(too_few$observed_sum), 0)
    expect_identical(too_few$conditioning, "not_applicable")

    zero <- reference_scaled_components(c(-0, 0, NA_real_))
    expect_identical(zero$semantic, "zero_total")
    expect_identical(zero$effective_size, 2L)
    expect_identical(scaled_ref_as_double(zero$score), 0)
    expect_identical(scaled_ref_as_double(zero$penalty), 0)
    expect_null(zero$balance)
    expect_identical(zero$conditioning, "not_applicable")

    one_positive <- reference_scaled_components(c(4, 0))
    expect_identical(one_positive$semantic, "scoreable")
    expect_identical(scaled_ref_as_double(one_positive$score), 0)
    expect_identical(scaled_ref_as_double(one_positive$penalty), 4)
    expect_identical(scaled_ref_as_double(one_positive$balance), 0)
    expect_null(one_positive$kappa)
    expect_identical(one_positive$conditioning, "ill_conditioned")

    subnormal <- reference_scaled_components(rep(2^-1074, 2L))
    expect_equal(
        scaled_ref_pair(subnormal$observed_sum),
        c(mantissa = 0.5, exponent = -1072)
    )
    expect_identical(scaled_ref_as_double(subnormal$score), 2^-1073)
    expect_identical(scaled_ref_as_double(subnormal$balance), 1)
})

test_that("scaled identities survive wide exponent ranges", {
    set.seed(18072026)
    for (index in seq_len(64L)) {
        size <- sample(2:16, 1L)
        exponents <- sample(-1074:1023, size, replace = TRUE)
        mantissas <- 0.5 + runif(size) / 2
        values <- vapply(
            seq_len(size),
            function(position) scaled_ref_pow2(
                mantissas[[position]],
                exponents[[position]]
            ),
            numeric(1L)
        )
        components <- reference_scaled_components(values)
        recomposed_sum <- scaled_ref_add(
            components$score,
            components$penalty
        )
        recomposed_score <- scaled_ref_multiply(
            components$observed_sum,
            components$balance
        )

        expect_true(
            scaled_ref_relative_error(
                recomposed_sum,
                components$observed_sum
            ) < 4e-29,
            info = paste("sum identity fixture", index)
        )
        expect_true(
            scaled_ref_relative_error(
                recomposed_score,
                components$score
            ) < 4e-29,
            info = paste("factorization fixture", index)
        )
        permuted <- reference_scaled_components(sample(values))
        expect_true(
            scaled_ref_relative_error(
                permuted$observed_sum,
                components$observed_sum
            ) < 4e-29,
            info = paste("permutation sum fixture", index)
        )
        expect_true(
            scaled_ref_relative_error(
                permuted$score,
                components$score
            ) < 4e-29,
            info = paste("permutation score fixture", index)
        )
    }
})
