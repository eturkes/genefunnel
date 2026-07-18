# Assisted-by: OpenAI Codex.

test_that("leave-one-observed-member-out deltas have the fixed sign", {
    values <- c(A = 1, B = 2, C = 7)
    deltas <- sensitivity_deltas_reference(values)
    summary <- sensitivity_summary_reference(values)

    expect_equal(deltas, c(A = 0.5, B = 2.5, C = 2.5))
    expect_identical(summary$effective_size, 3L)
    expect_identical(summary$largest_member, "B")
    expect_equal(summary$largest_absolute_delta, 2.5)
    expect_equal(summary$largest_delta, 2.5)
    expect_equal(summary$largest_delta_over_sum, 0.25)
    expect_equal(summary$median_absolute_delta, 2.5)
})

test_that("absolute ties retain sign and use canonical member order", {
    values <- c(A = 0, B = 1, C = 1)
    deltas <- sensitivity_deltas_reference(values)
    summary <- sensitivity_summary_reference(values)

    expect_equal(deltas, c(A = -1, B = 1, C = 1))
    expect_identical(summary$largest_member, "A")
    expect_equal(summary$largest_absolute_delta, 1)
    expect_equal(summary$largest_delta, -1)
    expect_equal(summary$largest_delta_over_sum, -0.5)
    expect_equal(summary$median_absolute_delta, 1)
})

test_that("zero totals and insufficient support remain distinct", {
    zeros <- sensitivity_summary_reference(c(A = 0, B = 0, C = 0))
    equal <- sensitivity_summary_reference(c(A = 4, B = 4, C = 4))
    pair <- sensitivity_summary_reference(c(A = 1, B = 2, C = NA_real_))

    expect_identical(zeros$largest_member, "A")
    expect_equal(zeros$largest_absolute_delta, 0)
    expect_equal(zeros$largest_delta, 0)
    expect_true(is.na(zeros$largest_delta_over_sum))
    expect_equal(zeros$median_absolute_delta, 0)

    expect_identical(equal$largest_member, "A")
    expect_equal(equal$largest_delta, 4)
    expect_equal(equal$largest_delta_over_sum, 1 / 3)
    expect_equal(equal$median_absolute_delta, 4)

    expect_identical(pair$effective_size, 2L)
    expect_true(is.na(pair$largest_member))
    expect_true(all(vapply(pair[-1L], is.na, logical(1))))
})

test_that("canonical support excludes duplicates unmatched and missing", {
    support <- sensitivity_members_reference(
        c("A", "B", "A", "Z", "C", "D", "E"),
        c("E", "D", "C", "B", "A"),
        c(E = NA_real_, D = 4, C = NaN, B = 2, A = 1)
    )

    expect_identical(support$members, c("A", "B", "D"))
    expect_identical(support$values, c(1, 2, 4))
    expect_named(
        sensitivity_deltas_reference(support$values, support$members),
        c("A", "B", "D")
    )
})

test_that("even-sized summaries use the ordinary midpoint median", {
    summary <- sensitivity_summary_reference(c(A = 1, B = 2, C = 3, D = 9))

    expect_identical(summary$largest_member, "C")
    expect_equal(summary$largest_delta, 3.5)
    expect_equal(summary$largest_delta_over_sum, 7 / 30)
    expect_equal(summary$median_absolute_delta, 2.75)
})

test_that("deltas are homogeneous bounded and shift equivariant", {
    set.seed(20260718)
    for (iteration in seq_len(100L)) {
        size <- sample(3:12, 1L)
        values <- stats::rexp(size, rate = 0.4)
        deltas <- sensitivity_deltas_reference(values, as.character(seq_len(size)))
        total <- sum(values)
        scale <- stats::runif(1L, min = 0.1, max = 20)
        shift <- stats::runif(1L, min = 0, max = 5)

        expect_equal(
            sensitivity_deltas_reference(scale * values),
            scale * unname(deltas),
            tolerance = 1e-12
        )
        expect_equal(
            sensitivity_deltas_reference(values + shift),
            unname(deltas) + shift,
            tolerance = 1e-12
        )
        expect_true(all(deltas >= -(total - values) - 1e-12 * total))
        expect_true(all(deltas <= total + 1e-12 * total))
        expect_true(all(abs(deltas / total) <= 1 + 1e-12))
    }
})

test_that("member permutations move attached deltas without changing values", {
    values <- c(A = 0.5, B = 3, C = 8, D = 1)
    permutation <- c(3L, 1L, 4L, 2L)
    original <- sensitivity_deltas_reference(values)
    permuted <- sensitivity_deltas_reference(values[permutation])

    expect_equal(unname(permuted), unname(original[permutation]))
    expect_equal(max(abs(permuted)), max(abs(original)))
    expect_equal(stats::median(abs(permuted)), stats::median(abs(original)))
})

test_that("deletion deltas are not additive score contributions", {
    values <- c(A = 1, B = 2, C = 7)
    score <- sensitivity_score_reference(values)
    deltas <- sensitivity_deltas_reference(values)

    expect_equal(score, 4.5)
    expect_equal(sum(deltas), 5.5)
    expect_false(isTRUE(all.equal(sum(deltas), score)))
})

test_that("the score is coordinatewise nondecreasing", {
    set.seed(20260719)
    for (iteration in seq_len(200L)) {
        values <- stats::rexp(sample(2:12, 1L), rate = 0.4)
        index <- sample.int(length(values), 1L)
        increment <- stats::runif(1L, min = 0, max = 20)
        increased <- values
        increased[[index]] <- increased[[index]] + increment

        expect_gte(
            sensitivity_score_reference(increased) + 1e-13 * sum(increased),
            sensitivity_score_reference(values)
        )
    }
})

test_that("non-negativity alone does not identify unknown members", {
    observed <- c(A = 1, B = 2, C = 4)
    total <- sum(observed)
    one_threshold <- 4 * max(observed) - total

    expect_equal(
        sensitivity_score_reference(c(observed, U = one_threshold)),
        4 * total / 3
    )
    expect_equal(
        sensitivity_score_reference(c(observed, U = 1000)),
        4 * total / 3,
        tolerance = 1e-12
    )

    for (limit in c(10, 100, 1000)) {
        expect_equal(
            sensitivity_score_reference(c(observed, U = limit, V = limit)),
            limit / 2 + 3 * total / 2,
            tolerance = 1e-12
        )
    }

    for (observed_size in 1:5) {
        observed <- seq_len(observed_size)^2
        total <- sum(observed)
        for (unknown_size in 1:5) {
            complete_size <- observed_size + unknown_size
            limit <- max(
                (complete_size * max(observed) - total) / unknown_size,
                total / observed_size
            ) + 1
            expected <- unknown_size * (unknown_size - 1) /
                (complete_size - 1) * limit +
                (observed_size + 2 * unknown_size - 1) /
                (complete_size - 1) * total
            expect_equal(
                sensitivity_score_reference(c(
                    observed, rep(limit, unknown_size)
                )),
                expected,
                tolerance = 1e-12
            )
        }
    }
})

test_that("finite member limits bound scores and complete-data deltas", {
    observed <- c(A = 1, B = 2, C = 4)
    lower <- c(U = 0, V = 0)
    upper <- c(U = 3, V = 5)
    score_bounds <- sensitivity_score_bounds_reference(observed, lower, upper)
    delta_bounds <- sensitivity_delta_enclosure_reference(
        observed, lower, upper
    )

    expect_equal(score_bounds, c(lower = 3, upper = 11.25))
    set.seed(20260720)
    for (iteration in seq_len(200L)) {
        unknown <- stats::runif(length(lower), min = lower, max = upper)
        names(unknown) <- names(lower)
        values <- c(observed, unknown)
        score <- sensitivity_score_reference(values)
        deltas <- sensitivity_deltas_reference(values)

        expect_gte(score, score_bounds[["lower"]] - 1e-14)
        expect_lte(score, score_bounds[["upper"]] + 1e-14)
        expect_true(all(deltas >= delta_bounds["lower", ] - 1e-14))
        expect_true(all(deltas <= delta_bounds["upper", ] + 1e-14))
    }
})
