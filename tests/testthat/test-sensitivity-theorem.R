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
