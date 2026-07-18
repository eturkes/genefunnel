# Assisted-by: OpenAI Codex.

sensitivity_big_as_double <- function(value) {
    base <- genefunnel:::.GF_BIG_BASE
    sum(value * base^(seq_along(value) - 1L))
}

sensitivity_integer_score_numerator <- function(values) {
    size <- length(values)
    total <- sum(values)
    below <- size * values < total
    (size - 1L - sum(below)) * total + size * sum(values[below])
}

sensitivity_integer_delta_numerators <- function(values) {
    size <- length(values)
    full <- sensitivity_integer_score_numerator(values)
    vapply(seq_along(values), function(index) {
        deleted <- sensitivity_integer_score_numerator(values[-index])
        full * (size - 2L) - deleted * (size - 1L)
    }, numeric(1L))
}

sensitivity_exact_score_double <- function(values) {
    dyadic <- genefunnel:::.gf_dyadic_vector(values)
    score <- genefunnel:::.gf_exact_score_numerator(dyadic$magnitudes)
    genefunnel:::.gf_big_ratio_double(
        score$numerator,
        genefunnel:::.gf_big_from_integer(length(values) - 1L),
        dyadic$exponent,
        1L,
        length(values)
    )
}

test_that("exact limbs implement integer arithmetic", {
    base <- genefunnel:::.GF_BIG_BASE
    set.seed(20260719)

    for (iteration in seq_len(200L)) {
        left <- floor(stats::runif(1L, 0, base))
        right <- floor(stats::runif(1L, 0, base))
        left_big <- genefunnel:::.gf_big_from_integer(left)
        right_big <- genefunnel:::.gf_big_from_integer(right)

        expect_identical(
            sensitivity_big_as_double(
                genefunnel:::.gf_big_add(left_big, right_big)
            ),
            left + right,
            info = iteration
        )
        expect_identical(
            sensitivity_big_as_double(
                genefunnel:::.gf_big_multiply(left_big, right_big)
            ),
            left * right,
            info = iteration
        )
        expect_identical(
            sensitivity_big_as_double(genefunnel:::.gf_big_subtract(
                if (left >= right) left_big else right_big,
                if (left >= right) right_big else left_big
            )),
            abs(left - right),
            info = iteration
        )
        expect_identical(
            genefunnel:::.gf_big_compare(left_big, right_big),
            as.integer(sign(left - right)),
            info = iteration
        )
    }

    random_big <- function() {
        genefunnel:::.gf_big_trim(floor(stats::runif(5L, 0, base)))
    }
    for (iteration in seq_len(50L)) {
        left <- random_big()
        right <- random_big()
        factor <- random_big()
        distributed <- genefunnel:::.gf_big_add(
            genefunnel:::.gf_big_multiply(left, factor),
            genefunnel:::.gf_big_multiply(right, factor)
        )
        expect_identical(
            genefunnel:::.gf_big_multiply(
                genefunnel:::.gf_big_add(left, right),
                factor
            ),
            distributed,
            info = iteration
        )
        expect_identical(
            genefunnel:::.gf_big_multiply(left, right),
            genefunnel:::.gf_big_multiply(right, left),
            info = iteration
        )
    }
})

test_that("dyadic decomposition round-trips the binary64 domain", {
    set.seed(20260720)
    exponents <- sample(-1074:1023, 1000L, replace = TRUE)
    randomized <- (1 + stats::runif(length(exponents))) * 2^exponents
    values <- c(
        0, -0, 2^-1074, 3 * 2^-1074, .Machine$double.xmin,
        .Machine$double.xmin * (1 + .Machine$double.eps),
        2^-1000, 2^-500, .Machine$double.eps, 1, pi,
        2^500, 2^1000, 2^1023, .Machine$double.xmax,
        randomized[is.finite(randomized) & randomized > 0]
    )

    for (index in seq_along(values)) {
        piece <- genefunnel:::.gf_dyadic_value(values[[index]])
        restored <- genefunnel:::.gf_scale_pow2(
            sensitivity_big_as_double(piece$magnitude),
            piece$exponent
        )
        expect_identical(restored, abs(values[[index]]), info = index)
    }
    zero <- genefunnel:::.gf_dyadic_value(-0)
    expect_identical(1 / sensitivity_big_as_double(zero$magnitude), Inf)
})

test_that("exact deletion integers match an independent integer oracle", {
    set.seed(20260721)
    for (iteration in seq_len(200L)) {
        size <- sample(3:12, 1L)
        values <- floor(stats::runif(size, 0, 1001))
        expected <- sensitivity_integer_delta_numerators(values)
        exact <- genefunnel:::.gf_exact_deltas(values)

        observed <- vapply(exact$deltas, function(delta) {
            delta$sign * genefunnel:::.gf_scale_pow2(
                sensitivity_big_as_double(delta$magnitude),
                exact$exponent
            )
        }, numeric(1L))
        expect_identical(observed, expected, info = iteration)
        expect_identical(
            sensitivity_big_as_double(exact$denominator),
            as.double((size - 1L) * (size - 2L)),
            info = iteration
        )
    }
})

test_that("sorted-prefix exact deltas are identical to the brute oracle", {
    smallest <- 2^-1074
    fixtures <- list(
        c(1, 2, 7),
        c(0, 1, 1),
        c(0, 0, 0, 0),
        c(smallest, smallest, smallest),
        c(smallest, 2 * smallest, 7 * smallest),
        rep(.Machine$double.xmax, 4L),
        c(smallest, 1, 2^500, .Machine$double.xmax)
    )
    set.seed(20260725)
    for (iteration in seq_len(250L)) {
        size <- sample(3:40, 1L)
        exponents <- sample(-1074:1000, size, replace = TRUE)
        values <- (1 + stats::runif(size)) * 2^exponents
        values[!is.finite(values)] <- .Machine$double.xmax
        values[sample.int(size, floor(size / 8L))] <- 0
        fixtures[[length(fixtures) + 1L]] <- values
    }

    for (iteration in seq_along(fixtures)) {
        expect_identical(
            genefunnel:::.gf_exact_deltas_sorted(fixtures[[iteration]]),
            genefunnel:::.gf_exact_deltas(fixtures[[iteration]]),
            info = iteration
        )
    }

    for (iteration in seq_len(10L)) {
        values <- stats::rexp(128L)
        values[sample.int(128L, 16L)] <- 0
        expect_identical(
            genefunnel:::.gf_exact_deltas_sorted(values),
            genefunnel:::.gf_exact_deltas(values),
            info = paste("size 128", iteration)
        )
    }
})

test_that("exact magnitude ordering is stable and the API bypasses brute", {
    magnitudes <- lapply(c(5, 1, 5, 0, 1), function(value) {
        genefunnel:::.gf_big_from_integer(value)
    })
    expect_identical(
        genefunnel:::.gf_big_order(magnitudes),
        c(4L, 2L, 5L, 1L, 3L)
    )

    testthat::local_mocked_bindings(
        .gf_exact_deltas = function(...) {
            stop("brute oracle invoked", call. = FALSE)
        },
        .package = "genefunnel"
    )
    expect_identical(
        genefunnel:::.sensitivity_cell(c(1, 2, 7), c("A", "B", "C"))$
            largest_member,
        "B"
    )
})

test_that("exact ordering defeats subtract-rounded tie drift", {
    values <- c(
        A = 0x1.53aa55218p+6,
        B = 0x1.632ad0b2cp+4,
        C = 0x1.2bf230e48p+3,
        D = 0x1.305706508p+6
    )
    direct <- sensitivity_deltas_reference(values)
    exact <- genefunnel:::.gf_exact_deltas(unname(values))

    expect_identical(names(which.max(abs(direct))), "D")
    expect_identical(
        genefunnel:::.gf_big_compare(
            exact$deltas[[1L]]$magnitude,
            exact$deltas[[4L]]$magnitude
        ),
        0L
    )
    expect_identical(
        genefunnel:::.sensitivity_cell(unname(values), names(values))$
            largest_member,
        "A"
    )
})

test_that("wide homogeneous exponents preserve values or explicit status", {
    scales <- 2^c(-1074, -1073, -1000, -500, 0, 500, 1020)
    observed <- lapply(scales, function(scale) {
        values <- c(A = scale, B = 2 * scale, C = 7 * scale)
        genefunnel:::.sensitivity_cell(values, names(values))
    })

    expect_identical(observed[[1L]]$delta_status, "unavailable")
    expect_true(is.na(observed[[1L]]$largest_absolute_delta))
    for (index in seq_along(observed)[-1L]) {
        expect_identical(observed[[index]]$largest_member, "B", info = index)
        expect_identical(observed[[index]]$delta_status, "ordinary", info = index)
        expect_identical(
            observed[[index]]$largest_absolute_delta,
            2.5 * scales[[index]],
            info = index
        )
        expect_identical(
            observed[[index]]$median_absolute_delta,
            2.5 * scales[[index]],
            info = index
        )
        expect_identical(
            observed[[index]]$largest_delta_over_sum,
            0.25,
            info = index
        )
    }
})

test_that("exact score numerators agree with the authoritative scorer", {
    set.seed(20260722)
    mat <- matrix(stats::rexp(8L * 40L, rate = 0.2), nrow = 8L)
    mat[sample(length(mat), 20L)] <- 0
    rownames(mat) <- LETTERS[seq_len(nrow(mat))]
    observed <- genefunnel(
        mat,
        list(set = rownames(mat)),
        BPPARAM = BiocParallel::SerialParam()
    )
    expected <- vapply(seq_len(ncol(mat)), function(column) {
        sensitivity_exact_score_double(mat[, column])
    }, numeric(1L))

    expect_equal(unname(observed[1L, ]), expected, tolerance = 1e-13)
})

test_that("malformed exact arithmetic inputs fail closed", {
    one <- genefunnel:::.gf_big_from_integer(1L)

    expect_error(genefunnel:::.gf_exact_deltas(c(1, 2)), "support")
    expect_error(genefunnel:::.gf_dyadic_value(-1), "dyadic input")
    expect_error(genefunnel:::.gf_dyadic_value(Inf), "dyadic input")
    expect_error(genefunnel:::.gf_big_subtract(one, c(2)), "negative")
    expect_error(
        genefunnel:::.gf_big_ratio_double(one, 0, 0L, 1L, 3L),
        "exact ratio"
    )
    expect_error(
        genefunnel:::.gf_big_ratio_double(one, one, 0L, 2L, 3L),
        "exact ratio"
    )
    expect_error(
        genefunnel:::.gf_big_ratio_double(one, one, 0L, 1L, 0L),
        "exact ratio"
    )
})
