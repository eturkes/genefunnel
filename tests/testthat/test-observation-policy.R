test_that("sample coverage counts unique declared members and observed zeros", {
    mat <- rbind(A = c(2, 0, NA), B = c(2, 0, NaN), C = c(NA, 0, 3))
    colnames(mat) <- c("partial", "zero", "missing")
    sets <- list(partial = c("A", "A", "B", "C", "absent"), pair = c("A", "B"))
    fit <- genefunnel(mat, sets, BiocParallel::SerialParam(), details = TRUE)

    expect_named(fit, c("scores", "observed", "coverage", "set_coverage"))
    expect_identical(fit$scores, genefunnel(mat, sets, BiocParallel::SerialParam()))
    expected <- rbind(partial = c(2L, 3L, 1L), pair = c(2L, 2L, 0L))
    colnames(expected) <- colnames(mat)
    expect_identical(fit$observed, expected)
    expect_equal(fit$coverage, expected / c(4, 2))
    expect_identical(fit$set_coverage, gene_set_coverage(sets, rownames(mat)))
    expect_true(is.na(fit$scores["partial", "missing"]))
    expect_equal(fit$scores[, "zero"], c(partial = 0, pair = 0))

    sparse <- methods::as(Matrix::Matrix(mat, sparse = TRUE), "dgCMatrix")
    expect_identical(
        genefunnel(sparse, sets, BiocParallel::SerialParam(), details = TRUE),
        fit
    )
})

test_that("coverage and observed-count thresholds mask rather than remove rows", {
    mat <- rbind(A = c(2, 2), B = c(2, 2), C = c(NA, 2))
    colnames(mat) <- c("two", "three")
    sets <- list(set = c("A", "B", "C", "unmatched"))
    fit <- genefunnel(
        mat, sets, BiocParallel::SerialParam(), min_coverage = 0.75,
        details = TRUE
    )
    expect_identical(dim(fit$scores), c(1L, 2L))
    expect_true(is.na(fit$scores[1, 1]))
    expect_equal(unname(fit$scores[1, 2]), 6)
    expect_equal(unname(fit$coverage), matrix(c(0.5, 0.75), 1))
    expect_identical(
        genefunnel(mat, sets, BiocParallel::SerialParam(), min_observed = 3),
        fit$scores
    )
    expect_true(all(is.na(genefunnel(
        mat, sets, BiocParallel::SerialParam(), min_coverage = 1
    ))))
    expect_true(all(is.na(genefunnel(
        mat, sets, BiocParallel::SerialParam(), min_observed = 4
    ))))
    for (input in list(mat, Matrix::Matrix(mat, sparse = TRUE))) {
        expect_identical(
            genefunnel(input, sets, BiocParallel::SerialParam(),
                min_coverage = 0.5),
            genefunnel(input, sets, BiocParallel::SerialParam())
        )
    }
})

test_that("coverage rejection precedes potentially overflowing arithmetic", {
    mat <- matrix(.Machine$double.xmax, 2, 1, dimnames = list(c("A", "B"), "s"))
    sets <- list(set = c("A", "B", "absent"))
    for (input in list(mat, Matrix::Matrix(mat, sparse = TRUE))) {
        fit <- genefunnel(input, sets, BiocParallel::SerialParam(),
            min_coverage = 1, details = TRUE)
        expect_true(is.na(fit$scores[1, 1]))
        expect_equal(fit$observed[1, 1], 2L)
    }
})

test_that("fractional coverage accepts its exact reported boundary", {
    members <- paste0("g", 1:50)
    mat <- matrix(1, 7, 1, dimnames = list(members[1:7], "s"))
    fit <- genefunnel(mat, list(set = members), BiocParallel::SerialParam(),
        min_coverage = 7 / 50, details = TRUE)
    expect_equal(unname(fit$scores[1, 1]), 7)
    expect_identical(unname(fit$coverage[1, 1]), 7 / 50)
})

test_that("coverage policies reject malformed scalar options", {
    mat <- matrix(1, 2, 1, dimnames = list(c("A", "B"), NULL))
    sets <- list(set = c("A", "B"))
    for (value in list(NULL, NA, NaN, Inf, -0.1, 1.1, TRUE,
        c(0, 1), matrix(0.5), structure(0.5, class = "custom"))) {
        expect_error(genefunnel(mat, sets, BiocParallel::SerialParam(),
            min_coverage = value), "min_coverage")
    }
    for (value in list(NULL, NA, Inf, 0, 1, 2.5, TRUE, c(2, 3),
        .Machine$integer.max + 1, matrix(2))) {
        expect_error(genefunnel(mat, sets, BiocParallel::SerialParam(),
            min_observed = value), "min_observed")
    }
    for (value in list(NULL, NA, 1, c(TRUE, FALSE), matrix(TRUE))) {
        expect_error(genefunnel(mat, sets, BiocParallel::SerialParam(),
            details = value), "details")
    }
})

test_that("empty retained sets keep diagnostic dimensions and unmatched coverage", {
    mat <- matrix(1, 2, 2, dimnames = list(c("A", "B"), c("s1", "s2")))
    sets <- list(none = "Z", singleton = "A", empty = character())
    expect_warning(
        fit <- genefunnel(mat, sets, BiocParallel::SerialParam(), details = TRUE),
        "skipped 3 gene sets"
    )
    for (field in c("scores", "observed", "coverage")) {
        expect_identical(dim(fit[[field]]), c(0L, 2L))
        expect_identical(colnames(fit[[field]]), colnames(mat))
    }
    expect_type(fit$observed, "integer")
    expect_identical(fit$set_coverage, gene_set_coverage(sets, rownames(mat)))
})

test_that("diagnostics survive reusable socket workers and sample permutation", {
    mat <- matrix(rep(c(0, 2, NA, 3, 1), 12), 5,
        dimnames = list(LETTERS[1:5], paste0("s", 1:12)))
    sets <- list(all = LETTERS[1:5], partial = c("A", "B", "C", "absent"))
    expected <- genefunnel(mat, sets, BiocParallel::SerialParam(),
        min_coverage = 0.5, details = TRUE)
    backend <- BiocParallel::SnowParam(2, type = "SOCK", progressbar = FALSE)
    BiocParallel::bpstart(backend)
    on.exit(BiocParallel::bpstop(backend))
    for (input in list(mat, Matrix::Matrix(mat, sparse = TRUE))) {
        expect_identical(genefunnel(input, sets, backend,
            min_coverage = 0.5, details = TRUE), expected)
    }
    permuted <- genefunnel(mat[, 12:1], sets, BiocParallel::SerialParam(),
        min_coverage = 0.5, details = TRUE)
    for (field in c("scores", "observed", "coverage")) {
        expect_identical(permuted[[field]], expected[[field]][, 12:1])
    }
})
