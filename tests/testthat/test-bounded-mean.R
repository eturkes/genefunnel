test_that("large equal finite cells do not fail after mean accumulation", {
    value <- 2 - .Machine$double.eps
    for (n in c(2049L, 4097L, 10000L, 65536L)) {
        mat <- matrix(value, n, 1L,
            dimnames = list(paste0("g", seq_len(n)), "equal"))
        sets <- list(set = rownames(mat))
        for (input in list(mat, Matrix::Matrix(mat, sparse = TRUE))) {
            fit <- genefunnel(input, sets, BiocParallel::SerialParam(),
                details = TRUE)
            expect_true(is.finite(fit$scores[[1L]]))
            expect_equal(fit$scores[[1L]], n * value, tolerance = 5e-15)
            expect_identical(fit$observed[[1L]], n)
            expect_identical(fit$coverage[[1L]], 1)
        }
    }
})

test_that("sparse equal cells remain valid with missing values and stored zeros", {
    n <- 2049L
    value <- 2 - .Machine$double.eps
    mat <- cbind(
        equal = c(rep(value, n), NA_real_, NaN),
        zero = c(rep(value, n), 0, 0),
        all_zero = rep(0, n + 2L)
    )
    rownames(mat) <- paste0("g", seq_len(nrow(mat)))
    sets <- list(reverse = rev(rownames(mat)),
        duplicate = c(rownames(mat), rownames(mat)[1L]))
    expected_counts <- matrix(rep(c(n, n + 2L, n + 2L), each = 2L),
        2L, dimnames = list(names(sets), colnames(mat)))
    for (input in list(Matrix::Matrix(mat, sparse = TRUE), mat)) {
        fit <- genefunnel(input, sets, BiocParallel::SerialParam(),
            details = TRUE)
        expect_identical(fit$observed, expected_counts)
        expect_equal(unname(fit$scores[, "equal"]), rep(n * value, 2L),
            tolerance = 5e-15)
        expect_identical(unname(fit$scores[, "all_zero"]), c(0, 0))
        expect_true(all(is.finite(fit$scores)))
        expect_identical(dimnames(fit$scores), dimnames(expected_counts))
    }
})

test_that("a bounded mean preserves masking and genuine result overflow", {
    n <- 2049L
    mat <- matrix(.Machine$double.xmax, n, 1L,
        dimnames = list(paste0("g", seq_len(n)), "overflow"))
    sets <- list(set = rownames(mat))
    for (input in list(mat, Matrix::Matrix(mat, sparse = TRUE))) {
        fit <- genefunnel(input, sets, BiocParallel::SerialParam(),
            min_observed = n + 1L, details = TRUE)
        expect_true(is.na(fit$scores[[1L]]))
        expect_identical(fit$observed[[1L]], n)
        expect_error(genefunnel(input, sets, BiocParallel::SerialParam()),
            "non-finite result")
    }
})
