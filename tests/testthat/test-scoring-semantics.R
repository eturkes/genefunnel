test_that("native scores match the oracle with sample-specific missingness", {
    values <- list(
        equal = c(4, 4, 4),
        zeros = c(0, 0, 0),
        unequal_pair = c(4, 0, NA_real_),
        one_nonzero = c(4, 0, 0),
        one_missing = c(1, 2, NA_real_),
        one_nan = c(1, 2, NaN),
        insufficient = c(4, NA_real_, NaN)
    )
    mat <- matrix(
        unlist(values, use.names = FALSE),
        nrow = 3L,
        dimnames = list(c("A", "B", "C"), names(values))
    )

    observed <- genefunnel(
        mat,
        list(triplet = c("A", "B", "C")),
        BPPARAM = BiocParallel::SerialParam()
    )
    expected <- vapply(values, reference_score, numeric(1))

    expect_equal(unname(observed[1L, ]), unname(expected))
    expect_identical(
        unname(observed[1L, c("equal", "zeros", "unequal_pair", "one_nonzero")]),
        c(12, 0, 0, 0)
    )
    expect_identical(
        unname(observed[1L, c("one_missing", "one_nan")]),
        c(2, 2)
    )
    expect_true(is.na(observed[1L, "insufficient"]))
})

test_that("zeros remain observed in dense and sparse inputs", {
    dense <- matrix(
        c(4, 0),
        nrow = 2L,
        dimnames = list(c("A", "B"), "sample")
    )
    sparse <- Matrix::sparseMatrix(
        i = 1L,
        j = 1L,
        x = 4,
        dims = c(2L, 1L),
        dimnames = dimnames(dense)
    )
    gene_sets <- list(pair = c("A", "B"))

    dense_score <- genefunnel(
        dense,
        gene_sets,
        BPPARAM = BiocParallel::SerialParam()
    )
    sparse_score <- genefunnel(
        sparse,
        gene_sets,
        BPPARAM = BiocParallel::SerialParam()
    )
    expected <- reference_score(c(4, 0))

    expect_identical(unname(dense_score), matrix(expected, nrow = 1L))
    expect_identical(sparse_score, dense_score)
})

test_that("stored sparse missing values are omitted without dropping zeros", {
    sparse <- Matrix::sparseMatrix(
        i = c(1L, 2L),
        j = c(1L, 1L),
        x = c(1, NA_real_),
        dims = c(3L, 1L),
        dimnames = list(c("A", "B", "C"), "sample")
    )

    observed <- genefunnel(
        sparse,
        list(triplet = c("A", "B", "C")),
        BPPARAM = BiocParallel::SerialParam()
    )
    expected <- reference_score(c(1, NA_real_, 0))

    expect_identical(unname(observed), matrix(expected, nrow = 1L))
    expect_identical(expected, 0)
})

test_that("scale-aware cleanup preserves legitimate small positive scores", {
    values <- list(
        equal_tiny = c(1e-12, 1e-12),
        imbalanced = c(1, 1e-12),
        near_precision = c(1, 1e-16)
    )
    mat <- matrix(
        unlist(values, use.names = FALSE),
        nrow = 2L,
        dimnames = list(c("A", "B"), names(values))
    )

    observed <- genefunnel(
        mat,
        list(pair = c("A", "B")),
        BPPARAM = BiocParallel::SerialParam()
    )
    expected <- vapply(values, reference_score, numeric(1))

    expect_equal(unname(observed[1L, ]), unname(expected))
    expect_true(all(observed > 0))
})

test_that("stable arithmetic preserves scores across extreme dynamic range", {
    values <- c(1e150, 1, 0)
    mat <- matrix(
        values,
        nrow = 3L,
        dimnames = list(c("A", "B", "C"), "sample")
    )

    observed <- genefunnel(
        mat,
        list(triplet = c("A", "B", "C")),
        BPPARAM = BiocParallel::SerialParam()
    )
    expected <- reference_score(values)

    expect_identical(unname(observed), matrix(expected, nrow = 1L))
    expect_identical(expected, 1.5)
})

test_that("the native boundary rejects escaped infinities with context", {
    mat <- Matrix::sparseMatrix(
        i = c(1L, 2L),
        j = c(1L, 1L),
        x = c(Inf, 1),
        dims = c(2L, 1L)
    )

    expect_error(
        genefunnel:::calculateScores(mat, list(c(1L, 2L))),
        "infinite.*row 1.*column 1.*gene set 1"
    )
})

test_that("the native boundary rejects invalid arithmetic results", {
    negative <- Matrix::sparseMatrix(
        i = c(1L, 2L),
        j = c(1L, 1L),
        x = c(-1, -1),
        dims = c(2L, 1L)
    )
    overflowing <- Matrix::sparseMatrix(
        i = c(1L, 2L),
        j = c(1L, 1L),
        x = rep(.Machine$double.xmax, 2L),
        dims = c(2L, 1L)
    )

    expect_error(
        genefunnel:::calculateScores(negative, list(c(1L, 2L))),
        "materially negative.*gene set 1.*column 1"
    )
    expect_error(
        genefunnel:::calculateScores(overflowing, list(c(1L, 2L))),
        "non-finite.*gene set 1.*column 1"
    )
})
