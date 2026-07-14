test_that("global and sample-specific minimum sizes remain distinct", {
    mat <- cbind(
        no_observed = c(NA_real_, NaN, NA_real_),
        one_observed = c(4, NA_real_, NaN),
        two_observed = c(1, 2, NA_real_),
        three_observed = c(1, 2, 3)
    )
    rownames(mat) <- c("A", "B", "C")
    gene_sets <- list(
        no_match = c("X", "Y"),
        one_match = c("A", "X"),
        pair = c("A", "B"),
        triplet = c("A", "B", "C")
    )

    expect_warning(
        observed <- genefunnel(
            mat,
            gene_sets,
            BPPARAM = BiocParallel::SerialParam()
        ),
        "Omitting 2 gene sets"
    )

    expect_identical(dim(observed), c(2L, 4L))
    expect_identical(rownames(observed), c("pair", "triplet"))
    expect_identical(colnames(observed), colnames(mat))
    expect_true(all(is.na(observed[, c("no_observed", "one_observed")])))
    expect_identical(unname(observed[, "two_observed"]), c(2, 2))
    expect_identical(unname(observed[, "three_observed"]), c(2, 4.5))
})

test_that("extreme finite values and signed zero preserve score semantics", {
    smallest <- 2^-1074
    large <- .Machine$double.xmax / 4
    mat <- cbind(
        signed_zero = c(-0, 0, -0, 0),
        tiny = c(smallest, smallest, 0, 0),
        huge = c(large, large, 0, 0),
        dynamic = c(large, 1, 0, 0)
    )
    rownames(mat) <- c("A", "B", "C", "D")
    gene_sets <- list(pair = c("A", "B"), quartet = rownames(mat))
    expected <- reference_scores(mat, gene_sets)

    dense <- genefunnel(
        mat,
        gene_sets,
        BPPARAM = BiocParallel::SerialParam()
    )
    sparse <- genefunnel(
        Matrix::Matrix(mat, sparse = TRUE),
        gene_sets,
        BPPARAM = BiocParallel::SerialParam()
    )

    expect_equal(dense, expected, tolerance = 0)
    expect_equal(sparse, expected, tolerance = 0)
    expect_true(all(is.finite(dense)))
    expect_identical(unname(dense[, "signed_zero"]), c(0, 0))
    reciprocals <- 1 / dense[, "signed_zero"]
    expect_true(all(is.infinite(reciprocals) & reciprocals > 0))
})

test_that("integer and zero-only sparse matrices retain numeric semantics", {
    integer_mat <- matrix(
        as.integer(c(
            .Machine$integer.max, .Machine$integer.max, 0,
            4, 0, 0,
            1, 2, NA
        )),
        nrow = 3L,
        dimnames = list(c("A", "B", "C"), c("large", "one", "missing"))
    )
    gene_sets <- list(pair = c("A", "B"), triplet = c("A", "B", "C"))
    observed_integer <- genefunnel(
        integer_mat,
        gene_sets,
        BPPARAM = BiocParallel::SerialParam()
    )
    expect_identical(typeof(integer_mat), "integer")
    expect_equal(observed_integer, reference_scores(integer_mat, gene_sets))

    sparse_zero <- Matrix::sparseMatrix(
        i = integer(),
        j = integer(),
        x = double(),
        dims = c(3L, 4L),
        dimnames = list(c("A", "B", "C"), paste0("sample_", 1:4))
    )
    observed_zero <- genefunnel(
        sparse_zero,
        gene_sets,
        BPPARAM = BiocParallel::SerialParam()
    )
    expect_length(methods::slot(sparse_zero, "x"), 0L)
    expect_identical(
        observed_zero,
        matrix(0, nrow = 2L, ncol = 4L, dimnames = dimnames(observed_zero))
    )
})

test_that("unnamed samples and omitted sets preserve stable matrix dimensions", {
    mat <- matrix(
        c(1, 2, 3),
        ncol = 1L,
        dimnames = list(c("A", "B", "C"), NULL)
    )
    gene_sets <- list(
        omitted_before = c("A", "absent"),
        retained = c("A", "B", "C"),
        omitted_after = c("X", "Y")
    )

    expect_warning(
        dense <- genefunnel(
            mat,
            gene_sets,
            BPPARAM = BiocParallel::SerialParam()
        ),
        "Omitting 2 gene sets"
    )
    expect_warning(
        sparse <- genefunnel(
            Matrix::Matrix(mat, sparse = TRUE),
            gene_sets,
            BPPARAM = BiocParallel::SerialParam()
        ),
        "Omitting 2 gene sets"
    )

    expect_identical(dim(dense), c(1L, 1L))
    expect_identical(dimnames(dense), list("retained", NULL))
    expect_identical(sparse, dense)

    expect_warning(
        empty <- genefunnel(
            mat,
            gene_sets[c("omitted_before", "omitted_after")],
            BPPARAM = BiocParallel::SerialParam()
        ),
        "Omitting 2 gene sets"
    )
    expect_identical(dim(empty), c(0L, 1L))
    expect_null(rownames(empty))
    expect_null(colnames(empty))
})

test_that("scoring does not mutate matrices or gene-set declarations", {
    dense <- matrix(
        c(1, 2, 0, NA_real_, 3, 4),
        nrow = 3L,
        dimnames = list(c("A", "B", "C"), c("sample_1", "sample_2"))
    )
    sparse <- Matrix::sparseMatrix(
        i = c(1L, 1L, 2L, 3L),
        j = c(1L, 1L, 1L, 2L),
        x = c(1, 2, 2, NA_real_),
        dims = c(3L, 2L),
        dimnames = dimnames(dense),
        repr = "T"
    )
    gene_sets <- list(
        duplicated = c("A", "B", "A"),
        partial = c("C", "B", "absent")
    )
    dense_snapshot <- serialize(dense, NULL)
    sparse_snapshot <- serialize(sparse, NULL)
    sets_snapshot <- serialize(gene_sets, NULL)

    genefunnel(dense, gene_sets, BPPARAM = BiocParallel::SerialParam())
    genefunnel(sparse, gene_sets, BPPARAM = BiocParallel::SerialParam())

    expect_identical(serialize(dense, NULL), dense_snapshot)
    expect_identical(serialize(sparse, NULL), sparse_snapshot)
    expect_identical(serialize(gene_sets, NULL), sets_snapshot)
})
