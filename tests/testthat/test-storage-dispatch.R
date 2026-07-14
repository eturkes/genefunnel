storage_fixture <- function() {
    matrix(
        c(
            1, 2, 3, 4, 5,
            0, 0, 0, 0, 0,
            1, NA_real_, 2, NaN, 0,
            8, 4, 2, 1, 0
        ),
        nrow = 5L,
        dimnames = list(
            c("A", "B", "C", "D", "E"),
            c("ordinary", "all_zero", "missing", "unequal")
        )
    )
}

storage_gene_sets <- function() {
    list(
        partial = c("A", "C", "absent"),
        overlap = c("E", "C", "A"),
        pair = c("D", "B")
    )
}

test_that("dense and sparse inputs reach distinct native kernels", {
    dense <- storage_fixture()
    sparse <- Matrix::Matrix(dense, sparse = TRUE)
    calls <- new.env(parent = emptyenv())
    calls$dense <- list()
    calls$sparse <- list()

    testthat::local_mocked_bindings(
        calculateScoresDense = function(orig_mat, gene_indices) {
            calls$dense[[length(calls$dense) + 1L]] <- orig_mat
            matrix(0, nrow = length(gene_indices), ncol = ncol(orig_mat))
        },
        calculateScoresSparse = function(orig_mat, gene_indices) {
            calls$sparse[[length(calls$sparse) + 1L]] <- orig_mat
            matrix(0, nrow = length(gene_indices), ncol = ncol(orig_mat))
        },
        .package = "genefunnel"
    )

    genefunnel(
        dense,
        storage_gene_sets(),
        BPPARAM = BiocParallel::SerialParam()
    )
    expect_length(calls$dense, ncol(dense))
    expect_length(calls$sparse, 0L)
    expect_true(all(vapply(calls$dense, is.matrix, logical(1))))
    expect_true(all(vapply(calls$dense, ncol, integer(1)) == 1L))

    calls$dense <- list()
    genefunnel(
        sparse,
        storage_gene_sets(),
        BPPARAM = BiocParallel::SerialParam()
    )
    expect_length(calls$dense, 0L)
    expect_length(calls$sparse, ncol(sparse))
    expect_true(all(vapply(
        calls$sparse,
        inherits,
        logical(1),
        what = "sparseMatrix"
    )))
    expect_true(all(vapply(calls$sparse, ncol, integer(1)) == 1L))
})

test_that("dense Matrix and sparse scores agree across storage semantics", {
    dense <- storage_fixture()
    matrix_dense <- Matrix::Matrix(dense, sparse = FALSE)
    sparse <- Matrix::Matrix(dense, sparse = TRUE)
    gene_sets <- storage_gene_sets()

    expected <- genefunnel(
        dense,
        gene_sets,
        BPPARAM = BiocParallel::SerialParam()
    )
    observed_dense <- genefunnel(
        matrix_dense,
        gene_sets,
        BPPARAM = BiocParallel::SerialParam()
    )
    observed_sparse <- genefunnel(
        sparse,
        gene_sets,
        BPPARAM = BiocParallel::SerialParam()
    )

    expect_identical(observed_dense, expected)
    expect_identical(observed_sparse, expected)
    expect_identical(unname(expected[, "all_zero"]), c(0, 0, 0))
    expect_identical(
        unname(expected[, "missing"]),
        c(
            reference_score(c(1, 2)),
            reference_score(c(0, 2, 1)),
            NA_real_
        )
    )
})

test_that("dense and sparse paths preserve reordered rows and columns", {
    dense <- storage_fixture()
    gene_sets <- storage_gene_sets()
    row_order <- c(5L, 2L, 4L, 1L, 3L)
    column_order <- c(4L, 2L, 1L, 3L)
    reordered <- dense[row_order, column_order, drop = FALSE]
    sparse_reordered <- Matrix::Matrix(reordered, sparse = TRUE)

    expected <- genefunnel(
        dense,
        gene_sets,
        BPPARAM = BiocParallel::SerialParam()
    )
    observed_dense <- genefunnel(
        reordered,
        gene_sets,
        BPPARAM = BiocParallel::SerialParam()
    )
    observed_sparse <- genefunnel(
        sparse_reordered,
        gene_sets,
        BPPARAM = BiocParallel::SerialParam()
    )

    expect_identical(observed_dense[, colnames(dense), drop = FALSE], expected)
    expect_identical(observed_sparse, observed_dense)
})
