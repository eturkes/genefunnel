test_that("complete finite dense input preserves shape, names, and scores", {
    mat <- matrix(
        c(
            4, 4, 4,
            0, 0, 0,
            4, 0, 0
        ),
        nrow = 3L,
        dimnames = list(
            c("A", "B", "C"),
            c("equal", "zero", "one_nonzero")
        )
    )
    gene_sets <- list(triplet = c("A", "B", "C"))

    observed <- genefunnel(
        mat,
        gene_sets,
        BPPARAM = BiocParallel::SerialParam()
    )
    expected <- vapply(
        seq_len(ncol(mat)),
        function(column) reference_score(mat[, column]),
        numeric(1)
    )

    expect_true(is.matrix(observed))
    expect_type(observed, "double")
    expect_identical(dim(observed), c(1L, 3L))
    expect_identical(
        dimnames(observed),
        list("triplet", c("equal", "zero", "one_nonzero"))
    )
    expect_equal(unname(observed), matrix(expected, nrow = 1L), tolerance = 0)
    expect_identical(unname(observed[1L, ]), c(12, 0, 0))
})
