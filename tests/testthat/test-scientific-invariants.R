# Assisted-by: OpenAI Codex.

scientific_fixture <- function() {
    mat <- cbind(
        ordinary = c(1, 2, 4, 8, 3, 5),
        zero_rich = c(0, 0, 2, 0, 6, 0),
        missing = c(NA_real_, 3, 1, NaN, 4, 0),
        equal = c(5, 5, 5, 5, 5, 5)
    )
    rownames(mat) <- LETTERS[1:6]
    mat
}

scientific_gene_sets <- function() {
    list(
        alpha = c("A", "B", "C"),
        beta = c("C", "D", "E"),
        partial = c("F", "A", "absent")
    )
}

serial_scores <- function(mat, gene_sets) {
    genefunnel(
        mat,
        gene_sets,
        BPPARAM = BiocParallel::SerialParam()
    )
}

test_that("member, row, and sample permutations preserve corresponding scores", {
    mat <- scientific_fixture()
    gene_sets <- scientific_gene_sets()
    expected <- serial_scores(mat, gene_sets)

    member_permuted <- lapply(gene_sets, rev)
    expect_equal(
        serial_scores(mat, member_permuted),
        expected,
        tolerance = 5e-15
    )

    row_order <- c(6L, 2L, 4L, 1L, 5L, 3L)
    expect_identical(
        serial_scores(mat[row_order, , drop = FALSE], gene_sets),
        expected
    )

    sample_order <- c(4L, 2L, 1L, 3L)
    expect_identical(
        serial_scores(mat[, sample_order, drop = FALSE], gene_sets),
        expected[, sample_order, drop = FALSE]
    )

    irrelevant <- rbind(
        irrelevant = c(100, 0, NA_real_, 2),
        mat
    )
    expect_identical(serial_scores(irrelevant, gene_sets), expected)
})

test_that("other samples cannot affect a retained sample", {
    mat <- scientific_fixture()
    gene_sets <- scientific_gene_sets()
    expected <- serial_scores(mat, gene_sets)

    extended <- cbind(
        before = c(100, 0, 1, 2, 3, 4),
        mat,
        after = c(0, 8, 0, 8, 0, 8)
    )
    expect_identical(
        serial_scores(extended, gene_sets)[, colnames(mat), drop = FALSE],
        expected
    )

    retained <- c("ordinary", "missing")
    expect_identical(
        serial_scores(mat[, retained, drop = FALSE], gene_sets),
        expected[, retained, drop = FALSE]
    )

    modified <- mat
    modified[, "zero_rich"] <- c(1000, 2, 3, 4, 5, 6)
    expect_identical(
        serial_scores(modified, gene_sets)[, "ordinary", drop = FALSE],
        expected[, "ordinary", drop = FALSE]
    )
})

test_that("overlapping and other gene sets cannot affect a retained set", {
    mat <- scientific_fixture()
    gene_sets <- scientific_gene_sets()
    expected <- serial_scores(mat, gene_sets)

    for (set_name in names(gene_sets)) {
        isolated <- serial_scores(mat, gene_sets[set_name])
        expect_identical(isolated, expected[set_name, , drop = FALSE])
    }

    with_extra <- c(
        list(extra = c("B", "D", "F")),
        gene_sets
    )
    expect_identical(
        serial_scores(mat, with_extra)[names(gene_sets), , drop = FALSE],
        expected
    )

    expect_identical(
        serial_scores(mat, gene_sets[c("partial", "alpha")]),
        expected[c("partial", "alpha"), , drop = FALSE]
    )
    expect_identical(
        serial_scores(mat, gene_sets[rev(names(gene_sets))]),
        expected[rev(names(gene_sets)), , drop = FALSE]
    )

    modified <- gene_sets
    modified$beta <- c("A", "E", "F")
    unaffected <- c("alpha", "partial")
    expect_identical(
        serial_scores(mat, modified)[unaffected, , drop = FALSE],
        expected[unaffected, , drop = FALSE]
    )
})

test_that("scores are positively homogeneous and translation-covariant", {
    mat <- scientific_fixture()
    gene_sets <- scientific_gene_sets()
    expected <- serial_scores(mat, gene_sets)

    for (constant in c(0, 1e-100, 0.25, 7, 1e100)) {
        expect_equal(
            serial_scores(mat * constant, gene_sets),
            expected * constant,
            tolerance = 5e-13
        )
    }

    shift <- 2.75
    shifted <- mat
    shifted[!is.na(shifted)] <- shifted[!is.na(shifted)] + shift
    observed_sizes <- reference_cells(mat, gene_sets, function(values) {
        n_observed <- sum(!is.na(values))
        if (n_observed < 2L) NA_real_ else as.double(n_observed)
    })
    expect_equal(
        serial_scores(shifted, gene_sets),
        expected + shift * observed_sizes,
        tolerance = 5e-15
    )
})
