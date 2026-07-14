coverage_fixture <- function() {
    list(
        full = c("A", "B"),
        half = c("A", "B", "C", "D"),
        one = c("A", "C"),
        duplicated = c("B", "A", "B", "missing", "A"),
        empty = character()
    )
}

test_that("coverage reports canonical set facts in input order", {
    observed <- gene_set_coverage(coverage_fixture(), c("A", "B"))

    expect_s3_class(observed, "data.frame")
    expect_identical(
        names(observed),
        c(
            "gene_set",
            "declared_size",
            "matched_size",
            "unmatched_size",
            "coverage",
            "duplicate_member_count",
            "scoreable"
        )
    )
    expect_identical(observed$gene_set, names(coverage_fixture()))
    expect_identical(observed$declared_size, c(2L, 4L, 2L, 3L, 0L))
    expect_identical(observed$matched_size, c(2L, 2L, 1L, 2L, 0L))
    expect_identical(observed$unmatched_size, c(0L, 2L, 1L, 1L, 0L))
    expect_equal(observed$coverage, c(1, 0.5, 0.5, 2 / 3, NA_real_))
    expect_identical(observed$duplicate_member_count, c(0L, 0L, 0L, 2L, 0L))
    expect_identical(observed$scoreable, c(TRUE, TRUE, FALSE, TRUE, FALSE))
})

test_that("coverage preserves multiple entirely empty sets", {
    observed <- gene_set_coverage(
        list(first = character(), second = character()),
        c("A", "B")
    )

    expect_identical(observed$gene_set, c("first", "second"))
    expect_identical(observed$declared_size, c(0L, 0L))
    expect_identical(observed$matched_size, c(0L, 0L))
    expect_true(all(is.na(observed$coverage)))
    expect_false(any(observed$scoreable))
})

test_that("coverage matching is exact and case-sensitive", {
    observed <- gene_set_coverage(
        list(case = c("A", "a"), unmatched = "missing"),
        c("A", "B")
    )

    expect_identical(observed$matched_size, c(1L, 0L))
    expect_equal(observed$coverage, c(0.5, 0))
    expect_false(any(observed$scoreable))
})

test_that("callers can apply strict and 50-percent coverage policies", {
    gene_sets <- coverage_fixture()
    coverage <- gene_set_coverage(gene_sets, c("A", "B"))

    strict <- gene_sets[which(coverage$coverage == 1)]
    half <- gene_sets[which(
        coverage$coverage >= 0.5 & coverage$matched_size >= 2L
    )]

    expect_identical(names(strict), "full")
    expect_identical(names(half), c("full", "half", "duplicated"))

    mat <- matrix(
        c(1, 2),
        nrow = 2L,
        dimnames = list(c("A", "B"), "sample")
    )
    strict_scores <- genefunnel(
        mat,
        strict,
        BPPARAM = BiocParallel::SerialParam()
    )
    half_scores <- genefunnel(
        mat,
        half,
        BPPARAM = BiocParallel::SerialParam()
    )
    expect_identical(rownames(strict_scores), "full")
    expect_identical(rownames(half_scores), c("full", "half", "duplicated"))
})

test_that("partial sets score matched unique members and retain set order", {
    mat <- matrix(
        c(1, 2, 4, 4),
        nrow = 2L,
        dimnames = list(c("A", "B"), c("unequal", "equal"))
    )
    gene_sets <- list(
        no_match = c("X", "Y"),
        partial = c("A", "B", "C"),
        one_match = c("A", "X"),
        duplicated = c("B", "A", "B", "missing", "A")
    )

    warnings <- character()
    observed <- withCallingHandlers(
        genefunnel(
            mat,
            gene_sets,
            BPPARAM = BiocParallel::SerialParam()
        ),
        warning = function(condition) {
            warnings <<- c(warnings, conditionMessage(condition))
            invokeRestart("muffleWarning")
        }
    )

    expect_length(warnings, 1L)
    expect_match(warnings, "Omitting 2 gene sets")
    expect_identical(rownames(observed), c("partial", "duplicated"))
    expect_identical(colnames(observed), c("unequal", "equal"))
    expect_identical(unname(observed), matrix(c(2, 2, 8, 8), nrow = 2L))

    prepared <- genefunnel:::.prepare_gene_sets(gene_sets, rownames(mat))
    expect_identical(prepared$indices$duplicated, c(2L, 1L))
})

test_that("all globally unscoreable sets yield an empty score matrix", {
    mat <- matrix(
        c(1, 2, 3, 4),
        nrow = 2L,
        dimnames = list(c("A", "B"), c("sample_1", "sample_2"))
    )

    expect_warning(
        observed <- genefunnel(
            mat,
            list(empty = character(), one = c("A", "missing")),
            BPPARAM = BiocParallel::SerialParam()
        ),
        "Omitting 2 gene sets"
    )
    expect_true(is.matrix(observed))
    expect_type(observed, "double")
    expect_identical(dim(observed), c(0L, 2L))
    expect_null(rownames(observed))
    expect_identical(colnames(observed), c("sample_1", "sample_2"))
})
