experiment_score_fixture <- function() {
    rbind(
        A = c(s2 = 0, s1 = 3, s3 = NA),
        B = c(s2 = 2, s1 = 3, s3 = 1),
        C = c(s2 = 2, s1 = NaN, s3 = 4)
    )
}

test_that("experiment inputs preserve matrix scores and diagnostics", {
    skip_if_not_installed("SummarizedExperiment")
    mat <- experiment_score_fixture()
    sets <- list(all = c("C", "A", "B", "unmatched"), pair = c("A", "B"))
    for (data in list(mat, Matrix::Matrix(mat, sparse = TRUE))) {
        se <- SummarizedExperiment::SummarizedExperiment(
            assays = list(signal = data),
            colData = data.frame(batch = 1:3, row.names = colnames(mat)),
            metadata = list(source = "unchanged")
        )
        before <- serialize(se, NULL)
        for (details in c(FALSE, TRUE)) {
            expected <- genefunnel(data, sets,
                BPPARAM = BiocParallel::SerialParam(),
                min_coverage = 0.7, details = details)
            expect_identical(genefunnel(se, sets,
                BPPARAM = BiocParallel::SerialParam(),
                min_coverage = 0.7, details = details), expected)
        }
        extracted <- genefunnel:::.score_input_matrix(se, NULL)
        expect_identical(extracted, data)
        expect_identical(serialize(se, NULL), before)
    }
})

test_that("assay selection is explicit and only selected values are validated", {
    skip_if_not_installed("SummarizedExperiment")
    mat <- experiment_score_fixture()
    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(invalid = -abs(mat), signal = mat)
    )
    score <- function(x, ...) genefunnel(x, list(all = rownames(mat)),
        BPPARAM = BiocParallel::SerialParam(), ...)
    expect_error(score(se), "multiple assays")
    expect_identical(score(se, assay = "signal"), score(mat))
    expect_error(score(se, assay = "invalid"), "negative values")
    expect_error(score(se, assay = "sign"), "exactly one assay name")
    expect_error(score(se, assay = "missing"), "exactly one assay name")
    for (selector in list(NA_character_, "", 1L, TRUE,
        c("signal", "invalid"), factor("signal"), matrix("signal"))) {
        expect_error(score(se, assay = selector), "non-empty assay name")
    }
    SummarizedExperiment::assayNames(se) <- c("signal", "signal")
    expect_error(score(se, assay = "signal"), "exactly one assay name")
    expect_error(score(mat, assay = "signal"), "NULL for matrix input")
    expect_error(score(Matrix::Matrix(mat), assay = "signal"),
        "NULL for matrix input")
})

test_that("container dimnames define exact feature matching and sample order", {
    skip_if_not_installed("SummarizedExperiment")
    mat <- experiment_score_fixture()
    se <- SummarizedExperiment::SummarizedExperiment(assays = list(mat))
    rownames(se) <- c("gene_c", "gene_a", "gene_b")
    colnames(se) <- c("sample_z", "sample_x", "sample_y")
    sets <- list(pathway = c("gene_a", "gene_c"))
    expected <- genefunnel(SummarizedExperiment::assay(se), sets,
        BPPARAM = BiocParallel::SerialParam(), details = TRUE)
    expect_identical(genefunnel(se, sets,
        BPPARAM = BiocParallel::SerialParam(), details = TRUE), expected)
    expect_identical(colnames(expected$scores), colnames(se))
    expect_identical(expected$set_coverage$matched_size, 2L)
})

test_that("empty and unsupported assays fail without implicit coercion", {
    skip_if_not_installed("SummarizedExperiment")
    score <- function(x) genefunnel(x, list(all = c("A", "B")),
        BPPARAM = BiocParallel::SerialParam())
    expect_error(score(SummarizedExperiment::SummarizedExperiment()),
        "at least one assay")
    bad <- experiment_score_fixture() > 0
    expect_error(score(SummarizedExperiment::SummarizedExperiment(
        assays = list(bad))), "numeric")
    bad <- experiment_score_fixture()
    rownames(bad) <- c("A", "A", "B")
    expect_error(score(SummarizedExperiment::SummarizedExperiment(
        assays = list(bad))), "duplicated identifiers")
})

test_that("subclasses and unnamed single assays use the same adapter", {
    skip_if_not_installed("SummarizedExperiment")
    mat <- experiment_score_fixture()
    se <- SummarizedExperiment::SummarizedExperiment(assays = list(mat))
    ranged <- methods::as(se, "RangedSummarizedExperiment")
    sets <- list(pair = c("A", "B"))
    for (input in list(se, ranged, se[, 1L, drop = FALSE])) {
        expected <- genefunnel(SummarizedExperiment::assay(input), sets,
            BPPARAM = BiocParallel::SerialParam())
        expect_identical(genefunnel(input, sets,
            BPPARAM = BiocParallel::SerialParam()), expected)
    }
    expect_warning(empty <- genefunnel(se, list(missing = c("X", "Y")),
        BPPARAM = BiocParallel::SerialParam(), details = TRUE),
        "fewer than two")
    expect_identical(dim(empty$scores), c(0L, 3L))
    expect_identical(dimnames(empty$observed), dimnames(empty$scores))
})
