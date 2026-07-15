# Assisted-by: OpenAI Codex.

parallel_fixture <- function() {
    matrix(
        c(
            1, 2, 3, 4,
            0, 0, 0, 0,
            2, NA_real_, 6, 8,
            8, 4, NaN, 1,
            5, 5, 5, 5
        ),
        nrow = 4L,
        dimnames = list(
            c("A", "B", "C", "D"),
            c("ordinary", "zero", "missing", "unequal", "equal")
        )
    )
}

parallel_gene_sets <- function() {
    list(
        pair = c("A", "B"),
        partial = c("D", "C", "absent"),
        overlap = c("A", "C", "D")
    )
}

test_that("column chunks are bounded, contiguous, and nonempty", {
    serial <- genefunnel:::.column_chunk_ranges(1000000L, 1L)
    many_workers <- genefunnel:::.column_chunk_ranges(2L, 8L)
    parallel <- genefunnel:::.column_chunk_ranges(1000000L, 3L)

    expect_identical(nrow(serial), 1L)
    expect_identical(serial$first, 1L)
    expect_identical(serial$last, 1000000L)

    expect_identical(nrow(many_workers), 2L)
    expect_true(all(many_workers$last >= many_workers$first))

    expect_lte(nrow(parallel), 6L)
    expect_identical(parallel$first[[1L]], 1L)
    expect_identical(parallel$last[[nrow(parallel)]], 1000000L)
    expect_identical(
        parallel$first[-1L],
        parallel$last[-nrow(parallel)] + 1L
    )
    expect_true(all(parallel$last >= parallel$first))
})

test_that("the manager iterator streams each matrix column exactly once", {
    mat <- parallel_fixture()
    ranges <- genefunnel:::.column_chunk_ranges(ncol(mat), 2L)
    iterator <- genefunnel:::.matrix_chunk_iterator(mat, ranges)
    tasks <- lapply(seq_len(nrow(ranges)), function(index) iterator())

    expect_identical(
        vapply(tasks, `[[`, integer(1), "id"),
        seq_len(nrow(ranges))
    )
    expect_true(all(vapply(tasks, function(task) ncol(task$mat), integer(1)) > 0L))
    expect_identical(do.call(cbind, lapply(tasks, `[[`, "mat")), mat)
    expect_null(iterator())
    expect_null(iterator())

    sparse_iterator <- genefunnel:::.matrix_chunk_iterator(
        Matrix::Matrix(mat, sparse = TRUE),
        ranges
    )
    expect_s4_class(sparse_iterator()$mat, "sparseMatrix")
})

test_that("chunk assembly restores input order after out-of-order completion", {
    first <- list(
        id = 1L,
        first = 1L,
        last = 2L,
        scores = matrix(1:4, nrow = 2L)
    )
    second <- list(
        id = 2L,
        first = 3L,
        last = 4L,
        scores = matrix(5:8, nrow = 2L)
    )

    observed <- genefunnel:::.assemble_score_chunks(
        list(second, first),
        n_gene_sets = 2L,
        n_columns = 4L
    )

    expect_identical(observed, matrix(as.double(1:8), nrow = 2L))
})

test_that("serial and reusable SOCK backends return identical ordered scores", {
    mat <- parallel_fixture()
    sparse <- Matrix::Matrix(mat, sparse = TRUE)
    gene_sets <- parallel_gene_sets()
    serial <- BiocParallel::SerialParam()
    expected <- genefunnel(mat, gene_sets, BPPARAM = serial)

    snow <- BiocParallel::bpstart(BiocParallel::SnowParam(
        workers = 2L,
        type = "SOCK"
    ))
    on.exit(BiocParallel::bpstop(snow), add = TRUE)

    expect_identical(
        genefunnel(mat, gene_sets, BPPARAM = snow),
        expected
    )
    expect_identical(
        genefunnel(sparse, gene_sets, BPPARAM = snow),
        expected
    )

    reordered <- c(5L, 2L, 4L, 1L, 3L)
    expect_identical(
        genefunnel(
            mat[, reordered, drop = FALSE],
            gene_sets,
            BPPARAM = snow
        ),
        expected[, reordered, drop = FALSE]
    )
    expect_identical(
        genefunnel(mat[, 1L, drop = FALSE], gene_sets, BPPARAM = snow),
        expected[, 1L, drop = FALSE]
    )
    expect_true(BiocParallel::bpisup(snow))

    invalid_task <- list(
        id = 7L,
        first = 5L,
        last = 5L,
        mat = matrix(
            c(Inf, 1),
            nrow = 2L,
            dimnames = list(c("A", "B"), "equal")
        )
    )
    condition <- expect_error(BiocParallel::bpiterate(
        list(invalid_task),
        genefunnel:::.score_matrix_task,
        gene_indices = list(c(1L, 2L)),
        storage = "dense",
        BPPARAM = snow
    ))
    message <- conditionMessage(condition)
    expect_match(message, "GeneFunnel scoring failed in chunk 7", fixed = TRUE)
    expect_match(message, "matrix column 5, sample \"equal\"", fixed = TRUE)
    expect_match(message, "infinite value at matrix row 1", fixed = TRUE)
})

test_that("worker failures report chunk positions and sample names", {
    mat <- parallel_fixture()
    testthat::local_mocked_bindings(
        .score_matrix_chunk = function(...) {
            stop("synthetic native failure", call. = FALSE)
        },
        .package = "genefunnel"
    )

    condition <- expect_error(genefunnel(
        mat,
        parallel_gene_sets(),
        BPPARAM = BiocParallel::SerialParam()
    ))
    message <- conditionMessage(condition)

    expect_match(message, "GeneFunnel scoring failed in chunk 1", fixed = TRUE)
    expect_match(message, "matrix columns 1-5", fixed = TRUE)
    expect_match(message, 'samples "ordinary" through "equal"', fixed = TRUE)
    expect_match(message, "synthetic native failure", fixed = TRUE)
})
