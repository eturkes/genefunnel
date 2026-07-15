# Assisted-by: OpenAI Codex.

random_score_fixture <- function(case_id) {
    set.seed(110000L + case_id)
    n_features <- sample.int(23L, 1L) + 1L
    n_samples <- sample.int(7L, 1L)
    features <- sprintf("feature_%02d", seq_len(n_features))

    mat <- matrix(
        stats::rexp(n_features * n_samples),
        nrow = n_features,
        dimnames = list(features, sprintf("sample_%02d", seq_len(n_samples)))
    )
    mat[stats::runif(length(mat)) < 0.45] <- 0
    mat <- sweep(
        mat,
        2L,
        10^sample(-12:12, n_samples, replace = TRUE),
        `*`
    )
    missing <- which(stats::runif(length(mat)) < 0.16)
    mat[missing] <- sample(
        c(NA_real_, NaN),
        length(missing),
        replace = TRUE
    )

    n_sets <- sample.int(8L, 1L)
    gene_sets <- setNames(lapply(seq_len(n_sets), function(set_index) {
        matched_size <- sample.int(min(n_features, 10L) - 1L, 1L) + 1L
        matched <- sample(features, matched_size)
        unmatched_size <- sample.int(3L, 1L) - 1L
        unmatched <- if (unmatched_size == 0L) {
            character()
        } else {
            sprintf(
                "absent_%02d_%02d_%02d",
                case_id,
                set_index,
                seq_len(unmatched_size)
            )
        }
        duplicate_size <- sample(0:2, 1L)
        duplicates <- if (duplicate_size == 0L) {
            character()
        } else {
            sample(matched, duplicate_size, replace = TRUE)
        }
        c(matched, unmatched, duplicates)
    }), sprintf("set_%02d", seq_len(n_sets)))

    list(mat = mat, gene_sets = gene_sets)
}

test_that("deterministic randomized dense and sparse scores match the oracle", {
    serial <- BiocParallel::SerialParam()

    for (case_id in seq_len(32L)) {
        fixture <- random_score_fixture(case_id)
        expected <- reference_scores(fixture$mat, fixture$gene_sets)
        dense <- genefunnel(
            fixture$mat,
            fixture$gene_sets,
            BPPARAM = serial
        )
        sparse <- genefunnel(
            Matrix::Matrix(fixture$mat, sparse = TRUE),
            fixture$gene_sets,
            BPPARAM = serial
        )

        expect_equal(dense, expected, tolerance = 5e-13)
        expect_equal(sparse, expected, tolerance = 5e-13)
        expect_identical(dimnames(dense), dimnames(expected))
        expect_identical(dimnames(sparse), dimnames(expected))
    }
})

test_that("native scores satisfy non-negativity and observed-sum bounds", {
    set.seed(220714L)
    n_features <- 64L
    n_samples <- 96L
    features <- sprintf("feature_%02d", seq_len(n_features))
    mat <- matrix(
        stats::rexp(n_features * n_samples),
        nrow = n_features,
        dimnames = list(features, sprintf("sample_%02d", seq_len(n_samples)))
    )
    mat[stats::runif(length(mat)) < 0.5] <- 0
    mat <- sweep(mat, 2L, 10^seq(-100, 100, length.out = n_samples), `*`)
    mat[stats::runif(length(mat)) < 0.08] <- NA_real_
    set_sizes <- c(2L, 3L, 4L, 8L, 16L, 32L, 64L)
    gene_sets <- setNames(
        lapply(set_sizes, function(size) sample(features, size)),
        paste0("size_", set_sizes)
    )

    observed <- genefunnel(
        mat,
        gene_sets,
        BPPARAM = BiocParallel::SerialParam()
    )
    observed_sums <- reference_cells(mat, gene_sets, function(values) {
        values <- values[!is.na(values)]
        if (length(values) < 2L) NA_real_ else sum(values)
    })
    available <- !is.na(observed)
    margin <- 128 * .Machine$double.eps * pmax(
        observed_sums[available],
        .Machine$double.xmin
    )

    expect_identical(is.na(observed), is.na(observed_sums))
    expect_true(all(is.finite(observed[available])))
    expect_true(all(observed[available] >= 0))
    expect_true(all(
        observed[available] <= observed_sums[available] + margin
    ))
})
