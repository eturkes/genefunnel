# Assisted-by: OpenAI Codex.

component_contract_fixture <- function() {
    matrix(
        c(
            2, 2, NA_real_,
            6, 2, NA_real_,
            4 / 3, 4 / 3, 4 / 3,
            4, 0, NA_real_,
            0, 0, NA_real_,
            1, 2, NA_real_,
            1, NA_real_, NaN
        ),
        nrow = 3L,
        dimnames = list(
            c("A", "B", "C"),
            c(
                "equal_pair",
                "unequal_pair",
                "equal_triple",
                "one_positive",
                "zero_total",
                "partial",
                "too_few"
            )
        )
    )
}

component_contract_sets <- function(include_omitted = FALSE) {
    sets <- list(
        pair = c("A", "B"),
        triple = c("A", "B", "C")
    )
    if (include_omitted) {
        sets$omitted <- c("A", "absent")
    }
    sets
}

component_leaf_matrices <- function(result) {
    c(
        result[c(
            "score",
            "observed_sum",
            "penalty",
            "balance",
            "effective_size",
            "observed_fraction"
        )],
        result$status,
        unlist(result$scaled, recursive = FALSE, use.names = FALSE)
    )
}

expect_component_representation <- function(
    result,
    field,
    column,
    expected,
    effective_size,
    row = 1L
) {
    ordinary <- scaled_ref_as_double(expected)
    observed <- result[[field]][[row, column]]
    status <- result$status[[field]][[row, column]]
    mantissa <- result$scaled[[field]]$mantissa[[row, column]]
    exponent <- result$scaled[[field]]$exponent[[row, column]]

    if (is.na(ordinary)) {
        pair <- scaled_ref_pair(expected)
        allowance <- 8 * 64 * (effective_size + 1) *
            (.Machine$double.eps / 2)
        expect_true(is.na(observed), info = paste(field, column))
        expect_identical(status, "scaled", info = paste(field, column))
        expect_equal(
            mantissa,
            unname(pair[["mantissa"]]),
            tolerance = allowance,
            info = paste(field, column)
        )
        expect_identical(
            exponent,
            as.integer(pair[["exponent"]]),
            info = paste(field, column)
        )
    } else {
        expect_identical(status, "ordinary", info = paste(field, column))
        expect_equal(
            observed,
            ordinary,
            tolerance = 4e-14,
            info = paste(field, column)
        )
        expect_true(is.na(mantissa), info = paste(field, column))
        expect_true(is.na(exponent), info = paste(field, column))
    }
}

component_result_as_scaled <- function(result, field, row, column) {
    status <- result$status[[field]][[row, column]]
    if (status == "ordinary") {
        return(scaled_ref_from_double(result[[field]][[row, column]]))
    }
    stopifnot(status == "scaled")
    list(
        hi = result$scaled[[field]]$mantissa[[row, column]],
        lo = 0,
        exponent = result$scaled[[field]]$exponent[[row, column]]
    )
}

test_that("component API returns the locked aligned schema", {
    mat <- component_contract_fixture()
    gene_sets <- component_contract_sets(include_omitted = TRUE)
    serial <- BiocParallel::SerialParam()

    result <- NULL
    expect_warning(
        result <- genefunnel_components(mat, gene_sets, BPPARAM = serial),
        "Omitting 1 gene set"
    )
    scores <- NULL
    expect_warning(
        scores <- genefunnel(mat, gene_sets, BPPARAM = serial),
        "Omitting 1 gene set"
    )

    expect_named(result, c(
        "score",
        "observed_sum",
        "penalty",
        "balance",
        "effective_size",
        "observed_fraction",
        "status",
        "scaled"
    ), ignore.order = FALSE)
    expect_named(result$status, c(
        "semantic",
        "observed_sum",
        "penalty",
        "balance",
        "conditioning"
    ), ignore.order = FALSE)
    expect_named(
        result$scaled,
        c("observed_sum", "penalty", "balance"),
        ignore.order = FALSE
    )
    for (field in names(result$scaled)) {
        expect_named(
            result$scaled[[field]],
            c("mantissa", "exponent"),
            ignore.order = FALSE
        )
    }

    leaves <- component_leaf_matrices(result)
    expect_true(all(vapply(leaves, is.matrix, logical(1L))))
    expect_true(all(vapply(
        leaves,
        function(value) identical(dimnames(value), dimnames(scores)),
        logical(1L)
    )))
    expect_identical(result$score, scores)
    expect_true(all(vapply(
        result[c("score", "observed_sum", "penalty", "balance",
            "observed_fraction")],
        typeof,
        character(1L)
    ) == "double"))
    expect_identical(typeof(result$effective_size), "integer")
    expect_true(all(vapply(result$status, typeof, character(1L)) == "character"))
    expect_true(all(vapply(
        lapply(result$scaled, `[[`, "mantissa"),
        typeof,
        character(1L)
    ) == "double"))
    expect_true(all(vapply(
        lapply(result$scaled, `[[`, "exponent"),
        typeof,
        character(1L)
    ) == "integer"))
})

test_that("components implement the committed interpretation and edge cases", {
    result <- genefunnel_components(
        component_contract_fixture(),
        component_contract_sets(),
        BPPARAM = BiocParallel::SerialParam()
    )

    expect_identical(result$score[["pair", "equal_pair"]], 4)
    expect_identical(result$observed_sum[["pair", "equal_pair"]], 4)
    expect_identical(result$balance[["pair", "equal_pair"]], 1)
    expect_identical(result$effective_size[["pair", "equal_pair"]], 2L)

    expect_identical(result$score[["pair", "unequal_pair"]], 4)
    expect_identical(result$observed_sum[["pair", "unequal_pair"]], 8)
    expect_identical(result$penalty[["pair", "unequal_pair"]], 4)
    expect_identical(result$balance[["pair", "unequal_pair"]], 0.5)

    expect_equal(result$score[["triple", "equal_triple"]], 4)
    expect_equal(result$observed_sum[["triple", "equal_triple"]], 4)
    expect_identical(result$balance[["triple", "equal_triple"]], 1)
    expect_identical(
        result$effective_size[["triple", "equal_triple"]],
        3L
    )

    expect_identical(result$score[["pair", "one_positive"]], 0)
    expect_identical(result$penalty[["pair", "one_positive"]], 4)
    expect_identical(result$balance[["pair", "one_positive"]], 0)
    expect_identical(
        result$status$conditioning[["pair", "one_positive"]],
        "ill_conditioned"
    )

    expect_identical(result$score[["pair", "zero_total"]], 0)
    expect_identical(result$observed_sum[["pair", "zero_total"]], 0)
    expect_identical(result$penalty[["pair", "zero_total"]], 0)
    expect_true(is.na(result$balance[["pair", "zero_total"]]))
    expect_identical(
        result$status$semantic[["pair", "zero_total"]],
        "zero_total"
    )
    expect_identical(
        result$status$balance[["pair", "zero_total"]],
        "unavailable"
    )
    expect_identical(
        result$status$conditioning[["pair", "zero_total"]],
        "not_applicable"
    )

    expect_identical(result$score[["triple", "partial"]], 2)
    expect_identical(result$observed_sum[["triple", "partial"]], 3)
    expect_identical(result$penalty[["triple", "partial"]], 1)
    expect_equal(result$balance[["triple", "partial"]], 2 / 3)
    expect_identical(result$effective_size[["triple", "partial"]], 2L)
    expect_equal(result$observed_fraction[["triple", "partial"]], 2 / 3)

    expect_true(is.na(result$score[["triple", "too_few"]]))
    expect_identical(result$observed_sum[["triple", "too_few"]], 1)
    expect_true(is.na(result$penalty[["triple", "too_few"]]))
    expect_true(is.na(result$balance[["triple", "too_few"]]))
    expect_identical(result$effective_size[["triple", "too_few"]], 1L)
    expect_equal(result$observed_fraction[["triple", "too_few"]], 1 / 3)
    expect_identical(
        result$status$semantic[["triple", "too_few"]],
        "too_few_observed"
    )
    expect_identical(
        result$status$penalty[["triple", "too_few"]],
        "unavailable"
    )
    expect_identical(
        result$status$conditioning[["triple", "too_few"]],
        "not_applicable"
    )
})

test_that("extreme diagnostics use ordinary or scaled representation honestly", {
    maximum <- .Machine$double.xmax
    minimum <- 2^-1074
    mat <- matrix(
        c(
            maximum, maximum / 2, NA_real_, NA_real_,
            maximum, maximum, 0, 0,
            2^1023, 2^-1022, 0, NA_real_,
            minimum, minimum, NA_real_, NaN
        ),
        nrow = 4L,
        dimnames = list(
            LETTERS[1:4],
            c(
                "sum_overflow",
                "penalty_overflow",
                "balance_underflow",
                "subnormal"
            )
        )
    )
    gene_sets <- list(all = LETTERS[1:4])
    result <- genefunnel_components(
        mat,
        gene_sets,
        BPPARAM = BiocParallel::SerialParam()
    )
    expected_score <- genefunnel(
        mat,
        gene_sets,
        BPPARAM = BiocParallel::SerialParam()
    )

    expect_identical(result$score, expected_score)
    for (column in seq_len(ncol(mat))) {
        expected <- reference_scaled_components(mat[, column])
        expect_identical(
            result$status$semantic[[1L, column]],
            expected$semantic
        )
        expect_identical(
            result$status$conditioning[[1L, column]],
            expected$conditioning
        )
        expect_identical(
            result$effective_size[[1L, column]],
            expected$effective_size
        )
        expect_equal(
            result$observed_fraction[[1L, column]],
            expected$observed_fraction
        )
        for (field in c("observed_sum", "penalty", "balance")) {
            expect_component_representation(
                result,
                field,
                column,
                expected[[field]],
                expected$effective_size
            )
        }
    }

    expect_identical(
        result$status$observed_sum[[1L, "sum_overflow"]],
        "scaled"
    )
    expect_identical(
        result$status$penalty[[1L, "penalty_overflow"]],
        "scaled"
    )
    expect_identical(
        result$status$balance[[1L, "balance_underflow"]],
        "scaled"
    )
    expect_identical(
        result$scaled$balance$exponent[[1L, "balance_underflow"]],
        -2044L
    )

    for (column in seq_len(ncol(mat))) {
        if (result$score[[1L, column]] <= 0) {
            next
        }
        observed_sum <- component_result_as_scaled(
            result,
            "observed_sum",
            1L,
            column
        )
        penalty <- component_result_as_scaled(
            result,
            "penalty",
            1L,
            column
        )
        balance <- component_result_as_scaled(
            result,
            "balance",
            1L,
            column
        )
        score <- scaled_ref_from_double(result$score[[1L, column]])
        allowance <- 8 * 64 *
            (result$effective_size[[1L, column]] + 1) *
            (.Machine$double.eps / 2)
        expect_lt(
            scaled_ref_relative_error(
                scaled_ref_add(score, penalty),
                observed_sum
            ),
            allowance
        )
        expect_lt(
            scaled_ref_relative_error(
                scaled_ref_multiply(observed_sum, balance),
                score
            ),
            allowance
        )
    }
})

test_that("ordinary component cells match the oracle and safe identities", {
    set.seed(18072027)
    n_features <- 24L
    n_samples <- 20L
    features <- sprintf("feature_%02d", seq_len(n_features))
    mat <- matrix(
        stats::rexp(n_features * n_samples),
        nrow = n_features,
        dimnames = list(
            features,
            sprintf("sample_%02d", seq_len(n_samples))
        )
    )
    mat[stats::runif(length(mat)) < 0.35] <- 0
    mat <- sweep(
        mat,
        2L,
        10^sample(-80:80, n_samples, replace = TRUE),
        `*`
    )
    mat[stats::runif(length(mat)) < 0.12] <- NA_real_
    gene_sets <- list(
        small = features[c(8L, 2L, 17L)],
        medium = features[c(24L, 1L, 12L, 6L, 20L, 9L, 3L, 15L)],
        large = sample(features, 18L)
    )
    result <- genefunnel_components(
        mat,
        gene_sets,
        BPPARAM = BiocParallel::SerialParam()
    )
    expected_score <- genefunnel(
        mat,
        gene_sets,
        BPPARAM = BiocParallel::SerialParam()
    )

    expect_identical(result$score, expected_score)
    ordinary_balance <- result$status$balance == "ordinary"
    expect_true(all(result$balance[ordinary_balance] >= 0))
    expect_true(all(result$balance[ordinary_balance] <= 1))
    memberships <- reference_memberships(gene_sets, rownames(mat))
    for (set_index in seq_along(memberships)) {
        for (sample_index in seq_len(ncol(mat))) {
            values <- mat[memberships[[set_index]], sample_index]
            expected <- reference_scaled_components(values)
            expect_identical(
                result$status$semantic[[set_index, sample_index]],
                expected$semantic
            )
            expect_identical(
                result$status$conditioning[[set_index, sample_index]],
                expected$conditioning
            )
            expect_identical(
                result$effective_size[[set_index, sample_index]],
                expected$effective_size
            )
            expect_equal(
                result$observed_fraction[[set_index, sample_index]],
                expected$observed_fraction
            )
            expect_component_representation(
                result,
                "observed_sum",
                sample_index,
                expected$observed_sum,
                expected$effective_size,
                row = set_index
            )
            if (expected$semantic != "scoreable") {
                next
            }
            expect_component_representation(
                result,
                "penalty",
                sample_index,
                expected$penalty,
                expected$effective_size,
                row = set_index
            )
            expect_component_representation(
                result,
                "balance",
                sample_index,
                expected$balance,
                expected$effective_size,
                row = set_index
            )
            if (expected$conditioning == "safe") {
                n <- expected$effective_size
                tolerance <- 8 * 64 * (n + 1) *
                    (.Machine$double.eps / 2) * max(
                        result$observed_sum[[set_index, sample_index]],
                        result$penalty[[set_index, sample_index]],
                        result$score[[set_index, sample_index]]
                    )
                expect_lte(
                    abs(
                        result$observed_sum[[set_index, sample_index]] -
                            result$penalty[[set_index, sample_index]] -
                            result$score[[set_index, sample_index]]
                    ),
                    tolerance
                )
                expect_lte(
                    abs(
                        result$observed_sum[[set_index, sample_index]] *
                            result$balance[[set_index, sample_index]] -
                            result$score[[set_index, sample_index]]
                    ),
                    tolerance
                )
            }
        }
    }
})

test_that("member permutations and signed zero preserve component facts", {
    values <- c(-0, 2^-1074, 0.25, 2, 7, 19)
    mat <- matrix(
        values,
        dimnames = list(LETTERS[seq_along(values)], "sample")
    )
    serial <- BiocParallel::SerialParam()
    expected <- genefunnel_components(
        mat,
        list(set = rownames(mat)),
        BPPARAM = serial
    )

    set.seed(18072028)
    for (index in seq_len(24L)) {
        members <- sample(rownames(mat))
        observed <- genefunnel_components(
            mat,
            list(set = members),
            BPPARAM = serial
        )
        expect_equal(observed$score, expected$score, tolerance = 4e-14)
        expect_equal(
            observed$observed_sum,
            expected$observed_sum,
            tolerance = 4e-14
        )
        expect_equal(observed$penalty, expected$penalty, tolerance = 4e-14)
        expect_equal(observed$balance, expected$balance, tolerance = 4e-14)
        expect_identical(observed$effective_size, expected$effective_size)
        expect_identical(
            observed$observed_fraction,
            expected$observed_fraction
        )
        expect_identical(observed$status, expected$status)
    }

    zero <- genefunnel_components(
        matrix(
            c(-0, 0),
            dimnames = list(c("A", "B"), "sample")
        ),
        list(set = c("A", "B")),
        BPPARAM = serial
    )
    expect_identical(1 / zero$score[[1L]], Inf)
    expect_identical(1 / zero$observed_sum[[1L]], Inf)
    expect_identical(1 / zero$penalty[[1L]], Inf)
    expect_identical(zero$status$semantic[[1L]], "zero_total")
})

test_that("component outputs are dense-sparse and serial-SOCK identical", {
    mat <- component_contract_fixture()
    sparse <- Matrix::Matrix(mat, sparse = TRUE)
    gene_sets <- component_contract_sets()
    serial <- BiocParallel::SerialParam()
    expected <- genefunnel_components(mat, gene_sets, BPPARAM = serial)

    expect_identical(
        genefunnel_components(sparse, gene_sets, BPPARAM = serial),
        expected
    )

    snow <- BiocParallel::bpstart(BiocParallel::SnowParam(
        workers = 2L,
        type = "SOCK"
    ))
    on.exit(BiocParallel::bpstop(snow), add = TRUE)
    expect_identical(
        genefunnel_components(mat, gene_sets, BPPARAM = snow),
        expected
    )
    expect_identical(
        genefunnel_components(sparse, gene_sets, BPPARAM = snow),
        expected
    )
})

test_that("component API preserves zero-row results and core validation", {
    mat <- component_contract_fixture()
    result <- NULL
    expect_warning(
        result <- genefunnel_components(
            mat,
            list(omitted = c("A", "absent")),
            BPPARAM = BiocParallel::SerialParam()
        ),
        "Omitting 1 gene set"
    )
    scores <- NULL
    expect_warning(
        scores <- genefunnel(
            mat,
            list(omitted = c("A", "absent")),
            BPPARAM = BiocParallel::SerialParam()
        ),
        "Omitting 1 gene set"
    )

    expect_true(all(vapply(
        component_leaf_matrices(result),
        function(value) identical(dim(value), c(0L, ncol(mat))),
        logical(1L)
    )))
    expect_true(all(vapply(
        component_leaf_matrices(result),
        function(value) identical(dimnames(value), dimnames(scores)),
        logical(1L)
    )))

    invalid <- mat
    invalid[[1L, 1L]] <- -1
    expect_error(
        genefunnel_components(
            invalid,
            component_contract_sets(),
            BPPARAM = BiocParallel::SerialParam()
        ),
        "negative"
    )
})

test_that("component native boundaries fail closed on malformed calls", {
    dense <- matrix(c(1, 2), nrow = 2L)
    sparse <- Matrix::Matrix(dense, sparse = TRUE)
    for (native_call in list(
        function() genefunnel:::calculateComponentsDense(
            dense,
            list(integer())
        ),
        function() genefunnel:::calculateComponentsSparse(
            sparse,
            list(integer())
        )
    )) {
        expect_error(native_call(), "support counts are invalid")
    }
    for (native_call in list(
        function() genefunnel:::calculateComponentsDense(
            dense,
            list(c(1L, 3L))
        ),
        function() genefunnel:::calculateComponentsSparse(
            sparse,
            list(c(1L, 3L))
        )
    )) {
        expect_error(native_call(), "gene-set index is invalid")
    }
})
