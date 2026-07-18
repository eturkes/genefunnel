# Assisted-by: OpenAI Codex.

sensitivity_contract_fixture <- function() {
    cbind(
        positive = c(A = 1, B = 2, C = 7),
        negative = c(A = 0, B = 1, C = 1),
        equal = c(A = 4, B = 4, C = 4),
        zero = c(A = -0, B = 0, C = -0),
        insufficient = c(A = 1, B = 2, C = NA_real_)
    )
}

sensitivity_contract_sets <- function(include_omitted = FALSE) {
    sets <- list(triple = c("A", "B", "C"), pair = c("A", "B"))
    if (include_omitted) {
        sets$omitted <- c("A", "absent")
    }
    sets
}

sensitivity_leaf_matrices <- function(result) {
    c(result[seq_len(6L)], result$status)
}

test_that("internal sensitivity returns the frozen aligned schema", {
    mat <- sensitivity_contract_fixture()
    serial <- BiocParallel::SerialParam()
    result <- NULL
    expect_warning(
        result <- genefunnel:::.gene_set_sensitivity(
            mat,
            sensitivity_contract_sets(include_omitted = TRUE),
            BPPARAM = serial
        ),
        "Omitting 1 gene set"
    )

    expect_false(".gene_set_sensitivity" %in% getNamespaceExports("genefunnel"))
    expect_named(
        result,
        c(
            "largest_member", "largest_absolute_delta", "largest_delta",
            "largest_delta_over_sum", "median_absolute_delta",
            "effective_size", "status"
        ),
        ignore.order = FALSE
    )
    expect_named(
        result$status,
        c("semantic", "delta", "normalized"),
        ignore.order = FALSE
    )
    leaves <- sensitivity_leaf_matrices(result)
    expect_true(all(vapply(leaves, is.matrix, logical(1L))))
    expect_true(all(vapply(
        leaves,
        function(value) identical(dimnames(value), list(
            c("triple", "pair"), colnames(mat)
        )),
        logical(1L)
    )))
    expect_identical(typeof(result$largest_member), "character")
    expect_true(all(vapply(result[2:5], typeof, character(1L)) == "double"))
    expect_identical(typeof(result$effective_size), "integer")
    expect_true(all(vapply(result$status, typeof, character(1L)) == "character"))

    expect_identical(
        unname(result$largest_member["triple", ]),
        c("B", "A", "A", "A", NA_character_)
    )
    expect_equal(
        unname(result$largest_absolute_delta["triple", ]),
        c(2.5, 1, 4, 0, NA_real_)
    )
    expect_equal(
        unname(result$largest_delta["triple", ]),
        c(2.5, -1, 4, 0, NA_real_)
    )
    expect_equal(
        unname(result$largest_delta_over_sum["triple", ]),
        c(0.25, -0.5, 1 / 3, NA_real_, NA_real_)
    )
    expect_equal(
        unname(result$median_absolute_delta["triple", ]),
        c(2.5, 1, 4, 0, NA_real_)
    )
    expect_identical(
        unname(result$effective_size["triple", ]),
        c(3L, 3L, 3L, 3L, 2L)
    )
    expect_identical(
        unname(result$status$semantic["triple", ]),
        c(rep("defined", 4L), "too_few_observed")
    )
    expect_identical(
        unname(result$status$delta["triple", ]),
        c(rep("ordinary", 4L), "not_applicable")
    )
    expect_identical(
        unname(result$status$normalized["triple", ]),
        c("ordinary", "ordinary", "ordinary", "zero_total", "not_applicable")
    )
    expect_true(all(result$effective_size["pair", ] == 2L))
    expect_true(all(is.na(result$largest_member["pair", ])))
    expect_true(all(result$status$delta["pair", ] == "not_applicable"))
    expect_identical(1 / result$largest_delta[["triple", "zero"]], Inf)
})

test_that("canonical support is stable across row order and missing encoding", {
    mat <- matrix(
        c(NA_real_, 4, NaN, 2, 1),
        ncol = 1L,
        dimnames = list(c("E", "D", "C", "B", "A"), "sample")
    )
    sets <- list(set = c("A", "B", "A", "Z", "C", "D", "E"))
    serial <- BiocParallel::SerialParam()
    observed <- genefunnel:::.gene_set_sensitivity(mat, sets, serial)
    reordered <- genefunnel:::.gene_set_sensitivity(
        mat[rev(seq_len(nrow(mat))), , drop = FALSE],
        sets,
        serial
    )
    missing_na <- mat
    missing_na[["C", "sample"]] <- NA_real_

    expect_identical(reordered, observed)
    expect_identical(
        genefunnel:::.gene_set_sensitivity(missing_na, sets, serial),
        observed
    )
    support <- sensitivity_members_reference(
        sets$set,
        rownames(mat),
        mat[, 1L]
    )
    expected <- genefunnel:::.sensitivity_cell(
        support$values,
        support$members
    )
    expect_identical(observed$effective_size[[1L]], 3L)
    expect_identical(observed$largest_member[[1L]], expected$largest_member)
    expect_equal(
        observed$largest_delta[[1L]],
        expected$largest_delta,
        tolerance = 0
    )
})

test_that("compact summaries agree with direct randomized recomputation", {
    set.seed(20260723)
    mat <- matrix(stats::rexp(12L * 25L, rate = 0.3), nrow = 12L)
    dimnames(mat) <- list(paste0("g", seq_len(nrow(mat))), paste0("s", 1:25))
    mat[sample(length(mat), 30L)] <- 0
    missing <- sample(length(mat), 35L)
    mat[missing[seq(1L, length(missing), by = 2L)]] <- NA_real_
    mat[missing[seq(2L, length(missing), by = 2L)]] <- NaN
    sets <- list(
        alpha = c(paste0("g", 1:8), "g3", "absent"),
        beta = paste0("g", 12:4),
        pair = c("g1", "g2")
    )
    observed <- genefunnel:::.gene_set_sensitivity(
        mat,
        sets,
        BiocParallel::SerialParam()
    )

    for (set_name in names(sets)) {
        for (column in seq_len(ncol(mat))) {
            support <- sensitivity_members_reference(
                sets[[set_name]],
                rownames(mat),
                mat[, column]
            )
            size <- length(support$values)
            expect_identical(
                observed$effective_size[[set_name, column]],
                as.integer(size)
            )
            if (size < 3L) {
                expect_identical(
                    observed$status$semantic[[set_name, column]],
                    "too_few_observed"
                )
                expect_true(is.na(observed$largest_member[[set_name, column]]))
                next
            }

            deltas <- sensitivity_deltas_reference(
                support$values,
                support$members
            )
            selected <- observed$largest_member[[set_name, column]]
            selected_delta <- unname(deltas[[selected]])
            scale <- max(1, max(abs(deltas)))
            expect_lte(
                abs(abs(selected_delta) - max(abs(deltas))),
                2e-12 * scale
            )
            expect_equal(
                observed$largest_absolute_delta[[set_name, column]],
                abs(selected_delta),
                tolerance = 2e-12
            )
            expect_equal(
                observed$largest_delta[[set_name, column]],
                selected_delta,
                tolerance = 2e-12
            )
            expect_equal(
                observed$largest_delta_over_sum[[set_name, column]],
                selected_delta / sum(support$values),
                tolerance = 2e-12
            )
            expect_equal(
                observed$median_absolute_delta[[set_name, column]],
                stats::median(abs(deltas)),
                tolerance = 2e-12
            )
            expect_identical(observed$status$delta[[set_name, column]], "ordinary")
            expect_identical(
                observed$status$normalized[[set_name, column]],
                "ordinary"
            )
        }
    }
})

test_that("extreme doubles use explicit ordinary and unavailable states", {
    smallest <- 2^-1074
    maximum <- .Machine$double.xmax
    mat <- cbind(
        tiny_equal = c(A = smallest, B = smallest, C = smallest),
        tiny_fraction = c(A = smallest, B = 2 * smallest, C = 7 * smallest),
        maximum_equal = c(A = maximum, B = maximum, C = maximum),
        signed_zero = c(A = -0, B = 0, C = -0)
    )
    result <- genefunnel:::.gene_set_sensitivity(
        mat,
        list(set = rownames(mat)),
        BiocParallel::SerialParam()
    )

    expect_identical(result$largest_absolute_delta[[1L, "tiny_equal"]], smallest)
    expect_identical(result$largest_delta[[1L, "tiny_equal"]], smallest)
    expect_equal(result$largest_delta_over_sum[[1L, "tiny_equal"]], 1 / 3)
    expect_identical(result$status$delta[[1L, "tiny_equal"]], "ordinary")

    expect_true(is.na(result$largest_member[[1L, "tiny_fraction"]]))
    expect_true(all(vapply(result[2:5], function(field) {
        is.na(field[[1L, "tiny_fraction"]])
    }, logical(1L))))
    expect_identical(result$status$semantic[[1L, "tiny_fraction"]], "defined")
    expect_identical(result$status$delta[[1L, "tiny_fraction"]], "unavailable")
    expect_identical(
        result$status$normalized[[1L, "tiny_fraction"]],
        "unavailable"
    )

    expect_identical(
        result$largest_absolute_delta[[1L, "maximum_equal"]],
        maximum
    )
    expect_identical(result$largest_delta[[1L, "maximum_equal"]], maximum)
    expect_equal(result$largest_delta_over_sum[[1L, "maximum_equal"]], 1 / 3)
    expect_identical(result$status$delta[[1L, "maximum_equal"]], "ordinary")

    expect_identical(result$largest_member[[1L, "signed_zero"]], "A")
    expect_identical(1 / result$largest_delta[[1L, "signed_zero"]], Inf)
    expect_true(is.na(result$largest_delta_over_sum[[1L, "signed_zero"]]))
    expect_identical(
        result$status$normalized[[1L, "signed_zero"]],
        "zero_total"
    )
})

test_that("dense sparse integer and SOCK results are byte-identical", {
    mat <- cbind(
        ordinary = c(A = 1, B = 2, C = 7, D = 0, E = 4),
        zeros = c(A = -0, B = 0, C = 0, D = 0, E = 0),
        missing = c(A = 1, B = NA, C = 3, D = NaN, E = 5),
        mixed = c(A = 9, B = 2, C = 0, D = 4, E = 1),
        partial = c(A = NA, B = 2, C = 8, D = 1, E = NA),
        equal = c(A = 6, B = 6, C = 6, D = 6, E = 6)
    )
    sets <- list(
        all = rownames(mat),
        reverse = rev(rownames(mat)),
        triple = c("A", "C", "E"),
        pair = c("A", "B")
    )
    serial <- BiocParallel::SerialParam()
    expected <- genefunnel:::.gene_set_sensitivity(mat, sets, serial)
    dense_matrix <- Matrix::Matrix(mat, sparse = FALSE)
    sparse_matrix <- Matrix::Matrix(mat, sparse = TRUE)
    integer_mat <- matrix(
        as.integer(mat),
        nrow = nrow(mat),
        dimnames = dimnames(mat)
    )

    expect_identical(
        genefunnel:::.gene_set_sensitivity(dense_matrix, sets, serial),
        expected
    )
    expect_identical(
        genefunnel:::.gene_set_sensitivity(sparse_matrix, sets, serial),
        expected
    )
    expect_identical(
        genefunnel:::.gene_set_sensitivity(integer_mat, sets, serial),
        expected
    )
    set_order <- rev(seq_along(sets))
    sample_order <- rev(seq_len(ncol(mat)))
    reordered <- genefunnel:::.gene_set_sensitivity(
        mat[, sample_order, drop = FALSE],
        sets[set_order],
        serial
    )
    for (field in names(expected)[1:6]) {
        expect_identical(
            reordered[[field]],
            expected[[field]][set_order, sample_order, drop = FALSE],
            info = field
        )
    }
    for (field in names(expected$status)) {
        expect_identical(
            reordered$status[[field]],
            expected$status[[field]][set_order, sample_order, drop = FALSE],
            info = field
        )
    }

    snow <- BiocParallel::bpstart(BiocParallel::SnowParam(
        workers = 2L,
        type = "SOCK"
    ))
    on.exit(BiocParallel::bpstop(snow), add = TRUE)
    expect_identical(
        genefunnel:::.gene_set_sensitivity(mat, sets, snow),
        expected
    )
    expect_identical(
        genefunnel:::.gene_set_sensitivity(sparse_matrix, sets, snow),
        expected
    )
    expect_true(BiocParallel::bpisup(snow))
})

test_that("canonical tie policy and permutations remain explicit", {
    mat <- cbind(
        tie = c(A = 0, B = 1, C = 1, D = NA_real_),
        unique = c(A = 2, B = 1, C = 1, D = 0)
    )
    sets <- list(
        tie_forward = c("A", "B", "C"),
        tie_reverse = c("C", "B", "A"),
        unique_forward = c("A", "B", "C", "D"),
        unique_permuted = c("D", "B", "A", "C")
    )
    result <- genefunnel:::.gene_set_sensitivity(
        mat,
        sets,
        BiocParallel::SerialParam()
    )

    expect_identical(result$largest_member[["tie_forward", "tie"]], "A")
    expect_identical(result$largest_member[["tie_reverse", "tie"]], "C")
    expect_identical(result$largest_delta[["tie_forward", "tie"]], -1)
    expect_identical(result$largest_delta[["tie_reverse", "tie"]], 1)
    expect_identical(
        result$largest_absolute_delta[["tie_forward", "tie"]],
        result$largest_absolute_delta[["tie_reverse", "tie"]]
    )
    expect_identical(
        result$median_absolute_delta[["tie_forward", "tie"]],
        result$median_absolute_delta[["tie_reverse", "tie"]]
    )
    for (field in names(result)[1:6]) {
        expect_identical(
            result[[field]][["unique_forward", "unique"]],
            result[[field]][["unique_permuted", "unique"]],
            info = field
        )
    }
})

test_that("zero-row results validation and worker errors fail closed", {
    mat <- sensitivity_contract_fixture()
    result <- NULL
    expect_warning(
        result <- genefunnel:::.gene_set_sensitivity(
            mat,
            list(omitted = c("A", "absent")),
            BiocParallel::SerialParam()
        ),
        "Omitting 1 gene set"
    )
    expect_true(all(vapply(sensitivity_leaf_matrices(result), function(value) {
        identical(dim(value), c(0L, ncol(mat))) &&
            identical(dimnames(value), list(NULL, colnames(mat)))
    }, logical(1L))))

    good <- genefunnel:::.gene_set_sensitivity(
        mat,
        sensitivity_contract_sets(),
        BiocParallel::SerialParam()
    )
    expect_silent(genefunnel:::.validate_sensitivity_chunk(good, c(2L, 5L)))

    bad <- good
    names(bad)[[1L]] <- "wrong"
    expect_error(genefunnel:::.validate_sensitivity_chunk(bad, c(2L, 5L)), "schema")
    bad <- good
    class(bad) <- "corrupt"
    expect_error(genefunnel:::.validate_sensitivity_chunk(bad, c(2L, 5L)), "schema")
    bad <- good
    bad$status$semantic[[1L]] <- "impossible"
    expect_error(genefunnel:::.validate_sensitivity_chunk(bad, c(2L, 5L)), "status")
    bad <- good
    storage.mode(bad$largest_absolute_delta) <- "integer"
    expect_error(genefunnel:::.validate_sensitivity_chunk(bad, c(2L, 5L)), "invalid")

    invalid <- mat
    invalid[[1L]] <- -1
    expect_error(
        genefunnel:::.gene_set_sensitivity(
            invalid,
            sensitivity_contract_sets(),
            BiocParallel::SerialParam()
        ),
        "negative"
    )
    expect_error(
        genefunnel:::.gene_set_sensitivity(
            mat,
            sensitivity_contract_sets(),
            BPPARAM = list()
        ),
        "BiocParallelParam"
    )

    testthat::local_mocked_bindings(
        .sensitivity_matrix_chunk = function(...) {
            stop("synthetic sensitivity failure", call. = FALSE)
        },
        .package = "genefunnel"
    )
    condition <- expect_error(genefunnel:::.gene_set_sensitivity(
        mat,
        sensitivity_contract_sets(),
        BiocParallel::SerialParam()
    ))
    expect_match(
        conditionMessage(condition),
        "GeneFunnel sensitivity calculation failed in chunk 1",
        fixed = TRUE
    )
    expect_match(conditionMessage(condition), "synthetic sensitivity failure")
})
