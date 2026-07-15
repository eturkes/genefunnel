# Assisted-by: OpenAI Codex.

valid_contract_matrix <- function() {
    matrix(
        c(1, 2, 3, 4),
        nrow = 2L,
        dimnames = list(c("A", "B"), c("sample_1", "sample_2"))
    )
}

score_contract_matrix <- function(mat, gene_sets = list(set = c("A", "B"))) {
    genefunnel(
        mat,
        gene_sets,
        BPPARAM = BiocParallel::SerialParam()
    )
}

test_that("supported base and Matrix numeric classes are accepted", {
    base_double <- valid_contract_matrix()
    base_integer <- base_double
    storage.mode(base_integer) <- "integer"
    matrix_dense <- Matrix::Matrix(base_double, sparse = FALSE)
    matrix_sparse <- Matrix::Matrix(base_double, sparse = TRUE)

    observed <- lapply(
        list(base_double, base_integer, matrix_dense, matrix_sparse),
        score_contract_matrix
    )

    expect_true(all(vapply(observed, is.matrix, logical(1))))
    expect_true(all(vapply(observed, typeof, character(1)) == "double"))
    for (candidate in observed[-1L]) {
        expect_identical(candidate, observed[[1L]])
    }
})

test_that("coercion-prone and non-numeric matrix inputs are rejected", {
    logical_mat <- matrix(
        TRUE,
        nrow = 2L,
        dimnames = list(c("A", "B"), "sample")
    )
    character_mat <- logical_mat
    storage.mode(character_mat) <- "character"
    complex_mat <- logical_mat
    storage.mode(complex_mat) <- "complex"
    logical_sparse <- Matrix::Matrix(logical_mat, sparse = TRUE)
    data_frame <- data.frame(sample = c(1, 2), row.names = c("A", "B"))
    three_dimensional <- array(
        1,
        dim = c(2L, 1L, 1L),
        dimnames = list(c("A", "B"), "sample", "layer")
    )

    for (candidate in list(
        logical_mat,
        character_mat,
        complex_mat,
        logical_sparse,
        data_frame,
        three_dimensional
    )) {
        expect_error(
            score_contract_matrix(candidate),
            "numeric or integer matrix"
        )
    }
})

test_that("matrix dimensions and feature identifiers are strict", {
    zero_rows <- matrix(
        numeric(),
        nrow = 0L,
        ncol = 1L,
        dimnames = list(character(), "sample")
    )
    zero_columns <- matrix(
        numeric(),
        nrow = 2L,
        ncol = 0L,
        dimnames = list(c("A", "B"), character())
    )
    expect_error(score_contract_matrix(zero_rows), "at least one row")
    expect_error(score_contract_matrix(zero_columns), "at least one column")

    missing_names <- valid_contract_matrix()
    rownames(missing_names) <- NULL
    expect_error(score_contract_matrix(missing_names), "row names")

    invalid_names <- list(
        missing = c("A", NA_character_),
        empty = c("A", ""),
        duplicated = c("A", "A")
    )
    for (kind in names(invalid_names)) {
        candidate <- valid_contract_matrix()
        rownames(candidate) <- invalid_names[[kind]]
        expect_error(score_contract_matrix(candidate), kind)
    }
})

test_that("negative and infinite values are rejected before scoring", {
    invalid_values <- list(
        negative = -1,
        positive_infinity = Inf,
        negative_infinity = -Inf
    )

    for (kind in names(invalid_values)) {
        candidate <- valid_contract_matrix()
        candidate[1L, 1L] <- invalid_values[[kind]]
        expect_error(
            score_contract_matrix(candidate),
            if (kind == "negative") "negative" else "infinite"
        )

        sparse_candidate <- Matrix::sparseMatrix(
            i = 1L,
            j = 1L,
            x = invalid_values[[kind]],
            dims = c(2L, 2L),
            dimnames = dimnames(valid_contract_matrix())
        )
        expect_error(
            score_contract_matrix(sparse_candidate),
            if (kind == "negative") "negative" else "infinite"
        )
    }

    missing <- valid_contract_matrix()
    missing[1L, 1L] <- NA_real_
    missing[2L, 2L] <- NaN
    expect_no_error(score_contract_matrix(missing))

    sparse_missing <- Matrix::sparseMatrix(
        i = c(1L, 2L),
        j = c(1L, 2L),
        x = c(NA_real_, NaN),
        dims = c(2L, 2L),
        dimnames = dimnames(valid_contract_matrix())
    )
    expect_no_error(score_contract_matrix(sparse_missing))
})

test_that("Matrix validation uses represented rather than raw stored values", {
    duplicate_triplet <- Matrix::sparseMatrix(
        i = c(1L, 1L, 2L),
        j = c(1L, 1L, 1L),
        x = c(-1, 2, 1),
        dims = c(2L, 1L),
        dimnames = list(c("A", "B"), "sample"),
        repr = "T"
    )
    triangular_dense <- methods::new(
        "dtrMatrix",
        Dim = c(2L, 2L),
        Dimnames = list(c("A", "B"), c("sample_1", "sample_2")),
        x = c(1, -99, 2, 1),
        uplo = "U",
        diag = "N"
    )

    expect_identical(
        score_contract_matrix(duplicate_triplet),
        matrix(2, nrow = 1L, dimnames = list("set", "sample"))
    )
    expect_identical(
        score_contract_matrix(triangular_dense),
        matrix(
            c(0, 2),
            nrow = 1L,
            dimnames = list("set", c("sample_1", "sample_2"))
        )
    )
})

test_that("gene-set structures use one strict shared contract", {
    mat <- valid_contract_matrix()
    coverage_call <- function(gene_sets) {
        gene_set_coverage(gene_sets, rownames(mat))
    }
    scoring_call <- function(gene_sets) {
        score_contract_matrix(mat, gene_sets)
    }
    calls <- list(coverage_call, scoring_call)

    malformed <- list(
        empty_list = list(),
        data_frame = data.frame(set = c("A", "B")),
        unnamed = list(c("A", "B")),
        missing_name = stats::setNames(list(c("A", "B")), NA_character_),
        empty_name = stats::setNames(list(c("A", "B")), ""),
        duplicated_name = list(set = c("A", "B"), set = c("A", "B")),
        non_character = list(set = 1:2),
        matrix_member = list(set = matrix(c("A", "B"), nrow = 1L)),
        missing_member = list(set = c("A", NA_character_)),
        empty_member = list(set = c("A", ""))
    )
    patterns <- c(
        empty_list = "at least one",
        data_frame = "named list",
        unnamed = "names",
        missing_name = "missing",
        empty_name = "empty",
        duplicated_name = "duplicated",
        non_character = "character vector",
        matrix_member = "character vector",
        missing_member = "missing",
        empty_member = "empty"
    )

    for (call in calls) {
        for (kind in names(malformed)) {
            expect_error(call(malformed[[kind]]), patterns[[kind]])
        }
    }
})

test_that("coverage feature identifiers are validated without coercion", {
    gene_sets <- list(set = c("A", "B"))
    malformed <- list(
        factor(c("A", "B")),
        matrix(c("A", "B"), nrow = 1L),
        c("A", NA_character_),
        c("A", ""),
        c("A", "A")
    )

    for (features in malformed) {
        expect_error(gene_set_coverage(gene_sets, features), "features")
    }
    expect_no_error(gene_set_coverage(gene_sets, character()))
})

test_that("BPPARAM has a qualified default and rejects non-backends", {
    expect_identical(
        formals(genefunnel)$BPPARAM,
        quote(BiocParallel::bpparam())
    )

    for (candidate in list(NULL, "serial", list(BiocParallel::SerialParam()))) {
        expect_error(
            genefunnel(
                valid_contract_matrix(),
                list(set = c("A", "B")),
                candidate
            ),
            "BiocParallelParam"
        )
    }
})

test_that("the registered BiocParallel backend is usable by default", {
    previous <- BiocParallel::registered()
    on.exit({
        for (parameter in rev(previous)) {
            BiocParallel::register(parameter)
        }
    }, add = TRUE)
    BiocParallel::register(BiocParallel::SerialParam(), default = TRUE)

    observed <- genefunnel(
        valid_contract_matrix(),
        list(set = c("A", "B"))
    )
    expect_identical(unname(observed), matrix(c(2, 6), nrow = 1L))
})
