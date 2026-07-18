# Assisted-by: OpenAI Codex.

catalogue_fixture <- function() {
    list(
        full = c("A", "B", "A"),
        partial = c("C", "absent"),
        overlap = c("B", "C", "D"),
        empty = character()
    )
}

catalogue_matrix <- function() {
    matrix(
        c(
            4, 4, 2, 0,
            0, 0, 2, 8,
            4, NA_real_, 2, 1,
            1, 3, NaN, 5,
            0, 0, 0, 0
        ),
        nrow = 4L,
        dimnames = list(
            c("A", "B", "C", "D"),
            c("equal", "unequal", "missing", "mixed", "zero")
        )
    )
}

test_that("compiled catalogue stores the fixed canonical state", {
    gene_sets <- catalogue_fixture()
    features <- c(AliasA = "A", AliasB = "B", AliasC = "C", AliasD = "D")
    observed <- genefunnel:::.compile_gene_sets(gene_sets, features)
    state <- unclass(observed)

    expect_s3_class(observed, "genefunnel_catalogue", exact = TRUE)
    expect_named(
        state,
        c(
            "schema_version", "formula_version", "features", "members",
            "duplicate_member_count", "member_index", "indices", "coverage",
            "row_adjacency", "fingerprints"
        )
    )
    expect_identical(state$schema_version, 1L)
    expect_identical(state$formula_version, "GF-1")
    expect_identical(state$features, unname(features))
    expect_identical(state$members$full, c("A", "B"))
    expect_identical(state$duplicate_member_count, c(1L, 0L, 0L, 0L))
    expect_identical(state$member_index$full, c(1L, 2L))
    expect_identical(state$member_index$partial, c(3L, 0L))
    expect_identical(state$indices$partial, 3L)
    expect_identical(
        state$coverage,
        gene_set_coverage(gene_sets, unname(features))
    )
    expect_identical(state$row_adjacency$offsets, c(0L, 1L, 3L, 5L, 6L))
    expect_identical(state$row_adjacency$set_indices, c(1L, 1L, 3L, 2L, 3L, 3L))
    expect_identical(
        names(state$fingerprints),
        c("algorithm", "encoding", "catalogue", "features", "content")
    )
    expect_identical(unname(state$fingerprints[1:2]), c("SHA-256", "GFCAT-1"))
    expect_identical(
        unname(state$fingerprints[c("catalogue", "features", "content")]),
        c(
            "73b2135ea1a239e87967410718a14f0b6688b8689d90bf7f18fa408226eae89c",
            "58244d86ca48c66d23a1a558a34638fd6b3f014e661528339bf49d23b116463e",
            "e57cc2f8e0ead5d737f7c361ce9733dc944795d7794178bbf17275901a3f418e"
        )
    )
    expect_invisible(genefunnel:::.validate_compiled_catalogue(observed))

    gene_sets$full[[1L]] <- "changed"
    features[[1L]] <- "changed"
    expect_identical(state$members$full, c("A", "B"))
    expect_identical(state$features, c("A", "B", "C", "D"))
})

test_that("compiled constructor retains the named-list validation boundary", {
    expect_error(
        genefunnel:::.compile_gene_sets(unname(catalogue_fixture()), LETTERS[1:4]),
        "must have names"
    )
    expect_error(
        genefunnel:::.compile_gene_sets(catalogue_fixture(), c("A", "A")),
        "duplicated identifiers"
    )
    spoofed <- structure(c("A", "B"), class = "spoof")
    expect_error(
        genefunnel:::.compile_gene_sets(list(set = spoofed), c("A", "B")),
        "unclassed character vector"
    )
})

test_that("compiled scoring equals authoritative dense and sparse calls", {
    mat <- catalogue_matrix()
    gene_sets <- catalogue_fixture()
    compiled <- genefunnel:::.compile_gene_sets(gene_sets, rownames(mat))
    serial <- BiocParallel::SerialParam()

    expect_warning(
        expected <- genefunnel(mat, gene_sets, BPPARAM = serial),
        "Omitting 2 gene sets"
    )
    expect_warning(
        dense <- genefunnel:::.genefunnel_compiled(mat, compiled, serial),
        "Omitting 2 gene sets"
    )
    expect_warning(
        sparse <- genefunnel:::.genefunnel_compiled(
            Matrix::Matrix(mat, sparse = TRUE),
            compiled,
            serial
        ),
        "Omitting 2 gene sets"
    )
    expect_identical(dense, expected)
    expect_identical(sparse, expected)

    no_score <- genefunnel:::.compile_gene_sets(
        list(empty = character(), one = "A"),
        rownames(mat)
    )
    expect_warning(
        empty <- genefunnel:::.genefunnel_compiled(mat, no_score, serial),
        "Omitting 2 gene sets"
    )
    expect_identical(dim(empty), c(0L, ncol(mat)))
    expect_null(rownames(empty))
    expect_identical(colnames(empty), colnames(mat))
})

test_that("randomized compiled calls preserve list-path scores and bytes", {
    set.seed(18072028)
    serial <- BiocParallel::SerialParam()
    for (iteration in seq_len(12L)) {
        feature_count <- sample(4:18, 1L)
        sample_count <- sample(1:9, 1L)
        set_count <- sample(1:8, 1L)
        features <- sprintf("feature_%02d", seq_len(feature_count))
        pool <- c(features, sprintf("absent_%02d", seq_len(4L)))
        gene_sets <- lapply(seq_len(set_count), function(index) {
            sample(pool, sample(0:12, 1L), replace = TRUE)
        })
        names(gene_sets) <- sprintf("set_%02d", seq_len(set_count))
        mat <- matrix(
            stats::rexp(feature_count * sample_count),
            nrow = feature_count,
            dimnames = list(
                features,
                sprintf("sample_%02d", seq_len(sample_count))
            )
        )
        mat[sample(length(mat), floor(length(mat) / 7L))] <- 0
        mat[sample(length(mat), floor(length(mat) / 11L))] <- NA_real_

        compiled <- genefunnel:::.compile_gene_sets(gene_sets, features)
        expected <- suppressWarnings(genefunnel(mat, gene_sets, BPPARAM = serial))
        dense <- suppressWarnings(
            genefunnel:::.genefunnel_compiled(mat, compiled, serial)
        )
        sparse <- suppressWarnings(genefunnel:::.genefunnel_compiled(
            Matrix::Matrix(mat, sparse = TRUE),
            compiled,
            serial
        ))
        restored <- genefunnel:::.unserialize_compiled_catalogue(
            genefunnel:::.serialize_compiled_catalogue(compiled)
        )

        expect_identical(dense, expected, info = iteration)
        expect_identical(sparse, expected, info = iteration)
        expect_identical(restored, compiled, info = iteration)
    }
})

test_that("compiled scoring rejects every stale feature universe", {
    mat <- catalogue_matrix()
    compiled <- genefunnel:::.compile_gene_sets(catalogue_fixture(), rownames(mat))
    serial <- BiocParallel::SerialParam()

    reordered <- mat[c(2L, 1L, 3L, 4L), , drop = FALSE]
    renamed <- mat
    rownames(renamed)[[4L]] <- "changed"

    expect_error(
        genefunnel:::.genefunnel_compiled(reordered, compiled, serial),
        "feature identities and order"
    )
    expect_error(
        genefunnel:::.genefunnel_compiled(renamed, compiled, serial),
        "feature identities and order"
    )
    expect_error(
        genefunnel:::.genefunnel_compiled(unname(mat), compiled, serial),
        "must have row names"
    )

    invalid <- mat
    invalid[1L, 1L] <- Inf
    expect_error(
        genefunnel:::.genefunnel_compiled(invalid, compiled, "not a backend"),
        "infinite values"
    )
    expect_error(
        genefunnel:::.genefunnel_compiled(mat, compiled, "not a backend"),
        "BiocParallelParam"
    )
})

test_that("all compiled state mutations fail before scoring", {
    compiled <- genefunnel:::.compile_gene_sets(
        catalogue_fixture(),
        rownames(catalogue_matrix())
    )
    mutations <- list(
        schema = function(x) { x$schema_version <- 2L; x },
        formula = function(x) { x$formula_version <- "changed"; x },
        features = function(x) { x$features[[1L]] <- "changed"; x },
        members = function(x) { x$members$full[[1L]] <- "changed"; x },
        duplicate_count = function(x) { x$duplicate_member_count[[1L]] <- 0L; x },
        member_index = function(x) { x$member_index$full[[1L]] <- 2L; x },
        indices = function(x) { x$indices$full <- rev(x$indices$full); x },
        coverage = function(x) { x$coverage$matched_size[[1L]] <- 1L; x },
        adjacency = function(x) { x$row_adjacency$set_indices[[1L]] <- 3L; x },
        fingerprint = function(x) { x$fingerprints[["content"]] <- strrep("0", 64L); x }
    )

    for (label in names(mutations)) {
        candidate <- mutations[[label]](compiled)
        expect_error(
            genefunnel:::.validate_compiled_catalogue(candidate),
            "malformed",
            info = label
        )
    }

    candidate <- unclass(compiled)
    expect_error(
        genefunnel:::.validate_compiled_catalogue(candidate),
        "genefunnel_catalogue"
    )
})

test_that("GFCAT bytes round-trip deterministically with encoding marks", {
    latin <- "fran\xe7ais"
    Encoding(latin) <- "latin1"
    utf8 <- enc2utf8("\u6771\u4eac")
    byte_id <- rawToChar(as.raw(c(0x72, 0x61, 0x77, 0xff)))
    Encoding(byte_id) <- "bytes"
    features <- c(latin, utf8, byte_id)
    gene_sets <- list(encoded = c(latin, utf8, byte_id, latin), empty = character())
    compiled <- genefunnel:::.compile_gene_sets(gene_sets, features)

    first <- genefunnel:::.serialize_compiled_catalogue(compiled)
    second <- genefunnel:::.serialize_compiled_catalogue(compiled)
    restored <- genefunnel:::.unserialize_compiled_catalogue(first)
    restored_state <- unclass(restored)

    expect_type(first, "raw")
    expect_identical(first, second)
    expect_identical(restored, compiled)
    expect_identical(Encoding(restored_state$features), Encoding(features))
    expect_identical(
        lapply(restored_state$features, charToRaw),
        lapply(features, charToRaw)
    )
    expect_identical(
        restored_state$fingerprints[["content"]],
        unname(tools::sha256sum(bytes = first[seq_len(length(first) - 32L)]))
    )
})

test_that("GFCAT parser rejects corrupt, truncated, trailing, and bounded invalid bytes", {
    compiled <- genefunnel:::.compile_gene_sets(catalogue_fixture(), LETTERS[1:4])
    encoded <- genefunnel:::.serialize_compiled_catalogue(compiled)

    corrupt <- encoded
    corrupt[[10L]] <- as.raw(bitwXor(as.integer(corrupt[[10L]]), 1L))
    expect_error(
        genefunnel:::.unserialize_compiled_catalogue(corrupt),
        "fingerprint"
    )
    expect_error(
        genefunnel:::.unserialize_compiled_catalogue(encoded[-length(encoded)]),
        "fingerprint|truncated"
    )
    expect_error(
        genefunnel:::.unserialize_compiled_catalogue(c(encoded, as.raw(0L))),
        "fingerprint"
    )
    expect_error(
        genefunnel:::.unserialize_compiled_catalogue(structure(encoded, class = "rawish")),
        "unclassed raw vector"
    )

    body <- encoded[seq_len(length(encoded) - 32L)]
    invalid_tag <- body
    invalid_tag[[10L]] <- as.raw(9L)
    invalid_tag <- genefunnel:::.catalogue_append_digest(invalid_tag)
    expect_error(
        genefunnel:::.unserialize_compiled_catalogue(invalid_tag),
        "encoding tag"
    )

    impossible_length <- body
    impossible_length[11:14] <- as.raw(rep(0xff, 4L))
    impossible_length <- genefunnel:::.catalogue_append_digest(impossible_length)
    expect_error(
        genefunnel:::.unserialize_compiled_catalogue(impossible_length),
        "length|integer range"
    )

    trailing_body <- genefunnel:::.catalogue_append_digest(c(body, as.raw(0L)))
    expect_error(
        genefunnel:::.unserialize_compiled_catalogue(trailing_body),
        "trailing"
    )

    tampered_state <- unclass(compiled)
    tampered_state$member_index$partial[[2L]] <- 1L
    self_consistent_digest <- genefunnel:::.catalogue_append_digest(
        genefunnel:::.catalogue_content_body(tampered_state)
    )
    expect_error(
        genefunnel:::.unserialize_compiled_catalogue(self_consistent_digest),
        "mappings or adjacency"
    )
})

test_that("custom serialization and scoring survive fresh and reused SOCK workers", {
    mat <- catalogue_matrix()
    compiled <- genefunnel:::.compile_gene_sets(catalogue_fixture(), rownames(mat))
    encoded <- genefunnel:::.serialize_compiled_catalogue(compiled)
    serial <- BiocParallel::SerialParam()
    expect_warning(
        expected <- genefunnel:::.genefunnel_compiled(mat, compiled, serial),
        "Omitting 2 gene sets"
    )

    fresh <- BiocParallel::SnowParam(workers = 2L, type = "SOCK")
    fresh_bytes <- BiocParallel::bplapply(
        list(encoded, encoded),
        function(payload) {
            object <- genefunnel:::.unserialize_compiled_catalogue(payload)
            genefunnel:::.serialize_compiled_catalogue(object)
        },
        BPPARAM = fresh
    )
    expect_identical(fresh_bytes, list(encoded, encoded))

    reused <- BiocParallel::bpstart(BiocParallel::SnowParam(
        workers = 2L,
        type = "SOCK"
    ))
    on.exit(BiocParallel::bpstop(reused), add = TRUE)
    expect_warning(
        first <- genefunnel:::.genefunnel_compiled(mat, compiled, reused),
        "Omitting 2 gene sets"
    )
    expect_warning(
        second <- genefunnel:::.genefunnel_compiled(
            Matrix::Matrix(mat, sparse = TRUE),
            compiled,
            reused
        ),
        "Omitting 2 gene sets"
    )
    expect_identical(first, expected)
    expect_identical(second, expected)
    expect_true(BiocParallel::bpisup(reused))
})
