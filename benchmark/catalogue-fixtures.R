# Assisted-by: OpenAI Codex.

CATALOGUE_PROTOCOL_VERSION <- "C-1.0.2"
CATALOGUE_PROTOCOL_SHA256 <-
    "0bbab0adaf4264eed31c4ec81673dd5d0a3ab2ecc1e159695681c2c603920b21"

catalogue_read_protocol <- function(path, mode) {
    if (!identical(
        unname(tools::sha256sum(path)),
        CATALOGUE_PROTOCOL_SHA256
    )) {
        stop("Catalogue protocol bytes changed; version the protocol.", call. = FALSE)
    }
    protocol <- utils::read.delim(
        path,
        stringsAsFactors = FALSE,
        check.names = FALSE,
        na.strings = "NA"
    )
    required <- c(
        "protocol_version", "scenario_id", "storage", "overlap", "features",
        "samples", "sets", "canonical_set_size", "matrix_seed", "set_seed",
        "stored_density", "dense_zero_fraction", "dense_missing_fraction",
        "sparse_stored_missing_fraction", "unmatched_stride",
        "duplicate_stride", "batches", "workers", "repeats", "order_code"
    )
    if (!identical(names(protocol), required) ||
        !all(protocol$protocol_version == CATALOGUE_PROTOCOL_VERSION) ||
        anyDuplicated(protocol$scenario_id)) {
        stop("Catalogue protocol schema or identity is invalid.", call. = FALSE)
    }
    integer_fields <- c(
        "features", "samples", "sets", "canonical_set_size", "matrix_seed",
        "set_seed", "unmatched_stride", "duplicate_stride", "batches",
        "workers", "repeats"
    )
    valid_integer <- vapply(protocol[integer_fields], function(value) {
        all(is.finite(value) & value >= 1 & value == floor(value))
    }, logical(1))
    if (!all(valid_integer) ||
        any(protocol$features < protocol$sets * protocol$canonical_set_size) ||
        any(protocol$unmatched_stride > protocol$sets) ||
        any(protocol$duplicate_stride > protocol$sets) ||
        any(protocol$batches != 5L) || any(protocol$workers != 1L) ||
        any(protocol$repeats != 20L)) {
        stop("Catalogue protocol dimensions are invalid.", call. = FALSE)
    }
    valid_storage <- protocol$storage %in% c("dense", "sparse")
    valid_overlap <- protocol$overlap %in% c("low", "high")
    fraction_fields <- c(
        "stored_density", "dense_zero_fraction", "dense_missing_fraction",
        "sparse_stored_missing_fraction"
    )
    valid_fractions <- vapply(protocol[fraction_fields], function(value) {
        all(is.na(value) | (is.finite(value) & value >= 0 & value <= 1))
    }, logical(1))
    dense <- protocol$storage == "dense"
    valid_missingness <-
        all(is.na(protocol$stored_density[dense])) &&
        all(is.na(protocol$sparse_stored_missing_fraction[dense])) &&
        all(!is.na(protocol$dense_zero_fraction[dense])) &&
        all(!is.na(protocol$dense_missing_fraction[dense])) &&
        all(!is.na(protocol$stored_density[!dense])) &&
        all(!is.na(protocol$sparse_stored_missing_fraction[!dense])) &&
        all(is.na(protocol$dense_zero_fraction[!dense])) &&
        all(is.na(protocol$dense_missing_fraction[!dense]))
    if (!all(valid_storage) || !all(valid_overlap) ||
        !all(valid_fractions) || !valid_missingness) {
        stop("Catalogue protocol storage/fraction fields are invalid.", call. = FALSE)
    }
    valid_order <- nchar(protocol$order_code) == 20L &
        grepl("^[LC]+$", protocol$order_code) &
        vapply(strsplit(protocol$order_code, ""), function(order) {
            sum(order == "L") == 10L && sum(order == "C") == 10L
        }, logical(1))
    if (!all(valid_order)) {
        stop("Catalogue protocol order schedule is invalid.", call. = FALSE)
    }
    if (mode == "smoke") {
        dense <- protocol$storage == "dense"
        protocol$features <- 500L
        protocol$samples <- ifelse(dense, 8L, 16L)
        protocol$sets <- 40L
        protocol$canonical_set_size <- 8L
        protocol$batches <- 2L
        protocol$repeats <- 4L
        protocol$order_code <- "LCLC"
    }
    protocol
}

catalogue_gene_sets <- function(scenario, features) {
    sets <- benchmark_gene_sets(
        features,
        scenario$sets,
        scenario$canonical_set_size,
        scenario$overlap,
        scenario$set_seed
    )
    unmatched <- seq.int(
        scenario$unmatched_stride,
        scenario$sets,
        by = scenario$unmatched_stride
    )
    for (index in unmatched) {
        sets[[index]][[scenario$canonical_set_size]] <- sprintf(
            "absent_%06d",
            index
        )
    }
    duplicated <- seq.int(
        scenario$duplicate_stride,
        scenario$sets,
        by = scenario$duplicate_stride
    )
    for (index in duplicated) {
        sets[[index]] <- c(sets[[index]], sets[[index]][[1L]])
    }
    sets
}

catalogue_matrix_batch <- function(scenario, batch) {
    seed <- scenario$matrix_seed + batch - 1L
    switch(
        scenario$storage,
        dense = benchmark_dense_matrix(
            scenario$features,
            scenario$samples,
            scenario$dense_zero_fraction,
            scenario$dense_missing_fraction,
            seed
        ),
        sparse = benchmark_sparse_matrix(
            scenario$features,
            scenario$samples,
            scenario$stored_density,
            scenario$sparse_stored_missing_fraction,
            seed
        ),
        stop("Unknown catalogue storage kind.", call. = FALSE)
    )
}

catalogue_fixture <- function(scenario) {
    batches <- lapply(seq_len(scenario$batches), function(batch) {
        catalogue_matrix_batch(scenario, batch)
    })
    features <- rownames(batches[[1L]])
    gene_sets <- catalogue_gene_sets(scenario, features)
    canonical_count <- sum(lengths(lapply(gene_sets, unique)))
    expected_count <- scenario$sets * scenario$canonical_set_size
    if (!identical(as.double(canonical_count), as.double(expected_count))) {
        stop("Catalogue fixture canonical membership count is invalid.", call. = FALSE)
    }
    if (!all(vapply(batches, function(mat) {
        identical(rownames(mat), features)
    }, logical(1)))) {
        stop("Catalogue fixture feature universes disagree.", call. = FALSE)
    }
    list(batches = batches, gene_sets = gene_sets)
}

catalogue_score_digest <- function(scores) {
    bytes <- serialize(scores, NULL, ascii = FALSE, xdr = TRUE, version = 3L)
    unname(tools::sha256sum(bytes = bytes))
}
