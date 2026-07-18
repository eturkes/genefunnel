# Assisted-by: OpenAI Codex.

AGGREGATION_PROTOCOL_VERSION <- "B-1.0.0"
AGGREGATION_PROTOCOL_SHA256 <-
    "2fb9e937c6c6ac49952e2a7351ec25821c787a786a58d21f4c2183dd246d5810"
AGGREGATION_DATA_SHA256 <-
    "52bf82d58aa4dd2843f97b47aaec23f496262c5968564633b7899dc1544ee70f"

aggregation_read_registry <- function(path) {
    if (!identical(
        unname(tools::sha256sum(path)),
        AGGREGATION_PROTOCOL_SHA256
    )) {
        stop("Aggregation protocol bytes changed; version it.", call. = FALSE)
    }
    registry <- utils::read.delim(
        path,
        stringsAsFactors = FALSE,
        check.names = FALSE,
        colClasses = "character",
        na.strings = character()
    )
    required <- c("protocol_version", "section", "key", "value")
    if (!identical(names(registry), required) ||
        !all(registry$protocol_version == AGGREGATION_PROTOCOL_VERSION) ||
        any(!nzchar(registry$section)) || any(!nzchar(registry$key)) ||
        any(!nzchar(registry$value))) {
        stop("Aggregation protocol schema or values are invalid.", call. = FALSE)
    }
    keys <- paste(registry$section, registry$key, sep = "::")
    if (anyDuplicated(keys)) {
        stop("Aggregation protocol keys must be unique.", call. = FALSE)
    }
    registry
}

aggregation_registry_value <- function(registry, section, key) {
    selected <- registry$section == section & registry$key == key
    if (sum(selected) != 1L) {
        stop(
            "Aggregation protocol key is absent or duplicated: ",
            section,
            "::",
            key,
            call. = FALSE
        )
    }
    registry$value[selected]
}

aggregation_registry_vector <- function(registry, section, key) {
    strsplit(
        aggregation_registry_value(registry, section, key),
        ",",
        fixed = TRUE
    )[[1L]]
}

aggregation_read_data_manifest <- function(path) {
    if (!identical(
        unname(tools::sha256sum(path)),
        AGGREGATION_DATA_SHA256
    )) {
        stop("Aggregation data manifest bytes changed; version it.", call. = FALSE)
    }
    manifest <- utils::read.delim(
        path,
        stringsAsFactors = FALSE,
        check.names = FALSE,
        colClasses = "character",
        na.strings = character()
    )
    required <- c(
        "protocol_version", "data_id", "role", "url", "bytes", "sha256"
    )
    bytes <- suppressWarnings(as.numeric(manifest$bytes))
    if (!identical(names(manifest), required) ||
        !all(manifest$protocol_version == AGGREGATION_PROTOCOL_VERSION) ||
        anyDuplicated(manifest$data_id) || anyDuplicated(manifest$url) ||
        any(!nzchar(manifest$data_id)) || any(!nzchar(manifest$role)) ||
        any(!startsWith(manifest$url, "https://")) ||
        any(!is.finite(bytes) | bytes < 1 | bytes != floor(bytes)) ||
        any(!grepl("^[0-9a-f]{64}$", manifest$sha256))) {
        stop("Aggregation data manifest is invalid.", call. = FALSE)
    }
    manifest$bytes <- bytes
    manifest
}

aggregation_verify_data_files <- function(manifest, files) {
    if (!is.character(files) || is.null(names(files)) ||
        anyDuplicated(names(files)) || length(files) != nrow(manifest) ||
        !setequal(names(files), manifest$data_id)) {
        stop("Aggregation data files must be named by every data ID.", call. = FALSE)
    }
    files <- unname(files[manifest$data_id])
    if (any(!file.exists(files))) {
        stop("One or more aggregation data files are absent.", call. = FALSE)
    }
    bytes <- as.numeric(file.info(files)$size)
    hashes <- unname(tools::sha256sum(files))
    if (!identical(bytes, manifest$bytes) ||
        !identical(hashes, manifest$sha256)) {
        stop("Aggregation data bytes disagree with the manifest.", call. = FALSE)
    }
    invisible(data.frame(
        data_id = manifest$data_id,
        path = normalizePath(files),
        bytes = bytes,
        sha256 = hashes,
        stringsAsFactors = FALSE
    ))
}

aggregation_factorial_signs <- function() {
    signs <- expand.grid(
        A = c(-1L, 1L),
        B = c(-1L, 1L),
        C = c(-1L, 1L),
        D = c(-1L, 1L),
        E = c(-1L, 1L),
        KEEP.OUT.ATTRS = FALSE
    )
    signs$F <- signs$A * signs$B * signs$C
    signs$G <- signs$B * signs$C * signs$D

    for (degree in 1:3) {
        products <- utils::combn(names(signs), degree, simplify = FALSE)
        balanced <- vapply(products, function(fields) {
            sum(Reduce(`*`, signs[fields])) == 0L
        }, logical(1))
        if (!all(balanced)) {
            stop("Aggregation fractional factorial lost resolution IV.", call. = FALSE)
        }
    }
    signs
}

aggregation_map_sign <- function(sign, values) {
    if (length(values) != 2L) {
        stop("A binary aggregation factor needs two values.", call. = FALSE)
    }
    values[ifelse(sign < 0L, 1L, 2L)]
}

aggregation_synthetic_design <- function(registry) {
    integer_values <- function(key) {
        value <- suppressWarnings(as.integer(aggregation_registry_vector(
            registry,
            "synthetic",
            key
        )))
        if (anyNA(value)) {
            stop("Aggregation integer factor is invalid: ", key, call. = FALSE)
        }
        value
    }
    numeric_values <- function(key) {
        value <- suppressWarnings(as.numeric(aggregation_registry_vector(
            registry,
            "synthetic",
            key
        )))
        if (anyNA(value) || any(!is.finite(value))) {
            stop("Aggregation numeric factor is invalid: ", key, call. = FALSE)
        }
        value
    }

    core <- expand.grid(
        archetype = aggregation_registry_vector(
            registry,
            "synthetic",
            "archetype"
        ),
        member_count = integer_values("member_count"),
        unit_count = integer_values("unit_count"),
        weight_profile = aggregation_registry_vector(
            registry,
            "synthetic",
            "weight_profile"
        ),
        library_depth = integer_values("library_depth"),
        dropout_fraction = numeric_values("dropout_fraction"),
        complementarity = numeric_values("complementarity"),
        KEEP.OUT.ATTRS = FALSE,
        stringsAsFactors = FALSE
    )
    signs <- aggregation_factorial_signs()
    factorial <- data.frame(
        baseline_profile = aggregation_map_sign(
            signs$A,
            aggregation_registry_vector(
                registry,
                "synthetic",
                "baseline_profile"
            )
        ),
        dynamic_range = as.numeric(aggregation_map_sign(
            signs$B,
            numeric_values("dynamic_range")
        )),
        log_sd = as.numeric(aggregation_map_sign(
            signs$C,
            numeric_values("log_sd")
        )),
        correlation = as.numeric(aggregation_map_sign(
            signs$D,
            numeric_values("correlation")
        )),
        overlap_fraction = as.numeric(aggregation_map_sign(
            signs$E,
            numeric_values("overlap_fraction")
        )),
        outlier_multiplier = as.numeric(aggregation_map_sign(
            signs$F,
            numeric_values("outlier_multiplier")
        )),
        subunit_count = as.integer(aggregation_map_sign(
            signs$G,
            integer_values("subunit_count")
        )),
        stringsAsFactors = FALSE
    )

    core_rows <- rep(seq_len(nrow(core)), each = nrow(factorial))
    factorial_rows <- rep(seq_len(nrow(factorial)), times = nrow(core))
    latent <- cbind(
        core[core_rows, , drop = FALSE],
        factorial[factorial_rows, , drop = FALSE]
    )
    latent_id <- seq_len(nrow(latent))
    replicates <- aggregation_registry_vector(
        registry,
        "synthetic",
        "measurement_replicates"
    )
    if (!identical(replicates, c("A", "B"))) {
        stop("Aggregation measurement replicates are invalid.", call. = FALSE)
    }
    rows <- rep(latent_id, each = length(replicates))
    design <- latent[rows, , drop = FALSE]
    design$measurement_replicate <- rep(replicates, times = nrow(latent))
    base_seed <- as.integer(aggregation_registry_value(
        registry,
        "general",
        "simulation_seed"
    ))
    design$seed <- base_seed + rows + ifelse(
        design$measurement_replicate == "B",
        1000000L,
        0L
    )
    design$scenario_id <- sprintf(
        "B-S%06d-%s",
        rows,
        design$measurement_replicate
    )
    design <- design[c(
        "scenario_id", "seed", "measurement_replicate",
        setdiff(names(design), c(
            "scenario_id", "seed", "measurement_replicate"
        ))
    )]
    rownames(design) <- NULL

    expected <- as.integer(aggregation_registry_value(
        registry,
        "synthetic",
        "scenario_count"
    ))
    if (!identical(nrow(core), 1944L) ||
        !identical(nrow(factorial), 32L) ||
        !identical(nrow(latent), 62208L) ||
        !identical(nrow(design), expected) || anyDuplicated(design$scenario_id) ||
        any(design$seed > .Machine$integer.max)) {
        stop("Aggregation synthetic design dimensions are invalid.", call. = FALSE)
    }
    design
}

aggregation_weights <- function(unit_count, profile) {
    key <- paste(unit_count, profile, sep = "::")
    weights <- switch(
        key,
        `2::equal` = c(0.5, 0.5),
        `2::moderate` = c(0.25, 0.75),
        `2::dominant` = c(0.1, 0.9),
        `4::equal` = rep(0.25, 4L),
        `4::moderate` = c(0.1, 0.2, 0.3, 0.4),
        `4::dominant` = c(0.05, 0.05, 0.1, 0.8),
        stop("Unknown aggregation weight profile.", call. = FALSE)
    )
    stopifnot(abs(sum(weights) - 1) <= 4 * .Machine$double.eps)
    weights
}

aggregation_validate_protocol <- function(root = ".") {
    registry <- aggregation_read_registry(file.path(
        root,
        "benchmark",
        "aggregation-protocol.tsv"
    ))
    manifest <- aggregation_read_data_manifest(file.path(
        root,
        "benchmark",
        "aggregation-data.tsv"
    ))
    design <- aggregation_synthetic_design(registry)
    for (unit_count in c(2L, 4L)) {
        for (profile in c("equal", "moderate", "dominant")) {
            aggregation_weights(unit_count, profile)
        }
    }
    list(
        registry_rows = nrow(registry),
        data_files = nrow(manifest),
        synthetic_rows = nrow(design),
        latent_scenarios = nrow(design) / 2L
    )
}
