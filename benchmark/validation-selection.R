# Assisted-by: OpenAI Codex.

VALIDATION_SELECTION_SCHEMA <- "F-S-1.0.0"
VALIDATION_BYTE_SCHEMA <- "F-B-1.0.0"
VALIDATION_ACCESS_SCHEMA <- "F-A-1.0.0"
VALIDATION_SELECTION_VALIDATOR_PATH <- "benchmark/validation-selection.R"
VALIDATION_RUNTIME_PATHS <- c(
    "benchmark/validation-protocol.R",
    "benchmark/validation-protocol.tsv",
    "benchmark/validation-amendment.R",
    "benchmark/validation-amendment.tsv",
    "benchmark/validation-selection.R",
    "benchmark/validation-reactome97.tsv",
    "benchmark/generate-validation-reactome97.R"
)
VALIDATION_REACTOME97_SHA256 <-
    "915fce5e43b0e0fe57baf1a0b6cf646f5be8c4e25d325612c5ad6263ebc69202"
VALIDATION_REACTOME97_GENERATOR_SHA256 <-
    "b53399f9fe942fdda84d14eaffadfd547dc8bee0329f9db7e3f5d164b8b03b77"
VALIDATION_SELECTION_FILES <- c(
    "sources.tsv", "tasks.tsv", "task-design.tsv", "task-resampling.tsv",
    "task-filters.tsv", "assay-inputs.tsv", "objects-planned.tsv"
)
VALIDATION_SELECTION_FIELDS <- list(
    sources = c(
        "schema_version", "protocol_version", "protocol_sha256", "source_id",
        "source_order", "assay", "access_role", "source_record_id",
        "metadata_locator", "license_id", "license_locator",
        "retrieval_permission", "redistribution_permission",
        "independence_group"
    ),
    tasks = c(
        "schema_version", "protocol_version", "protocol_sha256", "task_id",
        "task_order", "source_id", "task_kind", "contrast_id",
        "population_id", "reference_condition_id",
        "intervention_condition_id", "exposure_id", "collection_time_id",
        "unit_frame_id", "dependence_group_id", "target_id",
        "expected_direction", "citation_id", "citation_locator",
        "null_assignment_mode", "label_seed", "permutation_seed"
    ),
    design = c(
        "schema_version", "protocol_version", "protocol_sha256", "task_id",
        "field_order", "field_role", "metadata_field", "field_type"
    ),
    resampling = c(
        "schema_version", "protocol_version", "protocol_sha256", "task_id",
        "pool_order", "resample_pool_id", "pool_role", "arm"
    ),
    filters = c(
        "schema_version", "protocol_version", "protocol_sha256", "task_id",
        "filter_order", "arm", "filter_role", "metadata_scope",
        "metadata_field", "metadata_value"
    ),
    inputs = c(
        "schema_version", "protocol_version", "protocol_sha256", "source_id",
        "values_object_id", "values_member", "values_container",
        "values_adapter",
        "matrix_orientation", "observation_id_field", "feature_id_field",
        "quantity_field", "quantity_measure", "quantity_selector",
        "observation_selection_rule",
        "sample_metadata_object_id", "sample_metadata_member",
        "sample_metadata_container", "sample_metadata_adapter",
        "sample_id_field", "technical_repeat_id_field",
        "technical_repeat_aggregation",
        "cell_metadata_object_id", "cell_metadata_member",
        "cell_metadata_container", "cell_metadata_adapter", "cell_id_field",
        "cell_sample_id_field", "quantity_level",
        "quantity_scale", "normalization_state", "imputation_state",
        "zero_semantics", "missing_semantics", "missing_tokens",
        "assignment_or_inference", "species_rule", "group_rule",
        "duplicate_accession_operator", "mapping_object_id", "mapping_member",
        "mapping_container", "mapping_adapter", "mapping_source_field",
        "mapping_target_field", "mapping_rule", "evidence_locator"
    ),
    objects = c(
        "schema_version", "protocol_version", "protocol_sha256", "object_id",
        "object_order", "source_id", "object_kind", "access_class",
        "request_locator"
    )
)
VALIDATION_BYTE_FIELDS <- c(
    "schema_version", "protocol_version", "protocol_sha256",
    "selection_commit", "selection_sha256", "object_id", "status",
    "resolved_locator", "terminal_http_status", "retrieved_utc", "bytes",
    "sha256", "failure_code"
)
VALIDATION_ACCESS_FIELDS <- c(
    "schema_version", "protocol_version", "protocol_sha256", "event_order",
    "event_type", "repository_head", "worktree_state", "selection_sha256",
    "byte_manifest_sha256", "implementation_commit",
    "execution_closure_sha256", "execution_attempt_id", "object_id",
    "event_utc", "bytes", "sha256"
)
VALIDATION_REACTOME97_FIELDS <- c("assay", "target_id", "declared_size")

validation_selection_require_amendment <- function(root) {
    if (!exists("validation_effective_protocol_validate", mode = "function",
            inherits = TRUE)) {
        stop(
            "Source benchmark/validation-amendment.R before the selection validator.",
            call. = FALSE
        )
    }
    effective <- validation_effective_protocol_validate(root)
    protocol <- validation_effective_protocol_read(root)
    registered_fields <- list(
        sources = "sources_fields", tasks = "tasks_fields",
        design = "design_fields", resampling = "resampling_fields",
        filters = "filters_fields", inputs = "inputs_fields",
        objects = "objects_fields"
    )
    fields_valid <- all(vapply(names(registered_fields), function(name) {
        value <- validation_effective_value(
            protocol, "selection", registered_fields[[name]]
        )
        identical(strsplit(value, ",", fixed = TRUE)[[1L]],
            VALIDATION_SELECTION_FIELDS[[name]])
    }, logical(1L))) && identical(strsplit(validation_effective_value(
        protocol, "bytes", "fields"
    ), ",", fixed = TRUE)[[1L]], VALIDATION_BYTE_FIELDS) &&
        identical(strsplit(validation_effective_value(
            protocol, "access", "fields"
        ), ",", fixed = TRUE)[[1L]], VALIDATION_ACCESS_FIELDS)
    validator_path <- validation_effective_value(
        protocol, "selection", "validator_path"
    )
    validator_sha256 <- validation_effective_value(
        protocol, "selection", "validator_sha256"
    )
    runtime_paths <- strsplit(validation_effective_value(
        protocol, "selection", "runtime_paths"
    ), ",", fixed = TRUE)[[1L]]
    live_validator <- file.path(root, validator_path)
    validator_valid <- identical(
        validator_path, VALIDATION_SELECTION_VALIDATOR_PATH
    ) && grepl("^[0-9a-f]{64}$", validator_sha256) &&
        validation_contract_is_type(live_validator, "regular file") &&
        !nzchar(Sys.readlink(live_validator)) &&
        identical(validation_contract_sha256(live_validator), validator_sha256)
    if (!identical(effective$protocol_version, VALIDATION_AMENDMENT_VERSION) ||
        !identical(effective$protocol_sha256, VALIDATION_AMENDMENT_SHA256) ||
        !identical(effective$selection_schema, VALIDATION_SELECTION_SCHEMA) ||
        !identical(effective$byte_schema, VALIDATION_BYTE_SCHEMA) ||
        !identical(effective$access_schema, VALIDATION_ACCESS_SCHEMA) ||
        !fields_valid || !validator_valid ||
        !identical(runtime_paths, VALIDATION_RUNTIME_PATHS)) {
        stop("Selection validator protocol identity is invalid.", call. = FALSE)
    }
    effective
}

validation_reactome97_read <- function(root = ".") {
    path <- file.path(root, "benchmark", "validation-reactome97.tsv")
    generator <- file.path(root, "benchmark", "generate-validation-reactome97.R")
    if (!identical(validation_contract_sha256(path),
            VALIDATION_REACTOME97_SHA256) ||
        !identical(validation_contract_sha256(generator),
            VALIDATION_REACTOME97_GENERATOR_SHA256)) {
        stop("Reactome 97 target roster identity is invalid.", call. = FALSE)
    }
    roster <- validation_contract_read_tsv(path, VALIDATION_REACTOME97_FIELDS)
    size <- validation_selection_integer(
        roster$declared_size, label = "declared Reactome size"
    )
    assays <- c("bulk_rnaseq", "pseudobulk_rnaseq", "bulk_proteomics")
    assay_order <- match(roster$assay, assays)
    ordered <- order(assay_order, roster$target_id, method = "radix")
    key <- paste(roster$assay, roster$target_id, sep = "\037")
    counts <- table(factor(roster$assay, levels = assays))
    rna <- roster$assay == "bulk_rnaseq"
    pseudo <- roster$assay == "pseudobulk_rnaseq"
    valid <- !anyNA(assay_order) && all(grepl(
        "^R-HSA-[0-9]+$", roster$target_id
    )) && all(size >= 10 & size <= 200) && !anyDuplicated(key) &&
        identical(ordered, seq_len(nrow(roster))) &&
        identical(as.integer(counts), c(1686L, 1686L, 1690L)) &&
        identical(roster$target_id[rna], roster$target_id[pseudo]) &&
        identical(roster$declared_size[rna], roster$declared_size[pseudo])
    if (!valid) stop("Reactome 97 target roster is invalid.", call. = FALSE)
    roster
}

validation_selection_integer <- function(value, zero = FALSE, label = "integer") {
    pattern <- if (zero) "^(0|[1-9][0-9]*)$" else "^[1-9][0-9]*$"
    valid <- grepl(pattern, value) & nchar(value, type = "bytes") <= 15L
    number <- suppressWarnings(as.numeric(value))
    valid <- valid & is.finite(number) & number <= 9007199254740991
    if (!all(valid)) stop("Noncanonical ", label, ".", call. = FALSE)
    number
}

validation_selection_normalize_percent <- function(value) {
    vapply(value, function(item) {
        characters <- strsplit(item, "", fixed = TRUE)[[1L]]
        output <- character()
        i <- 1L
        while (i <= length(characters)) {
            if (characters[[i]] != "%") {
                output <- c(output, characters[[i]])
                i <- i + 1L
                next
            }
            if (i + 2L > length(characters)) {
                stop("Invalid URL percent encoding.", call. = FALSE)
            }
            digits <- paste0(characters[[i + 1L]], characters[[i + 2L]])
            if (!grepl("^[0-9A-Fa-f]{2}$", digits)) {
                stop("Invalid URL percent encoding.", call. = FALSE)
            }
            code <- strtoi(digits, base = 16L)
            decoded <- intToUtf8(code)
            output <- c(output, if (grepl("^[A-Za-z0-9._~-]$", decoded)) {
                decoded
            } else {
                paste0("%", toupper(digits))
            })
            i <- i + 3L
        }
        paste0(output, collapse = "")
    }, character(1L), USE.NAMES = FALSE)
}

validation_selection_url <- function(value, label, opaque_resource = FALSE) {
    remainder <- sub("^https://", "", value)
    authority <- sub("[/?#].*$", "", remainder)
    labels <- strsplit(authority, ".", fixed = TRUE)
    authority_valid <- vapply(seq_along(labels), function(i) {
        parts <- labels[[i]]
        length(parts) > 0L &&
            identical(paste(parts, collapse = "."), authority[[i]]) &&
            nchar(authority[[i]], type = "bytes") <= 253L &&
            !grepl("^[0-9]+$", parts[[length(parts)]]) &&
            all(grepl(
                "^[a-z0-9]([a-z0-9-]{0,61}[a-z0-9])?$", parts
            ))
    }, logical(1L))
    percent_valid <- !grepl("%", value, fixed = TRUE)
    suffix <- substring(remainder, nchar(authority) + 1L)
    path <- sub("[?#].*$", "", suffix)
    normalized_path <- validation_selection_normalize_percent(path)
    suffix_valid <- grepl(
        "^[A-Za-z0-9._~!$&'()*+,;=:@/?#-]*$", suffix
    ) & lengths(regmatches(suffix, gregexpr("#", suffix, fixed = TRUE))) <= 1L
    path_valid <- (path == "" | startsWith(path, "/")) &
        !grepl("//", normalized_path, fixed = TRUE) &
        !grepl("(^|/)[.]{1,2}(/|$)", normalized_path)
    valid <- grepl("^https://[^[:space:]]+$", value) &
        nzchar(authority) & authority_valid & authority == tolower(authority) &
        percent_valid & path_valid & suffix_valid & !grepl("\\\\", value)
    if (opaque_resource) {
        valid <- valid & startsWith(path, "/") & nchar(path) > 1L &
            !endsWith(path, "/") & !grepl("//", path, fixed = TRUE) &
            !grepl("(^|/)[.]{1,2}(/|$)", path) &
            !grepl("#", value, fixed = TRUE) & !grepl("[?]$", value) &
            grepl("^/[A-Za-z0-9._~!$&'()*+,;=:@/?-]+$", suffix)
    }
    if (!all(valid)) {
        stop("Invalid HTTPS ", label, ".", call. = FALSE)
    }
    invisible(value)
}

validation_selection_url_identity <- function(value) {
    vapply(value, function(item) {
        item <- sub("#.*$", "", item)
        remainder <- sub("^https://", "", item)
        authority <- sub("[/?].*$", "", remainder)
        suffix <- substring(remainder, nchar(authority) + 1L)
        if (suffix == "") suffix <- "/"
        if (startsWith(suffix, "?")) suffix <- paste0("/", suffix)
        query_start <- regexpr("?", suffix, fixed = TRUE)[[1L]]
        if (query_start < 0L) {
            path <- suffix
            query <- ""
        } else {
            path <- substring(suffix, 1L, query_start - 1L)
            query <- substring(suffix, query_start)
        }
        if (nchar(path) > 1L && endsWith(path, "/")) {
            path <- substring(path, 1L, nchar(path) - 1L)
        }
        paste0(
            "https://", authority,
            validation_selection_normalize_percent(paste0(path, query))
        )
    }, character(1L), USE.NAMES = FALSE)
}

validation_selection_ids <- function(value, label) {
    if (!all(grepl("^[a-z][a-z0-9_]{0,63}$", value)) ||
        any(value %in% c("none", "not_applicable"))) {
        stop("Invalid ", label, ".", call. = FALSE)
    }
    invisible(value)
}

validation_selection_identity <- function(table) {
    valid <- all(table$schema_version == VALIDATION_SELECTION_SCHEMA) &&
        all(table$protocol_version == VALIDATION_AMENDMENT_VERSION) &&
        all(table$protocol_sha256 == VALIDATION_AMENDMENT_SHA256)
    if (!valid) stop("Selection table identity is invalid.", call. = FALSE)
    invisible(table)
}

validation_selection_bundle_paths <- function(bundle_dir) {
    if (!validation_contract_is_type(bundle_dir, "directory") ||
        nzchar(Sys.readlink(bundle_dir))) {
        stop("Selection bundle is not a regular directory.", call. = FALSE)
    }
    observed <- sort(list.files(
        bundle_dir, all.files = TRUE, no.. = TRUE, recursive = TRUE,
        include.dirs = TRUE
    ), method = "radix")
    expected <- sort(VALIDATION_SELECTION_FILES, method = "radix")
    paths <- file.path(bundle_dir, VALIDATION_SELECTION_FILES)
    if (!identical(observed, expected) ||
        any(nzchar(Sys.readlink(paths))) ||
        !all(vapply(
            paths, validation_contract_is_type, logical(1L),
            expected = "regular file"
        ))) {
        stop("Selection bundle file closure is invalid.", call. = FALSE)
    }
    paths
}

validation_selection_bundle_sha256 <- function(bundle_dir) {
    paths <- validation_selection_bundle_paths(bundle_dir)
    scratch <- tempfile("validation-selection-hash-")
    connection <- file(scratch, open = "wb")
    on.exit({
        if (inherits(connection, "connection")) close(connection)
        unlink(scratch)
    }, add = TRUE)
    writeBin(charToRaw("GENEFUNNEL-F-S-1.0.0"), connection)
    writeBin(as.raw(0L), connection)
    for (i in seq_along(paths)) {
        before <- file.info(paths[[i]])
        size <- before$size[[1L]]
        input <- file(paths[[i]], open = "rb")
        bytes <- readBin(input, what = "raw", n = size)
        close(input)
        after <- file.info(paths[[i]])
        stable_fields <- intersect(
            c("size", "mtime", "ctime", "ino"), names(before)
        )
        if (length(bytes) != size ||
            !identical(before[stable_fields], after[stable_fields]) ||
            !validation_contract_is_type(paths[[i]], "regular file") ||
            nzchar(Sys.readlink(paths[[i]]))) {
            stop("Selection bundle changed while hashing.", call. = FALSE)
        }
        writeBin(charToRaw(VALIDATION_SELECTION_FILES[[i]]), connection)
        writeBin(as.raw(0L), connection)
        writeBin(charToRaw(sprintf("%015d", size)), connection)
        writeBin(as.raw(0L), connection)
        writeBin(bytes, connection)
        writeBin(as.raw(0L), connection)
    }
    close(connection)
    connection <- NULL
    hash <- validation_contract_sha256(scratch)
    unlink(scratch)
    hash
}

validation_selection_group_rows <- function(table, group, order_field) {
    split_rows <- split(seq_len(nrow(table)), group, drop = TRUE)
    all(vapply(split_rows, function(rows) {
        values <- validation_selection_integer(
            table[[order_field]][rows], label = order_field
        )
        all(values == seq_along(rows))
    }, logical(1L)))
}

validation_selection_validate_sources <- function(sources) {
    validation_selection_identity(sources)
    validation_selection_ids(sources$source_id, "source_id")
    validation_selection_ids(sources$independence_group, "independence_group")
    order <- validation_selection_integer(sources$source_order, label = "source_order")
    valid <- all(order == seq_len(nrow(sources))) &&
        !anyDuplicated(sources$source_id) &&
        !anyDuplicated(paste(sources$assay, sources$source_record_id, sep = "::")) &&
        !anyDuplicated(paste(
            sources$assay, sources$independence_group, sep = "::"
        )) &&
        all(sources$assay %in% c(
            "bulk_rnaseq", "pseudobulk_rnaseq", "bulk_proteomics"
        )) &&
        all(sources$access_role %in% c("development", "heldout")) &&
        all(grepl("^[A-Za-z0-9][A-Za-z0-9._:-]*$", sources$source_record_id)) &&
        all(grepl("^[A-Za-z0-9][A-Za-z0-9._:+-]*$", sources$license_id)) &&
        all(!sources$source_record_id %in% c("none", "not_applicable")) &&
        all(!sources$license_id %in% c("none", "not_applicable")) &&
        all(sources$retrieval_permission == "scripted_retrieval_and_analysis") &&
        all(sources$redistribution_permission %in% c(
            "allowed", "forbidden", "unknown"
        ))
    validation_selection_url(sources$metadata_locator, "metadata locator")
    validation_selection_url(sources$license_locator, "license locator")
    license_valid <- all(vapply(
        split(validation_selection_url_identity(sources$license_locator),
            sources$license_id),
        function(x) length(unique(x)) == 1L,
        logical(1L)
    ))
    record_identity <- paste(
        sources$source_record_id,
        validation_selection_url_identity(sources$metadata_locator), sep = "\037"
    )
    record_valid <- all(vapply(split(
        sources$independence_group, record_identity
    ), function(x) length(unique(x)) == 1L, logical(1L)))
    role_group <- sources$independence_group
    role_consistent <- all(vapply(split(sources$access_role, role_group), function(x) {
        length(unique(x)) == 1L
    }, logical(1L)))
    if (!valid || !role_consistent || !license_valid || !record_valid) {
        stop("Selection sources are invalid.", call. = FALSE)
    }
    invisible(sources)
}

validation_selection_validate_tasks <- function(tasks, sources, roster) {
    validation_selection_identity(tasks)
    validation_selection_ids(tasks$task_id, "task_id")
    validation_selection_ids(tasks$contrast_id, "contrast_id")
    validation_selection_ids(tasks$population_id, "population_id")
    validation_selection_ids(tasks$reference_condition_id,
        "reference_condition_id")
    validation_selection_ids(tasks$intervention_condition_id,
        "intervention_condition_id")
    validation_selection_ids(tasks$exposure_id, "exposure_id")
    validation_selection_ids(tasks$collection_time_id, "collection_time_id")
    validation_selection_ids(tasks$unit_frame_id, "unit_frame_id")
    validation_selection_ids(tasks$dependence_group_id, "dependence_group_id")
    order <- validation_selection_integer(tasks$task_order, label = "task_order")
    label_seed <- validation_selection_integer(
        tasks$label_seed, zero = TRUE, label = "label_seed"
    )
    permutation_seed <- validation_selection_integer(
        tasks$permutation_seed, zero = TRUE, label = "permutation_seed"
    )
    source_row <- match(tasks$source_id, sources$source_id)
    positive <- tasks$task_kind == "positive"
    null <- tasks$task_kind == "null"
    paired_null <- null & tasks$null_assignment_mode == "paired_swap"
    unpaired_null <- null & tasks$null_assignment_mode == "unpaired_balance"
    positive_key <- paste(
        sources$assay[source_row[positive]], tasks$target_id[positive],
        sep = "\037"
    )
    roster_key <- paste(roster$assay, roster$target_id, sep = "\037")
    valid <- all(order == seq_len(nrow(tasks))) &&
        !anyDuplicated(tasks$task_id) && !anyNA(source_row) &&
        setequal(unique(tasks$source_id), sources$source_id) &&
        all(positive | null) &&
        all(tasks$reference_condition_id != tasks$intervention_condition_id) &&
        all(grepl("^[A-Za-z0-9][A-Za-z0-9._:+-]*$", tasks$citation_id) |
            tasks$citation_id == "not_applicable") &&
        all(positive == grepl("^R-HSA-[0-9]+$", tasks$target_id)) &&
        !anyNA(match(positive_key, roster_key)) &&
        all(tasks$expected_direction[positive] %in% c("+1", "-1")) &&
        all(!tasks$citation_id[positive] %in% c("none", "not_applicable")) &&
        all(tasks$null_assignment_mode[positive] == "none") &&
        all(label_seed[positive] == 0) && all(permutation_seed[positive] == 0) &&
        all(tasks$target_id[null] == "all_eligible_pathways") &&
        all(tasks$expected_direction[null] == "0") &&
        all(tasks$null_assignment_mode[null] %in% c(
            "unpaired_balance", "paired_swap"
        )) &&
        all(tasks$citation_id[unpaired_null] == "not_applicable") &&
        all(tasks$citation_locator[unpaired_null] == "not_applicable") &&
        all(tasks$intervention_condition_id[unpaired_null] ==
            "artificial_labels") &&
        all(!tasks$citation_id[paired_null] %in% c(
            "none", "not_applicable"
        )) &&
        all(label_seed[unpaired_null] > 0) &&
        !anyDuplicated(label_seed[unpaired_null]) &&
        all(label_seed[paired_null] == 0) &&
        all(permutation_seed[null] > 0) &&
        !anyDuplicated(permutation_seed[null]) &&
        all(label_seed <= .Machine$integer.max) &&
        all(permutation_seed <= .Machine$integer.max) &&
        all(sources$access_role[source_row[null]] == "heldout")
    cited <- positive | paired_null
    if (any(cited)) validation_selection_url(
        tasks$citation_locator[cited], "task citation locator"
    )
    citation_locator <- validation_selection_url_identity(
        tasks$citation_locator[cited]
    )
    citation_forward <- split(citation_locator, tasks$citation_id[cited])
    citation_reverse <- split(tasks$citation_id[cited], citation_locator)
    citation_valid <- all(vapply(citation_forward, function(x) {
        length(unique(x)) == 1L
    }, logical(1L))) && all(vapply(citation_reverse, function(x) {
        length(unique(x)) == 1L
    }, logical(1L)))
    frame_group <- paste(tasks$source_id, tasks$unit_frame_id, sep = "::")
    frame_consistent <- all(vapply(
        split(tasks$dependence_group_id, frame_group),
        function(x) length(unique(x)) == 1L, logical(1L)
    ))
    dependence_forward <- all(vapply(split(
        sources$independence_group[source_row], tasks$dependence_group_id
    ), function(x) length(unique(x)) == 1L, logical(1L)))
    dependence_reverse <- all(vapply(split(
        tasks$dependence_group_id,
        sources$independence_group[source_row]
    ), function(x) length(unique(x)) == 1L, logical(1L)))
    dependence_consistent <- dependence_forward && dependence_reverse
    contrast_key <- paste(tasks$source_id, tasks$contrast_id, sep = "::")
    positive_fields <- c(
        "source_id", "contrast_id", "population_id", "reference_condition_id",
        "intervention_condition_id", "exposure_id", "collection_time_id",
        "unit_frame_id", "dependence_group_id", "task_kind",
        "null_assignment_mode", "label_seed", "permutation_seed"
    )
    contrast_consistent <- all(vapply(
        split(seq_len(nrow(tasks))[positive], contrast_key[positive]),
        function(rows) length(unique(apply(
            tasks[rows, positive_fields, drop = FALSE], 1L, paste,
            collapse = "\037"
        ))) == 1L,
        logical(1L)
    ))
    null_unique <- !anyDuplicated(contrast_key[null])
    task_key <- paste(
        tasks$source_id[positive], tasks$contrast_id[positive],
        tasks$target_id[positive], sep = "::"
    )
    positive_contrasts <- !duplicated(contrast_key) & positive
    contrast_tuple <- paste(
        tasks$source_id[positive_contrasts],
        tasks$population_id[positive_contrasts],
        tasks$reference_condition_id[positive_contrasts],
        tasks$intervention_condition_id[positive_contrasts],
        tasks$exposure_id[positive_contrasts],
        tasks$collection_time_id[positive_contrasts],
        tasks$unit_frame_id[positive_contrasts], sep = "::"
    )
    contrast_names_unique <- !any(contrast_key[positive] %in% contrast_key[null])
    if (!valid || !citation_valid || !frame_consistent ||
        !dependence_consistent ||
        !contrast_consistent || !null_unique ||
        anyDuplicated(task_key) || anyDuplicated(contrast_tuple) ||
        !contrast_names_unique) {
        stop("Selection tasks are invalid.", call. = FALSE)
    }

    assay <- sources$assay[source_row]
    role <- sources$access_role[source_row]
    independence <- sources$independence_group[source_row]
    for (stratum in c("bulk_rnaseq", "pseudobulk_rnaseq", "bulk_proteomics")) {
        stratum_positive <- positive & assay == stratum
        development <- stratum_positive & role == "development"
        heldout <- stratum_positive & role == "heldout"
        heldout_sources <- unique(tasks$source_id[heldout])
        source_contrasts <- vapply(heldout_sources, function(id) {
            length(unique(tasks$contrast_id[heldout & tasks$source_id == id]))
        }, integer(1L))
        stratum_null <- null & assay == stratum
        null_sources <- unique(tasks$source_id[stratum_null])
        null_groups <- unique(independence[stratum_null])
        null_per_group <- vapply(null_groups, function(id) {
            sum(stratum_null & independence == id)
        }, integer(1L))
        cardinality <- any(development) &&
            length(unique(independence[heldout])) >= 3L &&
            length(heldout_sources) >= 3L && all(source_contrasts >= 2L) &&
            sum(heldout) >= 6L &&
            length(unique(tasks$target_id[heldout])) >= 2L &&
            sum(stratum_null) >= 20L &&
            length(unique(independence[stratum_null])) >= 3L &&
            length(null_sources) >= 3L && all(null_per_group >= 4L)
        if (!cardinality) {
            stop("Selection task cardinality is invalid for ", stratum, ".",
                call. = FALSE)
        }
    }
    invisible(tasks)
}

validation_selection_validate_design <- function(design, tasks) {
    validation_selection_identity(design)
    task_row <- match(design$task_id, tasks$task_id)
    valid <- !anyNA(task_row) &&
        all(design$field_role %in% c("unit", "pair", "block", "covariate")) &&
        all(design$field_type %in% c(
            "identifier", "categorical", "numeric"
        )) &&
        all((design$field_role == "unit") ==
            (design$field_type == "identifier")) &&
        all(!design$field_role %in% c("pair", "block") |
            design$field_type == "categorical") &&
        all(design$field_role != "covariate" |
            design$field_type %in% c("categorical", "numeric")) &&
        all(!design$metadata_field %in% c("none", "not_applicable")) &&
        !any(grepl("[*?]", design$metadata_field)) &&
        !anyDuplicated(paste(
            design$task_id, design$field_role, design$metadata_field, sep = "::"
        )) && !anyDuplicated(paste(
            design$task_id, design$metadata_field, sep = "::"
        )) && validation_selection_group_rows(
            design, design$task_id, "field_order"
        )
    rows <- split(seq_len(nrow(design)), design$task_id, drop = TRUE)
    role_counts <- lapply(rows, function(i) table(factor(
        design$field_role[i], levels = c("unit", "pair", "block", "covariate")
    )))
    shape_valid <- all(vapply(role_counts, function(x) {
        x[["unit"]] == 1L && x[["pair"]] <= 1L && x[["block"]] <= 1L
    }, logical(1L))) && setequal(names(rows), tasks$task_id)
    task_pair_count <- vapply(tasks$task_id, function(id) {
        sum(design$task_id == id & design$field_role == "pair")
    }, integer(1L))
    paired_null <- tasks$task_kind == "null" &
        tasks$null_assignment_mode == "paired_swap"
    unpaired_null <- tasks$task_kind == "null" &
        tasks$null_assignment_mode == "unpaired_balance"
    mode_valid <- all(task_pair_count[paired_null] == 1L) &&
        all(task_pair_count[unpaired_null] == 0L)
    unit_rows <- design$field_role == "unit"
    unit_task <- match(design$task_id[unit_rows], tasks$task_id)
    frame_key <- paste(
        tasks$source_id[unit_task], tasks$unit_frame_id[unit_task], sep = "\037"
    )
    frame_valid <- all(vapply(split(
        design$metadata_field[unit_rows], frame_key
    ), function(x) length(unique(x)) == 1L, logical(1L)))
    if (!valid || !shape_valid || !mode_valid || !frame_valid) {
        stop("Selection task design is invalid.", call. = FALSE)
    }
    invisible(design)
}

validation_selection_pool_signature <- function(
    task, arm, design, filters
) {
    design_rows <- design$task_id == task$task_id & design$field_role == "unit"
    design_values <- sort(paste(
        design$field_role[design_rows], design$metadata_field[design_rows],
        sep = "\037"
    ))
    if (task$task_kind == "null" &&
        task$null_assignment_mode == "unpaired_balance") {
        filter_rows <- filters$task_id == task$task_id & filters$arm == "pool"
    } else if (arm == "joint") {
        filter_rows <- filters$task_id == task$task_id
    } else {
        filter_rows <- filters$task_id == task$task_id &
            filters$arm %in% c("both", arm)
    }
    filter_values <- sort(paste(
        filters$metadata_scope[filter_rows], filters$metadata_field[filter_rows],
        filters$metadata_value[filter_rows], sep = "\037"
    ))
    paste(
        task$source_id, task$unit_frame_id,
        paste(design_values, collapse = "\036"),
        paste(filter_values, collapse = "\036"), sep = "\035"
    )
}

validation_selection_validate_resampling <- function(
    resampling, tasks, design, filters
) {
    validation_selection_identity(resampling)
    validation_selection_ids(resampling$resample_pool_id, "resample_pool_id")
    task_row <- match(resampling$task_id, tasks$task_id)
    valid <- !anyNA(task_row) &&
        all(resampling$pool_role %in% c("joint", "separate")) &&
        all(resampling$arm %in% c(
            "joint", "control", "intervention", "pool"
        )) &&
        !anyDuplicated(paste(
            resampling$task_id, resampling$resample_pool_id,
            resampling$arm, sep = "::"
        )) && validation_selection_group_rows(
            resampling, resampling$task_id, "pool_order"
        ) && setequal(unique(resampling$task_id), tasks$task_id)
    pair_count <- vapply(tasks$task_id, function(id) {
        sum(design$task_id == id & design$field_role == "pair")
    }, integer(1L))
    shape_valid <- all(vapply(seq_len(nrow(tasks)), function(i) {
        rows <- resampling$task_id == tasks$task_id[[i]]
        selected <- resampling[rows, , drop = FALSE]
        if (tasks$task_kind[[i]] == "null" &&
            tasks$null_assignment_mode[[i]] == "unpaired_balance") {
            nrow(selected) == 1L && selected$pool_role == "joint" &&
                selected$arm == "pool"
        } else if (pair_count[[i]] == 1L) {
            nrow(selected) == 1L && selected$pool_role == "joint" &&
                selected$arm == "joint"
        } else {
            nrow(selected) == 2L &&
                all(selected$pool_role == "separate") &&
                setequal(selected$arm, c("control", "intervention"))
        }
    }, logical(1L)))
    signatures <- vapply(seq_len(nrow(resampling)), function(i) {
        task <- tasks[task_row[[i]], , drop = FALSE]
        validation_selection_pool_signature(
            task, resampling$arm[[i]], design, filters
        )
    }, character(1L))
    id_to_signature <- split(signatures, resampling$resample_pool_id)
    signature_to_id <- split(resampling$resample_pool_id, signatures)
    identity_valid <- all(vapply(id_to_signature, function(x) {
        length(unique(x)) == 1L
    }, logical(1L))) && all(vapply(signature_to_id, function(x) {
        length(unique(x)) == 1L
    }, logical(1L)))
    id_to_dependence <- split(
        tasks$dependence_group_id[task_row], resampling$resample_pool_id
    )
    dependence_valid <- all(vapply(id_to_dependence, function(x) {
        length(unique(x)) == 1L
    }, logical(1L)))
    if (!valid || !shape_valid || !identity_valid || !dependence_valid) {
        stop("Selection resampling pools are invalid.", call. = FALSE)
    }
    invisible(resampling)
}

validation_selection_validate_filters <- function(filters, tasks, sources) {
    validation_selection_identity(filters)
    task_row <- match(filters$task_id, tasks$task_id)
    source_row <- match(tasks$source_id, sources$source_id)
    valid <- !anyNA(task_row) &&
        all(filters$arm %in% c("both", "control", "intervention", "pool")) &&
        all(filters$filter_role %in% c(
            "population", "condition", "exposure", "collection_time", "subset"
        )) && all(filters$metadata_scope %in% c("sample", "cell")) &&
        all(!filters$metadata_field %in% c("none", "not_applicable")) &&
        !any(grepl("[*?]", filters$metadata_field)) &&
        !any(grepl("[*?]", filters$metadata_value)) &&
        !anyDuplicated(paste(
            filters$task_id, filters$arm, filters$filter_role,
            filters$metadata_scope, filters$metadata_field,
            filters$metadata_value, sep = "::"
        )) && !anyDuplicated(paste(
            filters$task_id, filters$arm, filters$metadata_scope,
            filters$metadata_field, filters$metadata_value, sep = "::"
        )) && validation_selection_group_rows(
            filters, filters$task_id, "filter_order"
        ) && setequal(unique(filters$task_id), tasks$task_id)
    if (!valid) stop("Selection task filters are invalid.", call. = FALSE)

    for (i in seq_len(nrow(tasks))) {
        rows <- filters$task_id == tasks$task_id[[i]]
        selected <- filters[rows, , drop = FALSE]
        effective_arms <- if (tasks$task_kind[[i]] == "null" &&
            tasks$null_assignment_mode[[i]] == "unpaired_balance") {
            "pool"
        } else {
            c("control", "intervention")
        }
        predicates_compatible <- all(vapply(effective_arms, function(arm) {
            active <- selected$arm %in% c("both", arm)
            groups <- split(
                selected$metadata_value[active],
                paste(
                    selected$metadata_scope[active],
                    selected$metadata_field[active], sep = "\037"
                )
            )
            all(vapply(groups, function(x) {
                length(unique(x)) == 1L
            }, logical(1L)))
        }, logical(1L)))
        predicate_signature <- function(table) {
            sort(paste(
                table$filter_role, table$metadata_scope, table$metadata_field,
                table$metadata_value, sep = "\037"
            ))
        }
        if (tasks$task_kind[[i]] == "positive") {
            condition <- selected$filter_role == "condition"
            arm_valid <- all(selected$arm %in% c(
                "both", "control", "intervention"
            )) && sum(condition) == 2L &&
                setequal(selected$arm[condition], c("control", "intervention")) &&
                sum(selected$arm == "control" &
                selected$filter_role == "condition") == 1L &&
                sum(selected$arm == "intervention" &
                    selected$filter_role == "condition") == 1L &&
                all(selected$arm[selected$filter_role != "condition"] == "both")
            control <- selected[selected$arm == "control", c(
                "filter_role", "metadata_scope", "metadata_field",
                "metadata_value"
            ), drop = FALSE]
            intervention <- selected[selected$arm == "intervention", c(
                "filter_role", "metadata_scope", "metadata_field",
                "metadata_value"
            ), drop = FALSE]
            arm_valid <- arm_valid &&
                identical(control$metadata_scope,
                    intervention$metadata_scope) &&
                identical(control$metadata_field,
                    intervention$metadata_field) &&
                !identical(predicate_signature(control),
                    predicate_signature(intervention))
        } else if (tasks$null_assignment_mode[[i]] == "unpaired_balance") {
            arm_valid <- all(selected$arm == "pool") &&
                sum(selected$filter_role == "condition") == 1L
        } else {
            condition <- selected$filter_role == "condition"
            arm_valid <- all(selected$arm %in% c(
                "both", "control", "intervention"
            )) && sum(condition) == 2L &&
                setequal(selected$arm[condition], c("control", "intervention")) &&
                sum(selected$arm == "control" &
                selected$filter_role == "condition") == 1L &&
                sum(selected$arm == "intervention" &
                    selected$filter_role == "condition") == 1L &&
                all(selected$arm[selected$filter_role != "condition"] == "both")
            control <- selected[selected$arm == "control", c(
                "filter_role", "metadata_scope", "metadata_field",
                "metadata_value"
            ), drop = FALSE]
            intervention <- selected[selected$arm == "intervention", c(
                "filter_role", "metadata_scope", "metadata_field",
                "metadata_value"
            ), drop = FALSE]
            arm_valid <- arm_valid &&
                identical(control$metadata_scope,
                    intervention$metadata_scope) &&
                identical(control$metadata_field,
                    intervention$metadata_field) &&
                !identical(predicate_signature(control),
                    predicate_signature(intervention))
        }
        pseudo <- sources$assay[[source_row[[i]]]] == "pseudobulk_rnaseq"
        population_arm <- if (tasks$task_kind[[i]] == "null" &&
            tasks$null_assignment_mode[[i]] == "unpaired_balance") {
            "pool"
        } else {
            "both"
        }
        population_rows <- selected$filter_role == "population"
        population_expected <- population_rows & selected$arm == population_arm
        exposure_rows <- selected$filter_role == "exposure" &
            selected$arm == population_arm
        time_rows <- selected$filter_role == "collection_time" &
            selected$arm == population_arm
        population_valid <- if (pseudo) {
            sum(population_rows) == 1L && sum(population_expected) == 1L &&
                tasks$population_id[[i]] != "all"
        } else {
            sum(population_rows) <= 1L &&
                (tasks$population_id[[i]] == "all") ==
                (sum(population_rows) == 0L) &&
                (sum(population_rows) == 0L || sum(population_expected) == 1L)
        }
        regimen_valid <- sum(exposure_rows) == 1L && sum(time_rows) == 1L
        scope_valid <- all(selected$metadata_scope == ifelse(
            pseudo & population_rows, "cell", "sample"
        ))
        if (!arm_valid || !population_valid || !regimen_valid || !scope_valid ||
            !predicates_compatible) {
            stop("Selection filters are invalid for task ", tasks$task_id[[i]],
                ".", call. = FALSE)
        }
    }
    population_rows <- filters$filter_role == "population"
    if (any(population_rows)) {
        population_task <- match(filters$task_id[population_rows], tasks$task_id)
        population_key <- paste(
            tasks$source_id[population_task], tasks$population_id[population_task],
            sep = "::"
        )
        population_signature <- paste(
            filters$metadata_scope[population_rows],
            filters$metadata_field[population_rows],
            filters$metadata_value[population_rows], sep = "\037"
        )
        forward <- split(population_signature, population_key)
        reverse <- split(population_key, paste(
            tasks$source_id[population_task], population_signature, sep = "::"
        ))
        binding_valid <- all(vapply(forward, function(x) {
            length(unique(x)) == 1L
        }, logical(1L))) && all(vapply(reverse, function(x) {
            length(unique(x)) == 1L
        }, logical(1L)))
        if (!binding_valid) {
            stop("Selection population bindings are invalid.", call. = FALSE)
        }
    }
    validate_binding <- function(rows, ids, label) {
        if (!any(rows)) return(invisible(TRUE))
        bound_task <- match(filters$task_id[rows], tasks$task_id)
        source <- tasks$source_id[bound_task]
        signature <- paste(
            filters$metadata_scope[rows], filters$metadata_field[rows],
            filters$metadata_value[rows], sep = "\037"
        )
        forward <- split(signature, paste(source, ids, sep = "::"))
        reverse <- split(paste(source, ids, sep = "::"), paste(
            source, signature, sep = "::"
        ))
        valid <- all(vapply(forward, function(x) {
            length(unique(x)) == 1L
        }, logical(1L))) && all(vapply(reverse, function(x) {
            length(unique(x)) == 1L
        }, logical(1L)))
        if (!valid) stop("Selection ", label, " bindings are invalid.",
            call. = FALSE)
        invisible(TRUE)
    }
    exposure_rows <- filters$filter_role == "exposure"
    exposure_task <- match(filters$task_id[exposure_rows], tasks$task_id)
    validate_binding(
        exposure_rows, tasks$exposure_id[exposure_task], "exposure"
    )
    time_rows <- filters$filter_role == "collection_time"
    time_task <- match(filters$task_id[time_rows], tasks$task_id)
    validate_binding(
        time_rows, tasks$collection_time_id[time_task], "collection-time"
    )
    condition_rows <- filters$filter_role == "condition"
    condition_task <- match(filters$task_id[condition_rows], tasks$task_id)
    condition_ids <- ifelse(
        filters$arm[condition_rows] == "intervention",
        tasks$intervention_condition_id[condition_task],
        tasks$reference_condition_id[condition_task]
    )
    validate_binding(condition_rows, condition_ids, "condition")
    invisible(filters)
}

validation_selection_validate_objects <- function(objects, sources) {
    validation_selection_identity(objects)
    validation_selection_ids(objects$object_id, "object_id")
    source_row <- match(objects$source_id, sources$source_id)
    order <- validation_selection_integer(objects$object_order, label = "object_order")
    kinds <- c(
        "molecular_values", "sample_metadata", "cell_metadata",
        "feature_metadata", "mapping_metadata", "combined_archive",
        "metadata_archive", "mapping_archive"
    )
    access <- c(
        molecular_values = "molecular_sealed",
        sample_metadata = "metadata_allowlisted",
        cell_metadata = "metadata_allowlisted",
        feature_metadata = "mapping_allowlisted",
        mapping_metadata = "mapping_allowlisted",
        combined_archive = "molecular_sealed",
        metadata_archive = "metadata_allowlisted",
        mapping_archive = "mapping_allowlisted"
    )
    valid <- all(order == seq_len(nrow(objects))) &&
        all(objects$object_id != "none") && !anyDuplicated(objects$object_id) &&
        !anyNA(source_row) && all(objects$object_kind %in% kinds) &&
        all(objects$access_class == unname(access[objects$object_kind])) &&
        setequal(unique(objects$source_id), sources$source_id)
    validation_selection_url(
        objects$request_locator, "planned object locator", opaque_resource = TRUE
    )
    locator_unique <- !anyDuplicated(validation_selection_url_identity(
        objects$request_locator
    ))
    molecular <- objects$object_kind %in% c(
        "molecular_values", "combined_archive"
    )
    coverage <- all(vapply(sources$source_id, function(id) {
        any(objects$source_id == id & molecular)
    }, logical(1L)))
    if (!valid || !locator_unique || !coverage) {
        stop("Selection planned objects are invalid.", call. = FALSE)
    }
    invisible(objects)
}

validation_selection_validate_inputs <- function(
    inputs, sources, objects, tasks, design, filters
) {
    validation_selection_identity(inputs)
    if (!identical(inputs$source_id, sources$source_id)) {
        stop("Selection assay-input source order is invalid.", call. = FALSE)
    }
    value_row <- match(inputs$values_object_id, objects$object_id)
    sample_row <- match(inputs$sample_metadata_object_id, objects$object_id)
    cell_present <- inputs$cell_metadata_object_id != "none"
    cell_row <- match(
        inputs$cell_metadata_object_id[cell_present], objects$object_id
    )
    mapping_present <- inputs$mapping_object_id != "none"
    mapping_row <- match(
        inputs$mapping_object_id[mapping_present], objects$object_id
    )
    source_match <- !anyNA(value_row) && !anyNA(sample_row) &&
        !anyNA(cell_row) && !anyNA(mapping_row) &&
        all(objects$source_id[value_row] == inputs$source_id) &&
        all(objects$source_id[sample_row] == inputs$source_id) &&
        all(objects$source_id[cell_row] == inputs$source_id[cell_present]) &&
        all(objects$source_id[mapping_row] == inputs$source_id[mapping_present]) &&
        all(objects$access_class[value_row] == "molecular_sealed") &&
        all(objects$access_class[sample_row] == "metadata_allowlisted") &&
        all(objects$access_class[cell_row] == "metadata_allowlisted") &&
        all(objects$access_class[mapping_row] == "mapping_allowlisted")
    member_valid <- function(member, container, rows, direct, archives) {
        kind <- objects$object_kind[rows]
        safe <- !grepl("(^|/)[.]{1,2}(/|$)|//|/$|^/|^[A-Za-z]:|\\\\",
            member) & !grepl("[*?]", member) & member != "whole_object" &
            member != "none"
        direct_valid <- kind %in% direct & member == "whole_object" &
            container %in% c("none", "gzip")
        archive_valid <- kind %in% archives & safe &
            container %in% c("zip", "tar", "tar_gzip")
        all(direct_valid | archive_valid)
    }
    members_valid <- member_valid(
        inputs$values_member, inputs$values_container, value_row,
        "molecular_values", "combined_archive"
    ) && member_valid(
        inputs$sample_metadata_member, inputs$sample_metadata_container,
        sample_row, "sample_metadata", "metadata_archive"
    ) && member_valid(
        inputs$cell_metadata_member[cell_present],
        inputs$cell_metadata_container[cell_present], cell_row,
        "cell_metadata", "metadata_archive"
    ) && member_valid(
        inputs$mapping_member[mapping_present],
        inputs$mapping_container[mapping_present], mapping_row,
        c("feature_metadata", "mapping_metadata"), "mapping_archive"
    )
    present_optional <- function(value, present) {
        all((value == "none") == !present)
    }
    no_wildcard <- c(
        inputs$observation_id_field, inputs$feature_id_field,
        inputs$quantity_field, inputs$sample_id_field,
        inputs$technical_repeat_id_field[
            inputs$technical_repeat_id_field != "none"
        ],
        inputs$cell_id_field[cell_present],
        inputs$cell_sample_id_field[cell_present],
        inputs$mapping_source_field[mapping_present],
        inputs$mapping_target_field[mapping_present]
    )
    adapters <- c(
        "delimited_wide", "delimited_long", "tenx_h5",
        "maxquant_protein_groups",
        "fragpipe_combined_protein"
    )
    metadata_adapters <- c(
        "tsv_header", "csv_header", "json_records", "parquet_table"
    )
    long <- inputs$matrix_orientation == "long_observation_feature"
    generic_wide <- inputs$values_adapter == "delimited_wide"
    tenx <- inputs$values_adapter == "tenx_h5"
    maxquant <- inputs$values_adapter == "maxquant_protein_groups"
    fragpipe <- inputs$values_adapter == "fragpipe_combined_protein"
    axis_selectors_valid <- all(
        inputs$observation_id_field != inputs$feature_id_field
    ) && all(!long |
        (inputs$observation_id_field != inputs$quantity_field &
            inputs$feature_id_field != inputs$quantity_field))
    quantity_selector_valid <- all(ifelse(
        long,
        inputs$quantity_selector == "exact_quantity_field",
        ifelse(
            generic_wide,
            inputs$quantity_selector == ifelse(
                inputs$matrix_orientation == "features_by_observations",
                "exact_observation_columns", "exact_observation_rows"
            ),
            ifelse(
                tenx,
                inputs$quantity_selector == "tenx_matrix",
                ifelse(
                    maxquant,
                    inputs$quantity_selector == ifelse(
                        inputs$quantity_measure == "maxlfq",
                        "maxquant_lfq_intensity_columns",
                        "maxquant_intensity_columns"
                    ),
                    inputs$quantity_selector == ifelse(
                        inputs$quantity_measure == "maxlfq",
                        "fragpipe_maxlfq_intensity_columns",
                        "fragpipe_intensity_columns"
                    )
                )
            )
        )
    ))
    scalar_valid <- source_match && members_valid &&
        axis_selectors_valid && quantity_selector_valid &&
        !any(grepl("[*?]", no_wildcard)) &&
        all(!no_wildcard %in% c("none", "not_applicable")) &&
        all(inputs$values_adapter %in% adapters) &&
        all(inputs$sample_metadata_adapter %in% metadata_adapters) &&
        all(inputs$cell_metadata_adapter[cell_present] %in%
            metadata_adapters) &&
        all(inputs$mapping_adapter[mapping_present] %in%
            metadata_adapters) &&
        all(inputs$matrix_orientation %in% c(
            "features_by_observations", "observations_by_features",
            "long_observation_feature"
        )) && all((inputs$values_adapter == "delimited_long") == long) &&
        all(!tenx | inputs$matrix_orientation == "features_by_observations") &&
        all(!(maxquant | fragpipe) |
            inputs$matrix_orientation == "features_by_observations") &&
        all((inputs$quantity_field == "matrix_values") == !long) &&
        all(inputs$observation_selection_rule %in% c(
            "exact_sample_ids_from_sample_metadata",
            "exact_cell_ids_from_cell_metadata"
        )) &&
        all(inputs$imputation_state == "none") &&
        all(inputs$quantity_level %in% c("gene", "protein")) &&
        all(inputs$quantity_scale %in% c(
            "raw_nonnegative_integer_counts", "linear_nonnegative"
        )) && all(inputs$normalization_state %in% c(
            "raw_counts", "raw_linear", "source_normalized_linear"
        ))
    validation_selection_url(inputs$evidence_locator, "assay-input evidence")
    assay <- sources$assay
    rna <- assay %in% c("bulk_rnaseq", "pseudobulk_rnaseq")
    pseudo <- assay == "pseudobulk_rnaseq"
    protein <- assay == "bulk_proteomics"
    tokens_valid <- vapply(seq_len(nrow(inputs)), function(i) {
        if (inputs$missing_semantics[[i]] == "missing_forbidden") {
            return(inputs$missing_tokens[[i]] == "none")
        }
        if (inputs$missing_semantics[[i]] != "explicit_tokens" ||
            inputs$missing_tokens[[i]] == "none") return(FALSE)
        tokens <- strsplit(inputs$missing_tokens[[i]], "|", fixed = TRUE)[[1L]]
        numeric_lexeme <- grepl(
            "^[+-]?(([0-9]+([.][0-9]*)?)|([.][0-9]+))([eE][+-]?[0-9]+)?$",
            tokens, perl = TRUE
        )
        numeric_value <- suppressWarnings(as.numeric(tokens))
        invalid_numeric_token <- numeric_lexeme &
            (!is.finite(numeric_value) | numeric_value == 0)
        identical(paste(tokens, collapse = "|"), inputs$missing_tokens[[i]]) &&
            all(nzchar(tokens)) && !anyDuplicated(tokens) &&
            !any(invalid_numeric_token)
    }, logical(1L))
    metadata_valid <- all(cell_present == pseudo) &&
        all(inputs$cell_metadata_member[!pseudo] == "none") &&
        all(inputs$cell_metadata_container[!pseudo] == "none") &&
        all(inputs$cell_metadata_adapter[!pseudo] == "none") &&
        all(inputs$cell_id_field[!pseudo] == "none") &&
        all(inputs$cell_sample_id_field[!pseudo] == "none") &&
        all(inputs$cell_id_field[pseudo] !=
            inputs$cell_sample_id_field[pseudo]) &&
        all(inputs$observation_selection_rule[pseudo] ==
            "exact_cell_ids_from_cell_metadata") &&
        all(inputs$observation_selection_rule[!pseudo] ==
            "exact_sample_ids_from_sample_metadata") &&
        all(!pseudo | paste(
            inputs$sample_metadata_object_id,
            inputs$sample_metadata_member, sep = "\037"
        ) != paste(
            inputs$cell_metadata_object_id,
            inputs$cell_metadata_member, sep = "\037"
        )) &&
        present_optional(inputs$mapping_member, mapping_present) &&
        all(inputs$mapping_container[!mapping_present] == "none") &&
        present_optional(inputs$mapping_adapter, mapping_present) &&
        present_optional(inputs$mapping_source_field, mapping_present) &&
        present_optional(inputs$mapping_target_field, mapping_present) &&
        all(!mapping_present |
            inputs$mapping_source_field != inputs$mapping_target_field)
    technical_present <- inputs$technical_repeat_id_field != "none"
    technical_valid <- all(technical_present ==
        (inputs$technical_repeat_aggregation != "none")) &&
        all(!technical_present |
            inputs$technical_repeat_id_field != inputs$sample_id_field) &&
        all(inputs$technical_repeat_aggregation[rna] %in% c(
            "none", "sum_raw_counts"
        )) && all(inputs$technical_repeat_aggregation[protein] %in% c(
            "none", "arithmetic_mean_linear",
            "source_nonoverlapping_linear_sum"
        ))
    technical_model_valid <- all(vapply(which(technical_present), function(i) {
        source_tasks <- tasks$task_id[tasks$source_id == inputs$source_id[[i]]]
        used_fields <- c(
            design$metadata_field[design$task_id %in% source_tasks],
            filters$metadata_field[
                filters$task_id %in% source_tasks &
                    filters$metadata_scope == "sample"
            ]
        )
        !inputs$technical_repeat_id_field[[i]] %in% used_fields
    }, logical(1L)))
    join_field_valid <- all(vapply(seq_len(nrow(inputs)), function(i) {
        source_tasks <- tasks$task_id[tasks$source_id == inputs$source_id[[i]]]
        model_fields <- design$metadata_field[
            design$task_id %in% source_tasks &
                design$field_role %in% c("pair", "block", "covariate")
        ]
        sample_filters <- filters$metadata_field[
            filters$task_id %in% source_tasks &
                filters$metadata_scope == "sample"
        ]
        cell_filters <- filters$metadata_field[
            filters$task_id %in% source_tasks & filters$metadata_scope == "cell"
        ]
        sample_valid <- !inputs$sample_id_field[[i]] %in% c(
            sample_filters, model_fields
        )
        cell_valid <- !cell_present[[i]] || !any(c(
            inputs$cell_id_field[[i]], inputs$cell_sample_id_field[[i]]
        ) %in% cell_filters)
        sample_valid && cell_valid
    }, logical(1L)))
    referenced <- c(
        inputs$values_object_id, inputs$sample_metadata_object_id,
        inputs$cell_metadata_object_id[cell_present],
        inputs$mapping_object_id[mapping_present]
    )
    reference_valid <- setequal(referenced, objects$object_id)
    assay_valid <- all(inputs$quantity_level[rna] == "gene") &&
        all(inputs$values_adapter[rna] %in% c(
            "delimited_wide", "delimited_long", "tenx_h5"
        )) &&
        all(inputs$quantity_measure[rna] == "raw_count") &&
        all(inputs$quantity_scale[rna] == "raw_nonnegative_integer_counts") &&
        all(inputs$normalization_state[rna] == "raw_counts") &&
        all(inputs$zero_semantics[rna] == "measured_zero") &&
        all(inputs$missing_semantics[rna] == "missing_forbidden") &&
        all(inputs$assignment_or_inference[rna] == "direct_gene_identifier") &&
        all(inputs$species_rule[rna] == "human_ensembl_only") &&
        all(inputs$mapping_rule[rna] == "strip_ensembl_version") &&
        all(inputs$group_rule[rna] == "not_applicable") &&
        all(inputs$duplicate_accession_operator[rna] == "not_applicable") &&
        all(inputs$quantity_level[protein] == "protein") &&
        all(inputs$values_adapter[protein] %in% c(
            "delimited_wide", "delimited_long",
            "maxquant_protein_groups", "fragpipe_combined_protein"
        )) &&
        all(inputs$quantity_measure[protein] %in% c(
            "intensity", "abundance", "area", "maxlfq"
        )) && all(!maxquant | inputs$quantity_measure %in% c(
            "intensity", "maxlfq"
        )) && all(!fragpipe | inputs$quantity_measure %in% c(
            "intensity", "maxlfq"
        )) && all(inputs$quantity_measure != "maxlfq" |
            inputs$normalization_state == "source_normalized_linear") &&
        all(inputs$quantity_scale[protein] == "linear_nonnegative") &&
        all(inputs$normalization_state[protein] %in% c(
            "raw_linear", "source_normalized_linear"
        )) && all(inputs$species_rule[protein] == "human_only_before_values") &&
        all(inputs$group_rule[protein] == "exclude_ambiguous_groups") &&
        all(inputs$mapping_rule[protein] ==
            "canonical_uniprot_strip_isoform") &&
        all(inputs$zero_semantics[protein] %in% c(
            "measured_zero", "zero_is_missing", "zero_forbidden"
        )) && all(inputs$missing_semantics[protein] %in% c(
            "missing_forbidden", "explicit_tokens"
        )) && all(inputs$assignment_or_inference[protein] %in% c(
            "direct_protein_accession", "source_protein_inference"
        )) && all(inputs$duplicate_accession_operator[protein] %in% c(
            "require_unique", "source_nonoverlapping_linear_sum"
        ))
    if (!scalar_valid || !metadata_valid || !technical_valid ||
        !technical_model_valid || !join_field_valid ||
        !reference_valid ||
        !all(tokens_valid) || !assay_valid) {
        stop("Selection assay-input contracts are invalid.", call. = FALSE)
    }
    invisible(inputs)
}

validation_selection_contrast_closure <- function(
    tasks, design, resampling, filters
) {
    design_key <- paste(
        design$task_id, "sample", design$metadata_field, sep = "\037"
    )
    filter_key <- paste(
        filters$task_id, filters$metadata_scope, filters$metadata_field,
        sep = "\037"
    )
    if (any(design_key %in% filter_key)) {
        stop("Selection model fields cannot be exact filter predicates.",
            call. = FALSE)
    }
    positive <- tasks$task_kind == "positive"
    groups <- split(
        tasks$task_id[positive],
        paste(tasks$source_id[positive], tasks$contrast_id[positive], sep = "::"),
        drop = TRUE
    )
    valid <- all(vapply(groups, function(ids) {
        design_signature <- vapply(ids, function(id) {
            rows <- design$task_id == id
            paste(paste(
                design$field_role[rows], design$metadata_field[rows],
                design$field_type[rows],
                sep = "\037"
            ), collapse = "\036")
        }, character(1L))
        resampling_signature <- vapply(ids, function(id) {
            rows <- resampling$task_id == id
            paste(sort(paste(
                resampling$resample_pool_id[rows], resampling$pool_role[rows],
                resampling$arm[rows], sep = "\037"
            )), collapse = "\036")
        }, character(1L))
        filter_signature <- vapply(ids, function(id) {
            rows <- filters$task_id == id
            paste(sort(paste(
                filters$arm[rows], filters$filter_role[rows],
                filters$metadata_scope[rows], filters$metadata_field[rows],
                filters$metadata_value[rows], sep = "\037"
            )), collapse = "\036")
        }, character(1L))
        length(unique(design_signature)) == 1L &&
            length(unique(resampling_signature)) == 1L &&
            length(unique(filter_signature)) == 1L
    }, logical(1L)))
    representatives <- vapply(groups, `[[`, character(1L), 1L)
    semantic <- vapply(representatives, function(id) {
        task <- tasks[tasks$task_id == id, , drop = FALSE]
        design_rows <- design$task_id == id &
            design$field_role %in% c("unit", "pair", "block")
        resampling_rows <- resampling$task_id == id
        filter_rows <- filters$task_id == id
        paste(
            task$source_id,
            paste(sort(paste(
                design$field_role[design_rows],
                design$metadata_field[design_rows],
                design$field_type[design_rows], sep = "\037"
            )), collapse = "\036"),
            paste(sort(paste(
                resampling$resample_pool_id[resampling_rows],
                resampling$pool_role[resampling_rows],
                resampling$arm[resampling_rows], sep = "\037"
            )), collapse = "\036"),
            paste(sort(paste(
                filters$arm[filter_rows], filters$filter_role[filter_rows],
                filters$metadata_scope[filter_rows],
                filters$metadata_field[filter_rows],
                filters$metadata_value[filter_rows], sep = "\037"
            )), collapse = "\036"), sep = "\035"
        )
    }, character(1L))
    task_membership_signature <- function(id) {
        task <- tasks[tasks$task_id == id, , drop = FALSE]
        unit_rows <- design$task_id == id & design$field_role == "unit"
        filter_rows <- filters$task_id == id
        paste(
            task$source_id,
            paste(sort(design$metadata_field[unit_rows]), collapse = "\036"),
            paste(sort(paste(
                filters$metadata_scope[filter_rows],
                filters$metadata_field[filter_rows],
                filters$metadata_value[filter_rows], sep = "\037"
            )), collapse = "\036"), sep = "\035"
        )
    }
    positive_membership <- vapply(
        representatives, task_membership_signature, character(1L)
    )
    paired_null <- tasks$task_kind == "null" &
        tasks$null_assignment_mode == "paired_swap"
    paired_semantic <- vapply(tasks$task_id[paired_null], function(id) {
        design_rows <- design$task_id == id &
            design$field_role %in% c("unit", "pair", "block")
        resampling_rows <- resampling$task_id == id
        filter_rows <- filters$task_id == id
        task <- tasks[tasks$task_id == id, , drop = FALSE]
        paste(
            task$source_id,
            paste(sort(paste(
                design$field_role[design_rows],
                design$metadata_field[design_rows],
                design$field_type[design_rows], sep = "\037"
            )), collapse = "\036"),
            paste(sort(paste(
                resampling$resample_pool_id[resampling_rows],
                resampling$pool_role[resampling_rows],
                resampling$arm[resampling_rows], sep = "\037"
            )), collapse = "\036"),
            paste(sort(paste(
                filters$arm[filter_rows], filters$filter_role[filter_rows],
                filters$metadata_scope[filter_rows],
                filters$metadata_field[filter_rows],
                filters$metadata_value[filter_rows], sep = "\037"
            )), collapse = "\036"), sep = "\035"
        )
    }, character(1L))
    paired_membership <- vapply(
        tasks$task_id[paired_null], task_membership_signature, character(1L)
    )
    if (!valid || anyDuplicated(semantic) ||
        anyDuplicated(positive_membership) ||
        anyDuplicated(paired_semantic) || anyDuplicated(paired_membership) ||
        any(semantic %in% paired_semantic) ||
        any(paired_membership %in% positive_membership)) {
        stop("Selection contrast closure is invalid.", call. = FALSE)
    }
    invisible(tasks)
}

validation_selection_locator_closure <- function(
    sources, tasks, inputs, objects
) {
    molecular <- objects$access_class == "molecular_sealed"
    molecular_locator <- validation_selection_url_identity(
        objects$request_locator[molecular]
    )
    cited <- tasks$citation_locator != "not_applicable"
    evidence_locator <- c(
        sources$metadata_locator, sources$license_locator,
        tasks$citation_locator[cited], inputs$evidence_locator
    )
    evidence_identity <- validation_selection_url_identity(evidence_locator)
    if (any(evidence_identity %in% molecular_locator)) {
        stop("Selection evidence aliases sealed molecular bytes.", call. = FALSE)
    }
    invisible(objects)
}

validation_selection_validate <- function(bundle_dir, root = ".") {
    validation_selection_require_amendment(root)
    roster <- validation_reactome97_read(root)
    paths <- validation_selection_bundle_paths(bundle_dir)
    sources <- validation_contract_read_tsv(
        paths[[1L]], VALIDATION_SELECTION_FIELDS$sources
    )
    tasks <- validation_contract_read_tsv(
        paths[[2L]], VALIDATION_SELECTION_FIELDS$tasks
    )
    design <- validation_contract_read_tsv(
        paths[[3L]], VALIDATION_SELECTION_FIELDS$design
    )
    resampling <- validation_contract_read_tsv(
        paths[[4L]], VALIDATION_SELECTION_FIELDS$resampling
    )
    filters <- validation_contract_read_tsv(
        paths[[5L]], VALIDATION_SELECTION_FIELDS$filters
    )
    inputs <- validation_contract_read_tsv(
        paths[[6L]], VALIDATION_SELECTION_FIELDS$inputs
    )
    objects <- validation_contract_read_tsv(
        paths[[7L]], VALIDATION_SELECTION_FIELDS$objects
    )
    validation_selection_validate_sources(sources)
    validation_selection_validate_tasks(tasks, sources, roster)
    validation_selection_validate_design(design, tasks)
    validation_selection_validate_filters(filters, tasks, sources)
    validation_selection_validate_resampling(
        resampling, tasks, design, filters
    )
    validation_selection_validate_objects(objects, sources)
    validation_selection_validate_inputs(
        inputs, sources, objects, tasks, design, filters
    )
    validation_selection_locator_closure(sources, tasks, inputs, objects)
    validation_selection_contrast_closure(
        tasks, design, resampling, filters
    )
    list(
        schema_version = VALIDATION_SELECTION_SCHEMA,
        protocol_version = VALIDATION_AMENDMENT_VERSION,
        protocol_sha256 = VALIDATION_AMENDMENT_SHA256,
        selection_sha256 = validation_selection_bundle_sha256(bundle_dir),
        eligible_targets = nrow(roster),
        sources = sources, tasks = tasks, design = design,
        resampling = resampling, filters = filters, inputs = inputs,
        objects = objects
    )
}

validation_byte_timestamp <- function(value) {
    pattern <- "^[0-9]{4}-[0-9]{2}-[0-9]{2}T[0-9]{2}:[0-9]{2}:[0-9]{2}Z$"
    parsed <- as.POSIXct(strptime(value, "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"))
    valid <- grepl(pattern, value) & !is.na(parsed) &
        format(parsed, "%Y-%m-%dT%H:%M:%SZ", tz = "UTC") == value
    if (!all(valid)) stop("Byte manifest UTC is invalid.", call. = FALSE)
    invisible(value)
}

validation_git_evidence <- function(
    root, arguments, stdout = FALSE, stderr = FALSE
) {
    command <- c(
        "-i",
        "LC_ALL=C", "LANG=C", "GIT_CONFIG_NOSYSTEM=1",
        "GIT_CONFIG_SYSTEM=/dev/null", "GIT_CONFIG_GLOBAL=/dev/null",
        "GIT_CONFIG_COUNT=0", "GIT_NO_REPLACE_OBJECTS=1",
        "GIT_ATTR_NOSYSTEM=1", "GIT_TERMINAL_PROMPT=0",
        "GIT_OPTIONAL_LOCKS=0", "GIT_NO_LAZY_FETCH=1",
        "/usr/bin/git", "--no-replace-objects",
        "-c", "core.fsmonitor=false", "-c", "core.untrackedCache=false",
        "-C", shQuote(root), vapply(
            arguments, shQuote, character(1L), USE.NAMES = FALSE
        )
    )
    result <- suppressWarnings(system2(
        "/usr/bin/env", command, stdout = stdout, stderr = stderr
    ))
    status <- attr(result, "status")
    if (is.null(status)) {
        status <- if (is.character(result)) 0L else result
    }
    list(
        status = as.integer(status),
        output = if (is.character(result)) unname(result) else character()
    )
}

validation_git_repository_validate <- function(root) {
    root <- normalizePath(root, mustWork = TRUE)
    top <- validation_git_evidence(
        root, c("rev-parse", "--show-toplevel"), stdout = TRUE
    )
    top_path <- if (top$status == 0L && length(top$output) == 1L) {
        suppressWarnings(normalizePath(top$output, mustWork = TRUE))
    } else {
        character()
    }
    shallow <- validation_git_evidence(
        root, c("rev-parse", "--is-shallow-repository"), stdout = TRUE
    )
    replacements <- validation_git_evidence(
        root,
        c("for-each-ref", "--format=%(refname)", "refs/replace/"),
        stdout = TRUE
    )
    unsafe_config <- validation_git_evidence(root, c(
        "config", "--local", "--get-regexp",
        "^(filter\\..*\\.(clean|smudge|process)|core\\.fsmonitor|extensions\\.(partialclone|worktreeconfig)|remote\\..*\\.promisor|include(if)?\\..*)$"
    ), stdout = TRUE)
    unsafe_paths <- vapply(c("info/grafts", "objects/info/alternates"),
        function(relative) {
            result <- validation_git_evidence(root, c(
                "rev-parse", "--path-format=absolute", "--git-path", relative
            ), stdout = TRUE)
            if (result$status != 0L || length(result$output) != 1L ||
                !startsWith(result$output, .Platform$file.sep)) return(TRUE)
            link <- Sys.readlink(result$output)
            file.exists(result$output) || (!is.na(link) && nzchar(link))
        }, logical(1L))
    valid <- length(top_path) == 1L && identical(top_path, root) &&
        shallow$status == 0L && identical(shallow$output, "false") &&
        replacements$status == 0L && length(replacements$output) == 0L &&
        unsafe_config$status %in% c(0L, 1L) &&
        length(unsafe_config$output) == 0L &&
        !any(unsafe_paths)
    if (!valid) {
        stop("Git evidence repository provenance is invalid.", call. = FALSE)
    }
    invisible(root)
}

validation_git_regular_blob <- function(root, commit, relative) {
    listing <- validation_git_evidence(root, c(
        "ls-tree", "--full-tree",
        "--format=%(objectmode)%x09%(objecttype)%x09%(path)",
        commit, "--", relative
    ), stdout = TRUE)
    listing$status == 0L && identical(
        listing$output, paste("100644", "blob", relative, sep = "\t")
    )
}

validation_byte_git_selection <- function(
    selection, bundle_dir, root, selection_commit
) {
    root <- normalizePath(root, mustWork = TRUE)
    validation_git_repository_validate(root)
    bundle <- normalizePath(bundle_dir, mustWork = TRUE)
    prefix <- paste0(root, .Platform$file.sep)
    if (!startsWith(bundle, prefix)) {
        stop("Selection bundle is outside the repository.", call. = FALSE)
    }
    relative <- substring(bundle, nchar(prefix) + 1L)
    commit_object <- paste0(selection_commit, "^{commit}")
    status <- validation_git_evidence(
        root, c("cat-file", "-e", commit_object)
    )
    ancestor <- validation_git_evidence(root, c(
        "merge-base", "--is-ancestor", selection_commit, "HEAD"
    ))
    if (status$status != 0L || ancestor$status != 0L) {
        stop("Selection commit is not an ancestor commit.", call. = FALSE)
    }
    tree_spec <- paste0(selection_commit, ":", relative)
    listing <- validation_git_evidence(root, c(
        "ls-tree", "--format=%(objectmode)%x09%(objecttype)%x09%(path)",
        tree_spec
    ), stdout = TRUE)
    expected_listing <- sort(paste(
        "100644", "blob", VALIDATION_SELECTION_FILES, sep = "\t"
    ), method = "radix")
    if (listing$status != 0L || !identical(
            sort(listing$output, method = "radix"), expected_listing
        )) {
        stop("Selection commit bundle tree closure is invalid.", call. = FALSE)
    }
    scratch <- tempfile("validation-selection-commit-")
    dir.create(scratch)
    on.exit(unlink(scratch, recursive = TRUE, force = TRUE), add = TRUE)
    for (name in VALIDATION_SELECTION_FILES) {
        spec <- paste0(selection_commit, ":", file.path(relative, name))
        output <- file.path(scratch, name)
        result <- validation_git_evidence(
            root, c("show", spec), stdout = output
        )
        if (result$status != 0L) {
            stop("Selection commit lacks registered bundle bytes.", call. = FALSE)
        }
    }
    observed <- validation_selection_bundle_sha256(scratch)
    if (!identical(observed, selection$selection_sha256)) {
        stop("Selection commit bytes do not match the manifest.", call. = FALSE)
    }
    invisible(selection_commit)
}

validation_byte_verify_objects <- function(selection, manifest, object_root) {
    root_info <- file.info(object_root)
    if (nrow(root_info) != 1L || is.na(root_info$isdir) || !root_info$isdir ||
        !validation_contract_is_type(object_root, "directory") ||
        nzchar(Sys.readlink(object_root))) {
        stop("Object root is not a regular directory.", call. = FALSE)
    }
    source_role <- setNames(
        selection$sources$access_role, selection$sources$source_id
    )
    object_source <- setNames(
        selection$objects$source_id, selection$objects$object_id
    )
    verified <- manifest$status == "verified"
    expected_files <- vapply(which(verified), function(i) {
        object_id <- manifest$object_id[[i]]
        role <- source_role[[object_source[[object_id]]]]
        file.path(role, paste0(object_id, ".bin"))
    }, character(1L))
    expected_dirs <- unique(dirname(expected_files))
    expected_dirs <- expected_dirs[expected_dirs != "."]
    expected_entries <- sort(c(expected_dirs, expected_files), method = "radix")
    observed_entries <- sort(list.files(
        object_root, all.files = TRUE, no.. = TRUE, recursive = TRUE,
        include.dirs = TRUE
    ), method = "radix")
    if (!identical(observed_entries, expected_entries)) {
        stop("Object root closure is invalid.", call. = FALSE)
    }
    if (length(observed_entries)) {
        paths <- file.path(object_root, observed_entries)
        info <- file.info(paths)
        expected_directory <- observed_entries %in% expected_dirs
        if (any(nzchar(Sys.readlink(paths))) || anyNA(info$isdir) ||
            !all(info$isdir == expected_directory) ||
            any(!expected_directory & !vapply(
                paths, validation_contract_is_type, logical(1L),
                expected = "regular file"
            )) || any(expected_directory & !vapply(
                paths, validation_contract_is_type, logical(1L),
                expected = "directory"
            ))) {
            stop("Object root contains an invalid entry.", call. = FALSE)
        }
    }
    for (i in seq_len(nrow(manifest))) {
        object_id <- manifest$object_id[[i]]
        role <- source_role[[object_source[[object_id]]]]
        path <- file.path(object_root, role, paste0(object_id, ".bin"))
        if (manifest$status[[i]] == "unavailable") {
            next
        }
        before <- file.info(path)
        if (!file.exists(path) || nrow(before) != 1L ||
            is.na(before$isdir) || before$isdir ||
            !validation_contract_is_type(path, "regular file") ||
            nzchar(Sys.readlink(path))) {
            stop("Verified object is not a regular file: ", object_id,
                call. = FALSE)
        }
        hash <- validation_contract_sha256(path)
        after <- file.info(path)
        stable_fields <- c("size", "mtime", "ctime", "ino")
        stable_fields <- intersect(stable_fields, names(before))
        stable <- identical(before[stable_fields], after[stable_fields])
        expected_size <- validation_selection_integer(
            manifest$bytes[[i]], zero = TRUE, label = "byte count"
        )
        if (!stable || before$size[[1L]] != expected_size ||
            !identical(hash, manifest$sha256[[i]])) {
            stop("Verified object bytes are invalid: ", object_id,
                call. = FALSE)
        }
    }
    invisible(manifest)
}

validation_byte_manifest_validate <- function(
    bundle_dir, manifest_path, root = ".", object_root = NULL,
    verify_git = TRUE
) {
    selection <- validation_selection_validate(bundle_dir, root)
    manifest <- validation_contract_read_tsv(manifest_path, VALIDATION_BYTE_FIELDS)
    valid_identity <- all(manifest$schema_version == VALIDATION_BYTE_SCHEMA) &&
        all(manifest$protocol_version == VALIDATION_AMENDMENT_VERSION) &&
        all(manifest$protocol_sha256 == VALIDATION_AMENDMENT_SHA256) &&
        length(unique(manifest$selection_commit)) == 1L &&
        all(grepl("^[0-9a-f]{40}$", manifest$selection_commit)) &&
        all(manifest$selection_sha256 == selection$selection_sha256) &&
        identical(manifest$object_id, selection$objects$object_id)
    if (!valid_identity) stop("Byte manifest identity is invalid.", call. = FALSE)
    validation_byte_timestamp(manifest$retrieved_utc)
    verified <- manifest$status == "verified"
    unavailable <- manifest$status == "unavailable"
    bytes <- validation_selection_integer(
        manifest$bytes, zero = TRUE, label = "byte count"
    )
    validation_selection_url(
        manifest$resolved_locator[verified], "resolved object locator",
        opaque_resource = TRUE
    )
    verified_valid <- all(
        grepl("^https://[^[:space:]]+$", manifest$resolved_locator[verified]) &
            manifest$terminal_http_status[verified] == "200" &
            bytes[verified] > 0 &
            grepl("^[0-9a-f]{64}$", manifest$sha256[verified]) &
            manifest$failure_code[verified] == "none"
    )
    failures <- c(
        "not_found", "access_denied", "license_blocked", "transport_error",
        "integrity_error"
    )
    terminal <- manifest$terminal_http_status[unavailable]
    terminal_valid <- terminal == "none" | grepl("^[1-5][0-9]{2}$", terminal)
    locator_valid <- ifelse(
        terminal == "none",
        manifest$resolved_locator[unavailable] == "none",
        grepl("^https://[^[:space:]]+$", manifest$resolved_locator[unavailable])
    )
    response_locator <- unavailable & manifest$terminal_http_status != "none"
    validation_selection_url(
        manifest$resolved_locator[response_locator], "failed object locator",
        opaque_resource = TRUE
    )
    failure <- manifest$failure_code[unavailable]
    failure_valid <-
        (terminal == "none" & failure %in% c(
            "license_blocked", "transport_error"
        )) |
        (terminal %in% c("404", "410") & failure == "not_found") |
        (terminal %in% c("401", "403") & failure %in% c(
            "access_denied", "license_blocked"
        )) |
        (terminal == "451" & failure == "license_blocked") |
        (grepl("^2[0-9]{2}$", terminal) & failure == "integrity_error") |
        ((grepl("^[135][0-9]{2}$", terminal) | terminal == "429") &
            failure == "transport_error") |
        (grepl("^4[0-9]{2}$", terminal) &
            !terminal %in% c("401", "403", "404", "410", "429", "451") &
            failure == "access_denied")
    unavailable_valid <- all(
        terminal_valid & locator_valid & failure_valid &
            bytes[unavailable] == 0 & manifest$sha256[unavailable] == "none" &
            manifest$failure_code[unavailable] %in% failures
    )
    resolved <- manifest$resolved_locator != "none"
    resolved_identity <- validation_selection_url_identity(
        manifest$resolved_locator[resolved]
    )
    resolved_unique <- !anyDuplicated(resolved_identity)
    planned_identity <- validation_selection_url_identity(
        selection$objects$request_locator
    )
    resolved_rows <- which(resolved)
    planned_match <- match(resolved_identity, planned_identity)
    terminal_planned_valid <- all(
        is.na(planned_match) | planned_match == resolved_rows
    )
    object_source_row <- match(
        selection$objects$source_id, selection$sources$source_id
    )
    object_role <- selection$sources$access_role[object_source_row]
    object_group <- selection$sources$independence_group[object_source_row]
    object_epoch <- ifelse(
        selection$objects$access_class == "molecular_sealed" &
            object_role == "heldout",
        "post_implementation", "pre_implementation"
    )
    hash_rows <- split(which(verified), manifest$sha256[verified])
    alias_valid <- all(vapply(hash_rows, function(rows) {
        classes <- selection$objects$access_class[rows]
        molecular <- classes == "molecular_sealed"
        class_safe <- !(any(molecular) && any(!molecular))
        epoch_safe <- length(unique(object_epoch[rows])) == 1L
        group_safe <- !any(molecular) ||
            length(unique(object_group[rows][molecular])) == 1L
        class_safe && epoch_safe && group_safe
    }, logical(1L)))
    if (!all(verified | unavailable) || !verified_valid || !unavailable_valid ||
        !resolved_unique || !terminal_planned_valid || !alias_valid) {
        stop("Byte manifest rows are invalid.", call. = FALSE)
    }
    if (any(verified) && is.null(object_root)) {
        stop("Verified byte rows require the exact object root.", call. = FALSE)
    }
    if (verify_git) validation_byte_git_selection(
        selection, bundle_dir, root, unique(manifest$selection_commit)
    )
    if (!is.null(object_root)) {
        validation_byte_verify_objects(selection, manifest, object_root)
    }
    list(
        schema_version = VALIDATION_BYTE_SCHEMA,
        protocol_version = VALIDATION_AMENDMENT_VERSION,
        selection_sha256 = selection$selection_sha256,
        selection_commit = unique(manifest$selection_commit),
        objects = nrow(manifest), verified = sum(verified),
        unavailable = sum(unavailable)
    )
}

validation_access_expected <- function(selection, manifest) {
    source_role <- setNames(
        selection$sources$access_role, selection$sources$source_id
    )
    object_source <- selection$objects$source_id
    role <- unname(source_role[object_source])
    verified <- manifest$status == "verified"
    allowlisted <- selection$objects$access_class %in% c(
        "metadata_allowlisted", "mapping_allowlisted"
    )
    molecular <- selection$objects$access_class == "molecular_sealed"
    first <- which(verified & (allowlisted | (molecular & role == "development")))
    heldout <- which(verified & molecular & role == "heldout")
    first_type <- ifelse(
        allowlisted[first], "metadata_first_access",
        "development_first_content_access"
    )
    types <- c(
        "acquisition_start",
        rep("object_stream_verification", nrow(manifest)),
        "byte_manifest_freeze", first_type, "implementation_freeze",
        rep("heldout_preopen_verification", length(heldout)),
        "execution_start",
        rep("heldout_first_content_access", length(heldout))
    )
    object_index <- c(
        NA_integer_, seq_len(nrow(manifest)), NA_integer_, first,
        NA_integer_, heldout, NA_integer_, heldout
    )
    object <- rep("none", length(types))
    bytes <- rep("0", length(types))
    sha256 <- rep("none", length(types))
    present <- !is.na(object_index)
    object[present] <- manifest$object_id[object_index[present]]
    bytes[present] <- manifest$bytes[object_index[present]]
    sha256[present] <- manifest$sha256[object_index[present]]
    list(
        event_type = types, object_id = object, bytes = bytes, sha256 = sha256,
        stream_rows = 1L + seq_len(nrow(manifest)),
        byte_freeze = nrow(manifest) + 2L,
        implementation_freeze = nrow(manifest) + 3L + length(first),
        first_rows = if (length(first)) {
            nrow(manifest) + 2L + seq_along(first)
        } else integer(),
        preopen_rows = if (length(heldout)) {
            nrow(manifest) + 3L + length(first) + seq_along(heldout)
        } else integer(),
        execution_start = nrow(manifest) + 4L + length(first) + length(heldout),
        heldout_rows = if (length(heldout)) {
            nrow(manifest) + 4L + length(first) + length(heldout) +
                seq_along(heldout)
        } else integer()
    )
}

validation_access_git_artifacts <- function(
    selection, bundle_dir, manifest_path, root, selection_commit,
    byte_commit, implementation_commit, byte_sha256, closure_sha256,
    object_root
) {
    root <- normalizePath(root, mustWork = TRUE)
    validation_git_repository_validate(root)
    live_head <- validation_git_evidence(
        root, c("rev-parse", "--verify", "HEAD"), stdout = TRUE
    )
    if (live_head$status != 0L || length(live_head$output) != 1L ||
        !identical(live_head$output, implementation_commit)) {
        stop("Production checkout is not the logged implementation commit.",
            call. = FALSE)
    }
    git_clean <- function() {
        status_output <- tempfile("validation-git-status-")
        on.exit(unlink(status_output), add = TRUE)
        status <- validation_git_evidence(root, c(
            "status", "--porcelain=v1", "--untracked-files=all",
            "--ignore-submodules=none"
        ), stdout = status_output)
        status_size <- file.info(status_output)$size[[1L]]
        status$status == 0L && !is.na(status_size) && status_size == 0
    }
    if (!git_clean()) {
        stop("Production access validation requires a clean worktree.",
            call. = FALSE)
    }
    manifest <- normalizePath(manifest_path, mustWork = TRUE)
    prefix <- paste0(root, .Platform$file.sep)
    if (!startsWith(manifest, prefix)) {
        stop("Byte manifest is outside the repository.", call. = FALSE)
    }
    manifest_relative <- substring(manifest, nchar(prefix) + 1L)
    strict_ancestor <- function(ancestor, descendant) {
        !identical(ancestor, descendant) &&
            validation_git_evidence(root, c(
                "merge-base", "--is-ancestor", ancestor, descendant
            ))$status == 0L
    }
    if (!strict_ancestor(selection_commit, byte_commit) ||
        !strict_ancestor(byte_commit, implementation_commit)) {
        stop("Access-event commit ancestry is invalid.", call. = FALSE)
    }
    protocol <- validation_effective_protocol_read(root)
    validator_relative <- validation_effective_value(
        protocol, "selection", "validator_path"
    )
    validator_sha256 <- validation_effective_value(
        protocol, "selection", "validator_sha256"
    )
    runtime_paths <- strsplit(validation_effective_value(
        protocol, "selection", "runtime_paths"
    ), ",", fixed = TRUE)[[1L]]
    commits <- c(selection_commit, byte_commit, implementation_commit)
    runtime_hashes <- vapply(runtime_paths, function(relative) {
        live_path <- file.path(root, relative)
        if (!validation_contract_is_type(live_path, "regular file") ||
            nzchar(Sys.readlink(live_path))) {
            stop("Semantic runtime path is invalid: ", relative, call. = FALSE)
        }
        committed <- vapply(commits, function(commit) {
            output <- tempfile("validation-runtime-commit-")
            on.exit(unlink(output), add = TRUE)
            if (!validation_git_regular_blob(root, commit, relative)) {
                return(NA_character_)
            }
            result <- validation_git_evidence(
                root, c("show", paste0(commit, ":", relative)), stdout = output
            )
            if (result$status != 0L) return(NA_character_)
            validation_contract_sha256(output)
        }, character(1L))
        live_hash <- validation_contract_sha256(live_path)
        if (anyNA(committed) || length(unique(committed)) != 1L ||
            !identical(live_hash, committed[[1L]])) {
            stop("Semantic runtime changed after selection freeze: ", relative,
                call. = FALSE)
        }
        live_hash
    }, character(1L))
    if (!identical(runtime_paths, VALIDATION_RUNTIME_PATHS) ||
        !identical(validator_relative, VALIDATION_SELECTION_VALIDATOR_PATH) ||
        !identical(
            unname(runtime_hashes[[match(validator_relative, runtime_paths)]]),
            validator_sha256
        )) {
        stop("Committed semantic runtime identity is invalid.", call. = FALSE)
    }
    bundle <- normalizePath(bundle_dir, mustWork = TRUE)
    if (!startsWith(bundle, prefix)) {
        stop("Selection bundle is outside the repository.", call. = FALSE)
    }
    bundle_relative <- substring(bundle, nchar(prefix) + 1L)
    selection_tree <- vapply(commits, function(commit) {
        result <- validation_git_evidence(root, c(
            "rev-parse", paste0(commit, ":", bundle_relative)
        ), stdout = TRUE)
        if (result$status != 0L || length(result$output) != 1L) {
            return(NA_character_)
        }
        result$output
    }, character(1L))
    if (anyNA(selection_tree) || length(unique(selection_tree)) != 1L) {
        stop("Selection bundle Git tree changed after selection freeze.",
            call. = FALSE)
    }
    for (commit in c(byte_commit, implementation_commit)) {
        validation_byte_git_selection(selection, bundle_dir, root, commit)
        scratch <- tempfile("validation-byte-commit-")
        on.exit(unlink(scratch), add = TRUE)
        spec <- paste0(commit, ":", manifest_relative)
        if (!validation_git_regular_blob(root, commit, manifest_relative)) {
            stop("Committed byte manifest is not a regular blob.", call. = FALSE)
        }
        status <- validation_git_evidence(
            root, c("show", spec), stdout = scratch
        )
        if (status$status != 0L ||
            !identical(validation_contract_sha256(scratch), byte_sha256)) {
            stop("Committed byte manifest identity is invalid.", call. = FALSE)
        }
        unlink(scratch)
    }
    closure_relative <- "benchmark/validation-execution-closure.R"
    closure_path <- file.path(root, closure_relative)
    if (!validation_contract_is_type(closure_path, "regular file") ||
        nzchar(Sys.readlink(closure_path)) ||
        !identical(validation_contract_sha256(closure_path), closure_sha256)) {
        stop("F-E-1.0.0 implementation closure validator is required.",
            call. = FALSE)
    }
    closure_commit <- tempfile("validation-execution-closure-")
    on.exit(unlink(closure_commit), add = TRUE)
    closure_spec <- paste0(implementation_commit, ":", closure_relative)
    if (!validation_git_regular_blob(
            root, implementation_commit, closure_relative
        )) {
        stop("Committed F-E-1.0.0 validator is not a regular blob.",
            call. = FALSE)
    }
    closure_status <- validation_git_evidence(
        root, c("show", closure_spec), stdout = closure_commit
    )
    if (closure_status$status != 0L ||
        !identical(validation_contract_sha256(closure_commit), closure_sha256)) {
        stop("Committed F-E-1.0.0 validator identity is invalid.", call. = FALSE)
    }
    prior_closure_final <- any(vapply(c(selection_commit, byte_commit),
        function(commit) {
            prior <- tempfile("validation-prior-execution-closure-")
            on.exit(unlink(prior), add = TRUE)
            result <- validation_git_evidence(root, c(
                "show", paste0(commit, ":", closure_relative)
            ), stdout = prior)
            result$status == 0L && identical(
                validation_contract_sha256(prior), closure_sha256
            )
        }, logical(1L)))
    if (prior_closure_final) {
        stop("Final F-E-1.0.0 closure predates implementation freeze.",
            call. = FALSE)
    }
    closure_environment <- new.env(parent = baseenv())
    sys.source(closure_commit, envir = closure_environment)
    if (!exists(
            "validation_execution_closure_validate",
            envir = closure_environment, mode = "function", inherits = FALSE
        )) {
        stop("Committed F-E-1.0.0 validator interface is invalid.", call. = FALSE)
    }
    closure_validator <- closure_environment$validation_execution_closure_validate
    formal_values <- formals(closure_validator)
    formals_valid <- identical(names(formal_values), c(
        "root", "implementation_commit", "validator_sha256"
    )) && all(vapply(formal_values, function(value) {
        identical(value, quote(expr = ))
    }, logical(1L)))
    if (!formals_valid) {
        stop("Committed F-E-1.0.0 validator formals are invalid.", call. = FALSE)
    }
    result <- closure_validator(
        root, implementation_commit, closure_sha256
    )
    result_valid <- identical(result, list(
        schema_version = "F-E-1.0.0",
        implementation_commit = implementation_commit,
        validator_sha256 = closure_sha256,
        executable = TRUE
    ))
    if (!result_valid) {
        stop("F-E-1.0.0 implementation closure result is invalid.", call. = FALSE)
    }
    validation_byte_verify_objects(
        selection,
        validation_contract_read_tsv(manifest_path, VALIDATION_BYTE_FIELDS),
        object_root
    )
    final_runtime_hashes <- vapply(runtime_paths, function(relative) {
        path <- file.path(root, relative)
        if (!validation_contract_is_type(path, "regular file") ||
            nzchar(Sys.readlink(path))) return(NA_character_)
        validation_contract_sha256(path)
    }, character(1L))
    artifacts_stable <- !anyNA(final_runtime_hashes) &&
        identical(unname(final_runtime_hashes), unname(runtime_hashes)) &&
        identical(validation_contract_sha256(manifest_path), byte_sha256) &&
        identical(
            validation_selection_bundle_sha256(bundle_dir),
            selection$selection_sha256
        ) && identical(
            validation_contract_sha256(closure_path), closure_sha256
        )
    if (!artifacts_stable) {
        stop("Production artifacts changed during access validation.",
            call. = FALSE)
    }
    validation_git_repository_validate(root)
    final_head <- validation_git_evidence(
        root, c("rev-parse", "--verify", "HEAD"), stdout = TRUE
    )
    if (!git_clean() || final_head$status != 0L ||
        !identical(final_head$output, implementation_commit)) {
        stop("Production checkout changed during access validation.",
            call. = FALSE)
    }
    invisible(implementation_commit)
}

validation_access_log_validate <- function(
    bundle_dir, manifest_path, access_path, root = ".", object_root = NULL,
    verify_git = TRUE
) {
    if (verify_git && is.null(object_root)) {
        stop("Production access validation requires the exact object root.",
            call. = FALSE)
    }
    validation_byte_manifest_validate(
        bundle_dir, manifest_path, root, object_root, verify_git
    )
    selection <- validation_selection_validate(bundle_dir, root)
    manifest <- validation_contract_read_tsv(manifest_path, VALIDATION_BYTE_FIELDS)
    access <- validation_contract_read_tsv(access_path, VALIDATION_ACCESS_FIELDS)
    expected <- validation_access_expected(selection, manifest)
    event_order <- validation_selection_integer(
        access$event_order, label = "access event order"
    )
    validation_byte_timestamp(access$event_utc)
    event_time <- as.numeric(as.POSIXct(
        strptime(access$event_utc, "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
    ))
    identity_valid <- nrow(access) == length(expected$event_type) &&
        all(access$schema_version == VALIDATION_ACCESS_SCHEMA) &&
        all(access$protocol_version == VALIDATION_AMENDMENT_VERSION) &&
        all(access$protocol_sha256 == VALIDATION_AMENDMENT_SHA256) &&
        all(event_order == seq_len(nrow(access))) &&
        identical(access$event_type, expected$event_type) &&
        identical(access$object_id, expected$object_id) &&
        identical(access$bytes, expected$bytes) &&
        identical(access$sha256, expected$sha256) &&
        all(access$selection_sha256 == selection$selection_sha256) &&
        all(access$worktree_state == "clean") &&
        length(unique(access$execution_attempt_id)) == 1L &&
        all(grepl("^attempt_[0-9a-f]{32}$", access$execution_attempt_id)) &&
        all(diff(event_time) >= 0) &&
        identical(
            access$event_utc[expected$stream_rows], manifest$retrieved_utc
        )
    if (!identity_valid) {
        stop("Access-event state closure is invalid.", call. = FALSE)
    }

    selection_commit <- unique(manifest$selection_commit)
    byte_sha256 <- validation_contract_sha256(manifest_path)
    byte_commit <- access$repository_head[[expected$byte_freeze]]
    implementation_commit <-
        access$repository_head[[expected$implementation_freeze]]
    commit_valid <- grepl("^[0-9a-f]{40}$", byte_commit) &&
        grepl("^[0-9a-f]{40}$", implementation_commit) &&
        !identical(selection_commit, byte_commit) &&
        !identical(byte_commit, implementation_commit)
    before_byte <- seq_len(expected$byte_freeze - 1L)
    from_byte <- expected$byte_freeze:nrow(access)
    before_implementation <- seq_len(expected$implementation_freeze - 1L)
    from_implementation <- expected$implementation_freeze:nrow(access)
    closure_sha256 <- unique(access$execution_closure_sha256[from_implementation])
    head_valid <- all(access$repository_head[before_byte] == selection_commit) &&
        all(access$repository_head[
            expected$byte_freeze:(expected$implementation_freeze - 1L)
        ] == byte_commit) &&
        all(access$repository_head[from_implementation] == implementation_commit)
    artifact_valid <- all(access$byte_manifest_sha256[before_byte] == "none") &&
        all(access$byte_manifest_sha256[from_byte] == byte_sha256) &&
        all(access$implementation_commit[before_implementation] == "none") &&
        all(access$implementation_commit[from_implementation] ==
            implementation_commit) &&
        all(access$execution_closure_sha256[before_implementation] == "none") &&
        length(closure_sha256) == 1L &&
        grepl("^[0-9a-f]{64}$", closure_sha256) &&
        all(access$execution_closure_sha256[from_implementation] == closure_sha256)
    if (!commit_valid || !head_valid || !artifact_valid) {
        stop("Access-event artifact closure is invalid.", call. = FALSE)
    }
    if (verify_git) validation_access_git_artifacts(
        selection, bundle_dir, manifest_path, root, selection_commit,
        byte_commit, implementation_commit, byte_sha256, closure_sha256,
        object_root
    )
    list(
        schema_version = VALIDATION_ACCESS_SCHEMA,
        protocol_version = VALIDATION_AMENDMENT_VERSION,
        access_sha256 = validation_contract_sha256(access_path),
        events = nrow(access), objects = nrow(manifest),
        attempt_id = unique(access$execution_attempt_id), executable = verify_git
    )
}
