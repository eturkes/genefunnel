# Assisted-by: OpenAI Codex.

VALIDATION_AMENDMENT_VERSION <- "F-2.1.0"
VALIDATION_AMENDMENT_SHA256 <-
    "7bda0577d077a5204165bd570af6bd6eaeeb2d581571e17ff54e832b4f2cf2a5"
VALIDATION_AMENDMENT_PARENT_VERSION <- "F-2.0.0"
VALIDATION_AMENDMENT_PARENT_SHA256 <-
    "00a3a9f5e1cd6f709fc3ecdbad111f938b80c849c1c6071d134db41e7c63361c"

validation_contract_sha256 <- function(path) {
    namespace <- asNamespace("tools")
    if (!exists("sha256sum", envir = namespace, inherits = FALSE)) {
        stop("Validation contracts require tools::sha256sum().", call. = FALSE)
    }
    unname(get("sha256sum", envir = namespace, inherits = FALSE)(path))
}

validation_contract_is_type <- function(path, expected) {
    observed <- suppressWarnings(system2(
        "/usr/bin/stat", c("--format=%F", "--", shQuote(path)), stdout = TRUE,
        stderr = FALSE, env = "LC_ALL=C"
    ))
    if (!identical(attr(observed, "status"), NULL) || length(observed) != 1L) {
        return(FALSE)
    }
    accepted <- if (expected == "regular file") {
        c("regular file", "regular empty file")
    } else {
        expected
    }
    unname(observed) %in% accepted
}

validation_contract_read_tsv <- function(path, fields, allow_empty = FALSE) {
    info <- file.info(path)
    if (nrow(info) != 1L || is.na(info$isdir) || info$isdir ||
        !validation_contract_is_type(path, "regular file") ||
        nzchar(Sys.readlink(path))) {
        stop("Validation contract is not a regular file: ", path,
            call. = FALSE)
    }
    size <- info$size[[1L]]
    connection <- file(path, open = "rb")
    on.exit(close(connection), add = TRUE)
    bytes <- readBin(connection, what = "raw", n = size)
    if (length(bytes) != size || size == 0 || bytes[[size]] != as.raw(10L) ||
        any(bytes == as.raw(0L)) || any(bytes == as.raw(13L)) ||
        (size >= 3L && identical(bytes[1:3], as.raw(c(239L, 187L, 191L))))) {
        stop("Validation contract bytes are noncanonical: ", path,
            call. = FALSE)
    }
    byte_values <- as.integer(bytes)
    forbidden_controls <- byte_values < 32L &
        !byte_values %in% c(9L, 10L)
    if (any(forbidden_controls) || any(byte_values == 127L)) {
        stop("Validation contract contains control bytes: ", path,
            call. = FALSE)
    }
    text <- rawToChar(bytes)
    if (!validUTF8(text) || grepl("\n\n", text, fixed = TRUE)) {
        stop("Validation contract text is noncanonical: ", path,
            call. = FALSE)
    }
    lines <- strsplit(sub("\n$", "", text), "\n", fixed = TRUE)[[1L]]
    tab_count <- nchar(lines, type = "bytes") - nchar(
        gsub("\t", "", lines, fixed = TRUE), type = "bytes"
    )
    if ((!allow_empty && length(lines) < 2L) ||
        any(tab_count != length(fields) - 1L)) {
        stop("Validation contract field closure is invalid: ", path,
            call. = FALSE)
    }
    table <- utils::read.delim(
        text = text, stringsAsFactors = FALSE, check.names = FALSE,
        colClasses = "character", na.strings = character(), quote = "",
        comment.char = "", row.names = NULL
    )
    if (!identical(names(table), fields) ||
        (!allow_empty && nrow(table) == 0L) || anyNA(table) ||
        nrow(table) != length(lines) - 1L ||
        !identical(row.names(table), as.character(seq_len(nrow(table)))) ||
        any(!nzchar(unlist(table, use.names = FALSE)))) {
        stop("Validation contract table is invalid: ", path, call. = FALSE)
    }
    table
}

validation_amendment_read <- function(root = ".") {
    path <- file.path(root, "benchmark", "validation-amendment.tsv")
    observed <- validation_contract_sha256(path)
    if (!identical(observed, VALIDATION_AMENDMENT_SHA256)) {
        stop("Validation amendment identity is invalid.", call. = FALSE)
    }
    amendment <- validation_contract_read_tsv(
        path, c("protocol_version", "operation", "section", "key", "value")
    )
    identity <- paste(amendment$section, amendment$key, sep = "::")
    expected_sections <- c(
        identity = 11L, environment = 2L, catalogue = 4L, assay = 4L,
        eligibility = 7L, split = 8L, preprocess_proteomics = 3L,
        task = 5L, primary = 14L, negative_control = 14L, output = 4L,
        selection = 25L, bytes = 8L, access = 14L, simulation = 1L,
        secondary = 1L, decision = 3L
    )
    section_counts <- table(factor(
        amendment$section, levels = names(expected_sections)
    ))
    valid <- nrow(amendment) == 128L &&
        all(amendment$protocol_version == VALIDATION_AMENDMENT_VERSION) &&
        identical(as.integer(table(factor(
            amendment$operation, levels = c("add", "replace")
        ))), c(82L, 46L)) &&
        identical(as.integer(section_counts), unname(expected_sections)) &&
        all(amendment$section %in% names(expected_sections)) &&
        all(grepl("^[a-z][a-z0-9_]*$", amendment$section)) &&
        all(grepl("^[A-Za-z][A-Za-z0-9_]*$", amendment$key)) &&
        !anyDuplicated(identity)
    if (!valid) stop("Validation amendment rows are invalid.", call. = FALSE)
    amendment
}

validation_amendment_parent_read <- function(root = ".") {
    environment <- new.env(parent = baseenv())
    sys.source(
        file.path(root, "benchmark", "validation-protocol.R"),
        envir = environment
    )
    if (!identical(environment$VALIDATION_PROTOCOL_VERSION,
            VALIDATION_AMENDMENT_PARENT_VERSION) ||
        !identical(environment$VALIDATION_PROTOCOL_SHA256,
            VALIDATION_AMENDMENT_PARENT_SHA256)) {
        stop("Validation amendment parent constants are invalid.", call. = FALSE)
    }
    environment$validation_protocol_read(root)
}

validation_effective_protocol_read <- function(root = ".") {
    parent <- validation_amendment_parent_read(root)
    amendment <- validation_amendment_read(root)
    parent_identity <- paste(parent$section, parent$key, sep = "::")
    amendment_identity <- paste(amendment$section, amendment$key, sep = "::")
    replace <- amendment$operation == "replace"
    add <- amendment$operation == "add"
    if (any(!amendment_identity[replace] %in% parent_identity) ||
        any(amendment_identity[add] %in% parent_identity)) {
        stop("Validation amendment operations are invalid.", call. = FALSE)
    }
    effective <- parent
    effective$protocol_version <- VALIDATION_AMENDMENT_VERSION
    for (i in which(replace)) {
        row <- match(amendment_identity[[i]], parent_identity)
        effective$value[[row]] <- amendment$value[[i]]
    }
    additions <- amendment[add, c(
        "protocol_version", "section", "key", "value"
    ), drop = FALSE]
    rownames(additions) <- NULL
    rownames(effective) <- NULL
    effective <- rbind(effective, additions)
    identity <- paste(effective$section, effective$key, sep = "::")
    if (nrow(effective) != 266L || anyDuplicated(identity) ||
        !all(effective$protocol_version == VALIDATION_AMENDMENT_VERSION)) {
        stop("Effective validation protocol is invalid.", call. = FALSE)
    }
    effective
}

validation_effective_value <- function(protocol, section, key) {
    selected <- protocol$section == section & protocol$key == key
    if (sum(selected) != 1L) {
        stop("Effective validation key is invalid: ", section, "/", key,
            call. = FALSE)
    }
    protocol$value[[which(selected)]]
}

validation_effective_protocol_validate <- function(root = ".") {
    protocol <- validation_effective_protocol_read(root)
    exact <- list(
        c("identity", "parent_protocol_version", VALIDATION_AMENDMENT_PARENT_VERSION),
        c("identity", "parent_protocol_sha256", VALIDATION_AMENDMENT_PARENT_SHA256),
        c("environment", "bootstrap_repeats",
            "1000000_shared_draws_for_all_methods_endpoints_and_dependent_tasks"),
        c("primary", "utility",
            "GeneFunnel_fixed_panel_stability_quantile_at_least_0.55"),
        c("primary", "families",
            "5_fixed_panel_gates_per_assay_times_3_assays_equals_15_nominal_tail_allocations"),
        c("primary", "estimand",
            "rule_derived_realized_eligible_subset_of_the_finite_prospectively_selected_candidate_source_contrast_target_and_pipeline_roster"),
        c("primary", "bound_interpretation",
            "type8_fixed_panel_positive_cluster_weight_stability_quantile_conditional_on_realized_QC_coverage_support_matching_score_matrices_and_controls;not_a_confidence_bound_or_error_rate_guarantee"),
        c("catalogue", "eligible_target_roster_sha256",
            "915fce5e43b0e0fe57baf1a0b6cf646f5be8c4e25d325612c5ad6263ebc69202"),
        c("catalogue", "eligible_target_generator_sha256",
            "b53399f9fe942fdda84d14eaffadfd547dc8bee0329f9db7e3f5d164b8b03b77"),
        c("selection", "schema_version", "F-S-1.0.0"),
        c("bytes", "schema_version", "F-B-1.0.0"),
        c("access", "schema_version", "F-A-1.0.0"),
        c("decision", "assay_promotion",
            "fixed_panel_pass_report_only_for_an_assay_whose_complete_primary_and_consistency_gates_pass_without_population_or_error_rate_inference")
    )
    exact_valid <- all(vapply(exact, function(row) {
        identical(validation_effective_value(
            protocol, row[[1L]], row[[2L]]
        ), row[[3L]])
    }, logical(1L)))
    selection_files <- strsplit(validation_effective_value(
        protocol, "selection", "bundle_files"
    ), ",", fixed = TRUE)[[1L]]
    if (!exact_valid || !identical(selection_files, c(
        "sources.tsv", "tasks.tsv", "task-design.tsv", "task-resampling.tsv",
        "task-filters.tsv", "assay-inputs.tsv", "objects-planned.tsv"
    ))) {
        stop("Effective validation semantics are invalid.", call. = FALSE)
    }
    list(
        protocol_version = VALIDATION_AMENDMENT_VERSION,
        protocol_sha256 = VALIDATION_AMENDMENT_SHA256,
        parent_version = VALIDATION_AMENDMENT_PARENT_VERSION,
        parent_sha256 = VALIDATION_AMENDMENT_PARENT_SHA256,
        amendment_rows = 128L,
        effective_rows = nrow(protocol),
        bootstrap_repeats = 1000000L,
        selection_schema = "F-S-1.0.0",
        byte_schema = "F-B-1.0.0",
        access_schema = "F-A-1.0.0",
        claim_scope = "fixed_panel",
        executable = FALSE
    )
}
