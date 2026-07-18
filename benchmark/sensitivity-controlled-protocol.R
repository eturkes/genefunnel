# Assisted-by: OpenAI Codex.

SENSITIVITY_CONTROLLED_VERSION <- "E-C-1.0.0"
SENSITIVITY_CONTROLLED_MD5 <- "1c119004ee749b4242f1b73ca6fd8c4c"
SENSITIVITY_CONTROLLED_RESULT_MD5 <- "107c7f8363097bdddd11bafe2d37ff83"
SENSITIVITY_CONTROLLED_CURVES_MD5 <- "0384f16f78b62101a0ead916d0a30100"

sensitivity_controlled_read_protocol <- function(root = ".") {
    path <- file.path(root, "benchmark", "sensitivity-controlled-protocol.tsv")
    if (!identical(unname(tools::md5sum(path)), SENSITIVITY_CONTROLLED_MD5)) {
        stop("Sensitivity controlled protocol identity is invalid.", call. = FALSE)
    }
    protocol <- utils::read.delim(
        path,
        stringsAsFactors = FALSE,
        check.names = FALSE,
        colClasses = "character",
        na.strings = character()
    )
    required <- c("protocol_version", "section", "key", "value")
    identity <- paste(protocol$section, protocol$key, sep = "::")
    valid <- identical(names(protocol), required) && nrow(protocol) == 45L &&
        !anyNA(protocol) && all(nzchar(unlist(protocol))) &&
        all(protocol$protocol_version == SENSITIVITY_CONTROLLED_VERSION) &&
        !anyDuplicated(identity)
    if (!valid) {
        stop("Sensitivity controlled protocol schema is invalid.", call. = FALSE)
    }
    protocol
}

sensitivity_controlled_protocol_value <- function(protocol, section, key) {
    selected <- protocol$section == section & protocol$key == key
    if (sum(selected) != 1L) {
        stop("Sensitivity controlled protocol key is invalid.", call. = FALSE)
    }
    protocol$value[[which(selected)]]
}

sensitivity_controlled_validate_protocol <- function(root = ".") {
    protocol <- sensitivity_controlled_read_protocol(root)
    parent <- sensitivity_controlled_protocol_value(
        protocol, "identity", "parent_protocol"
    )
    parent_md5 <- sensitivity_controlled_protocol_value(
        protocol, "identity", "parent_registry_md5"
    )
    observed <- unname(tools::md5sum(file.path(
        root, "benchmark", "sensitivity-protocol.tsv"
    )))
    if (!identical(parent, SENSITIVITY_PROTOCOL_VERSION) ||
        !identical(parent_md5, observed)) {
        stop("Sensitivity controlled parent identity is invalid.", call. = FALSE)
    }
    list(
        protocol_version = SENSITIVITY_CONTROLLED_VERSION,
        protocol_rows = nrow(protocol),
        protocol_md5 = SENSITIVITY_CONTROLLED_MD5,
        parent_protocol = parent,
        parent_registry_md5 = parent_md5
    )
}

sensitivity_controlled_read_result <- function(root = ".") {
    path <- file.path(root, "benchmark", "sensitivity-controlled-result.tsv")
    if (!identical(unname(tools::md5sum(path)), SENSITIVITY_CONTROLLED_RESULT_MD5)) {
        stop("Sensitivity controlled result identity is invalid.", call. = FALSE)
    }
    result <- utils::read.delim(
        path, stringsAsFactors = FALSE, check.names = FALSE,
        colClasses = "character", na.strings = character()
    )
    fields <- c(
        "parent_protocol", "execution_protocol", "parent_registry_md5",
        "execution_protocol_md5", "source_git_head", "generated_utc",
        "source_archive_sha256", "installed_manifest_md5",
        "artifact_manifest_sha256", "target", "endpoint", "estimate",
        "comparison", "threshold", "count", "passed",
        "controlled_gate_pass", "public_api_permitted"
    )
    if (!identical(names(result), fields) || nrow(result) != 4L ||
        anyNA(result) || any(!nzchar(unlist(result)))) {
        stop("Sensitivity controlled result schema is invalid.", call. = FALSE)
    }
    result
}

sensitivity_controlled_read_curves <- function(root = ".") {
    path <- file.path(root, "benchmark", "sensitivity-controlled-curves.tsv")
    if (!identical(unname(tools::md5sum(path)), SENSITIVITY_CONTROLLED_CURVES_MD5)) {
        stop("Sensitivity controlled curve identity is invalid.", call. = FALSE)
    }
    curves <- utils::read.delim(
        path, stringsAsFactors = FALSE, check.names = FALSE,
        colClasses = "character", na.strings = character()
    )
    fields <- c(
        "parent_protocol", "execution_protocol", "source_git_head",
        "removed_fraction", "mask_mechanism", "absence_mode", "rows",
        "median_feature_loss", "quantile90_feature_loss",
        "median_largest_absolute_delta_over_sum",
        "median_absolute_delta_over_sum", "median_declared_coverage",
        "median_observed_fraction", "zero_total_rows",
        "study_composition_dependent", "rescue_endpoint"
    )
    if (!identical(names(curves), fields) || nrow(curves) != 30L ||
        anyNA(curves) || any(!nzchar(unlist(curves)))) {
        stop("Sensitivity controlled curve schema is invalid.", call. = FALSE)
    }
    curves
}

sensitivity_controlled_validate_result_rows <- function(result) {
    target <- rep(c("feature", "technical"), each = 2L)
    endpoint <- rep(c(
        "median_fold_rmse_reduction", "bootstrap95_lower"
    ), times = 2L)
    estimate <- as.numeric(result$estimate)
    threshold <- as.numeric(result$threshold)
    count <- as.integer(result$count)
    valid <- identical(result$target, target) &&
        identical(result$endpoint, endpoint) &&
        identical(result$comparison, rep(">=", 4L)) &&
        identical(threshold, rep(c(0.1, 0.05), 2L)) &&
        identical(count, rep(c(10L, 2000L), 2L)) &&
        all(is.finite(estimate)) && all(estimate < threshold) &&
        all(result$passed == "FALSE") &&
        all(result$controlled_gate_pass == "FALSE") &&
        all(result$public_api_permitted == "FALSE")
    if (!valid) stop("Sensitivity controlled result rows are invalid.", call. = FALSE)
    invisible(result)
}

sensitivity_controlled_validate_curve_rows <- function(curves) {
    fractions <- rep(c("0.125", "0.25", "0.5"), each = 10L)
    mechanisms <- rep(rep(c(
        "uniform", "low_abundance", "high_abundance", "low_detection",
        "high_detection"
    ), each = 2L), times = 3L)
    modes <- rep(c("global_absence", "sample_missing"), times = 15L)
    pairs <- matrix(seq_len(nrow(curves)), nrow = 2L)
    paired_fields <- c(
        "median_feature_loss", "quantile90_feature_loss",
        "median_largest_absolute_delta_over_sum",
        "median_absolute_delta_over_sum", "zero_total_rows"
    )
    paired <- all(vapply(paired_fields, function(field) {
        identical(curves[[field]][pairs[1L, ]], curves[[field]][pairs[2L, ]])
    }, logical(1L)))
    valid <- identical(curves$removed_fraction, fractions) &&
        identical(curves$mask_mechanism, mechanisms) &&
        identical(curves$absence_mode, modes) &&
        all(as.integer(curves$rows) == 11520L) && paired &&
        all(curves$study_composition_dependent == "TRUE") &&
        all(curves$rescue_endpoint == "FALSE") &&
        all(curves$median_declared_coverage[pairs[2L, ]] == "1") &&
        all(curves$median_observed_fraction[pairs[1L, ]] == "1")
    if (!valid) stop("Sensitivity controlled curve rows are invalid.", call. = FALSE)
    invisible(curves)
}

sensitivity_controlled_validate_result <- function(root = ".") {
    result <- sensitivity_controlled_read_result(root)
    curves <- sensitivity_controlled_read_curves(root)
    shared <- list(result = result, curves = curves)
    valid <- all(vapply(shared, function(value) {
        all(value$parent_protocol == SENSITIVITY_PROTOCOL_VERSION) &&
            all(value$execution_protocol == SENSITIVITY_CONTROLLED_VERSION) &&
            length(unique(value$source_git_head)) == 1L &&
            grepl("^[0-9a-f]{40}$", value$source_git_head[[1L]])
    }, logical(1L))) &&
        identical(unique(result$source_git_head), unique(curves$source_git_head)) &&
        all(result$parent_registry_md5 == SENSITIVITY_PROTOCOL_MD5) &&
        all(result$execution_protocol_md5 == SENSITIVITY_CONTROLLED_MD5) &&
        all(grepl("^[0-9a-f]{64}$", result$source_archive_sha256)) &&
        all(grepl("^[0-9a-f]{64}$", result$artifact_manifest_sha256)) &&
        all(grepl("^[0-9a-f]{32}$", result$installed_manifest_md5))
    if (!valid) stop("Sensitivity controlled result provenance is invalid.", call. = FALSE)
    sensitivity_controlled_validate_result_rows(result)
    sensitivity_controlled_validate_curve_rows(curves)
    list(
        endpoints = nrow(result), curves = nrow(curves),
        controlled_gate_pass = FALSE, public_api_permitted = FALSE,
        source_git_head = unique(result$source_git_head)
    )
}
