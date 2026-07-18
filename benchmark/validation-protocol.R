# Assisted-by: OpenAI Codex.

VALIDATION_PROTOCOL_VERSION <- "F-2.0.0"
VALIDATION_PROTOCOL_SHA256 <-
    "00a3a9f5e1cd6f709fc3ecdbad111f938b80c849c1c6071d134db41e7c63361c"

validation_protocol_sha256 <- function(path) {
    namespace <- asNamespace("tools")
    if (!exists("sha256sum", envir = namespace, inherits = FALSE)) {
        stop("Validation protocol requires tools::sha256sum().", call. = FALSE)
    }
    unname(get("sha256sum", envir = namespace, inherits = FALSE)(path))
}

validation_protocol_read <- function(root = ".") {
    path <- file.path(root, "benchmark", "validation-protocol.tsv")
    observed <- validation_protocol_sha256(path)
    if (!identical(observed, VALIDATION_PROTOCOL_SHA256)) {
        stop("Validation protocol identity is invalid.", call. = FALSE)
    }
    protocol <- utils::read.delim(
        path, stringsAsFactors = FALSE, check.names = FALSE,
        colClasses = "character", na.strings = character()
    )
    required <- c("protocol_version", "section", "key", "value")
    identity <- paste(protocol$section, protocol$key, sep = "::")
    valid <- identical(names(protocol), required) && nrow(protocol) == 184L &&
        !anyNA(protocol) && all(nzchar(unlist(protocol))) &&
        all(protocol$protocol_version == VALIDATION_PROTOCOL_VERSION) &&
        !anyDuplicated(identity) &&
        all(grepl("^[a-z][a-z0-9_]*$", protocol$section)) &&
        all(grepl("^[A-Za-z][A-Za-z0-9_]*$", protocol$key))
    if (!valid) stop("Validation protocol rows are invalid.", call. = FALSE)
    protocol
}

validation_protocol_value <- function(protocol, section, key) {
    selected <- protocol$section == section & protocol$key == key
    if (sum(selected) != 1L) {
        stop("Validation protocol key is invalid: ", section, "/", key,
            call. = FALSE)
    }
    protocol$value[[which(selected)]]
}

validation_protocol_vector <- function(protocol, section, key) {
    strsplit(
        validation_protocol_value(protocol, section, key),
        ",", fixed = TRUE
    )[[1L]]
}

validation_protocol_number <- function(protocol, section, key) {
    raw <- validation_protocol_value(protocol, section, key)
    value <- suppressWarnings(as.numeric(sub("_.*$", "", raw)))
    if (length(value) != 1L || !is.finite(value)) {
        stop("Validation protocol number is invalid: ", section, "/", key,
            call. = FALSE)
    }
    value
}

validation_protocol_factorial <- function() {
    design <- expand.grid(
        A = c(-1L, 1L), B = c(-1L, 1L), C = c(-1L, 1L),
        D = c(-1L, 1L), E = c(-1L, 1L), F = c(-1L, 1L),
        G = c(-1L, 1L), KEEP.OUT.ATTRS = FALSE
    )
    design$H <- design$A * design$B * design$C
    design$I <- design$A * design$D * design$E
    design$J <- design$B * design$D * design$F
    design$K <- design$C * design$E * design$G
    balanced <- vapply(1:3, function(degree) {
        products <- utils::combn(names(design), degree, simplify = FALSE)
        all(vapply(products, function(fields) {
            sum(Reduce(`*`, design[fields])) == 0L
        }, logical(1L)))
    }, logical(1L))
    if (nrow(design) != 128L || !all(balanced)) {
        stop("Validation fractional factorial is below resolution IV.",
            call. = FALSE)
    }
    design
}

validation_protocol_validate <- function(root = ".") {
    protocol <- validation_protocol_read(root)
    expected_sections <- c(
        identity = 9L, environment = 7L, assay = 11L, eligibility = 11L,
        split = 8L, preprocess_rna = 8L, preprocess_proteomics = 8L,
        preprocess_common = 4L, catalogue = 14L, coverage = 5L,
        methods = 14L, task = 7L, controls = 12L, primary = 14L,
        secondary = 9L, negative_control = 7L, simulation = 16L,
        seeds = 11L, output = 5L, decision = 4L
    )
    observed_sections <- table(factor(
        protocol$section, levels = names(expected_sections)
    ))
    if (!identical(as.integer(observed_sections), unname(expected_sections)) ||
        any(!protocol$section %in% names(expected_sections))) {
        stop("Validation protocol section registry is invalid.", call. = FALSE)
    }

    exact <- list(
        c("identity", "parent_index", "F-I-1.0.0"),
        c("identity", "prior_execution_protocol", "1.0.0"),
        c("environment", "r_version", "4.6.1"),
        c("environment", "bioconductor_version", "3.23"),
        c("assay", "targets", "bulk_rnaseq,pseudobulk_rnaseq,bulk_proteomics"),
        c("coverage", "primary_target", "declared_size_10_to_200;matched_fraction_at_least_0.70;matched_size_at_least_10"),
        c("methods", "ordered", "genefunnel,sum,mean,singscore,gsva,ssgsea"),
        c("primary", "families", "5_hypotheses_per_claim_stratum_times_3_strata_equals_15_co_primary_hypotheses"),
        c("environment", "bootstrap_repeats", "10000"),
        c("controls", "count", "100_matched_random_sets_per_task"),
        c("controls", "proposal_count", "100000_seeded_draws_per_task")
    )
    exact_valid <- all(vapply(exact, function(row) {
        identical(validation_protocol_value(protocol, row[[1L]], row[[2L]]),
            row[[3L]])
    }, logical(1L)))

    catalog_hashes <- vapply(c("rna_sha256", "protein_sha256", "pathway_sha256"),
        function(key) validation_protocol_value(protocol, "catalogue", key),
        character(1L))
    package_sources <- vapply(c(
        "normalization_package", "singscore_source", "gsva_source"
    ), function(key) {
        section <- if (key == "normalization_package") "preprocess_rna" else "methods"
        validation_protocol_value(protocol, section, key)
    }, character(1L))
    source_hashes <- sub(".*_SHA256_", "", package_sources)
    source_commits <- package_sources[-1L]
    source_commits <- sub(".*_commit_", "", source_commits)
    source_commits <- sub("_SHA256_.*", "", source_commits)

    lower_probability <- validation_protocol_number(
        protocol, "primary", "lower_probability"
    )
    numeric_valid <- identical(validation_protocol_vector(
        protocol, "simulation", "member_count"
    ), c("16", "64", "192")) &&
        identical(validation_protocol_vector(
            protocol, "simulation", "magnitude_fold"
        ), c("1", "1.25", "2")) &&
        identical(validation_protocol_vector(
            protocol, "simulation", "balance_change"
        ), c("-0.2", "0", "0.2_full_cross_with_magnitude_fold")) &&
        abs(lower_probability - 0.05 / 15) < 1e-16

    seed_rows <- protocol$section == "seeds"
    seeds <- suppressWarnings(as.integer(protocol$value[seed_rows]))
    seed_valid <- !anyNA(seeds) && all(seeds > 0L) && !anyDuplicated(seeds)
    archive_valid <- all(grepl("^[0-9a-f]{64}$", catalog_hashes)) &&
        all(grepl("^[0-9a-f]{64}$", source_hashes)) &&
        all(grepl("^[0-9a-f]{40}$", source_commits))

    if (!exact_valid || !numeric_valid || !seed_valid || !archive_valid) {
        stop("Validation protocol semantics are invalid.", call. = FALSE)
    }
    factorial <- validation_protocol_factorial()
    simulation_rows <- 2L * 3L * 2L * 2L * 3L * 3L * 4L *
        nrow(factorial) * 2L
    list(
        protocol_version = VALIDATION_PROTOCOL_VERSION,
        protocol_sha256 = VALIDATION_PROTOCOL_SHA256,
        registry_rows = nrow(protocol),
        claim_strata = validation_protocol_vector(
            protocol, "assay", "targets"
        ),
        methods = validation_protocol_vector(protocol, "methods", "ordered"),
        co_primary_hypotheses = 15L,
        lower_probability = lower_probability,
        simulation_rows = simulation_rows,
        executable = FALSE
    )
}
