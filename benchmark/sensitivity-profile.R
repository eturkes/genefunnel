# Assisted-by: OpenAI Codex.

SENSITIVITY_PROFILE_VERSION <- "E-P-1.0.0"
SENSITIVITY_PROFILE_MD5 <- "ff86e032b7a766be20c5d61d940991ab"

sensitivity_profile_read <- function(root = ".") {
    path <- file.path(root, "benchmark", "sensitivity-profile-protocol.tsv")
    if (!identical(unname(tools::md5sum(path)), SENSITIVITY_PROFILE_MD5)) {
        stop("Sensitivity profile protocol identity is invalid.", call. = FALSE)
    }
    protocol <- utils::read.delim(
        path,
        stringsAsFactors = FALSE,
        check.names = FALSE,
        colClasses = "character"
    )
    required <- c("protocol_version", "section", "key", "value")
    identity <- paste(protocol$section, protocol$key, sep = "::")
    if (!identical(names(protocol), required) || nrow(protocol) != 26L ||
        anyNA(protocol) || any(!nzchar(unlist(protocol))) ||
        any(protocol$protocol_version != SENSITIVITY_PROFILE_VERSION) ||
        anyDuplicated(identity)) {
        stop("Sensitivity profile protocol schema is invalid.", call. = FALSE)
    }
    protocol
}

sensitivity_profile_value <- function(protocol, section, key) {
    selected <- protocol$section == section & protocol$key == key
    if (sum(selected) != 1L) {
        stop("Sensitivity profile protocol key is invalid.", call. = FALSE)
    }
    protocol$value[[which(selected)]]
}

sensitivity_profile_validate <- function(root = ".") {
    protocol <- sensitivity_profile_read(root)
    parent <- sensitivity_profile_value(protocol, "identity", "parent_protocol")
    parent_md5 <- sensitivity_profile_value(
        protocol,
        "identity",
        "parent_registry_md5"
    )
    observed_parent <- unname(tools::md5sum(file.path(
        root,
        "benchmark",
        "sensitivity-protocol.tsv"
    )))
    if (!identical(parent, "E-1.0.0") ||
        !identical(parent_md5, observed_parent)) {
        stop("Sensitivity profile parent identity is invalid.", call. = FALSE)
    }
    list(
        protocol_version = SENSITIVITY_PROFILE_VERSION,
        protocol_rows = nrow(protocol),
        protocol_md5 = SENSITIVITY_PROFILE_MD5,
        parent_protocol = parent,
        parent_registry_md5 = parent_md5
    )
}

sensitivity_profile_config <- function(registry, mode = "gate") {
    if (!mode %in% c("smoke", "gate")) {
        stop("Sensitivity profile mode is invalid.", call. = FALSE)
    }
    integer_value <- function(key) {
        as.integer(sensitivity_registry_value(registry, "profile", key))
    }
    config <- list(
        protocol_version = SENSITIVITY_PROFILE_VERSION,
        parent_protocol = "E-1.0.0",
        mode = mode,
        n_features = integer_value("n_features"),
        n_samples = integer_value("n_samples"),
        n_sets = integer_value("n_sets"),
        set_size = integer_value("set_size"),
        matrix_seed = integer_value("matrix_seed"),
        set_seed = integer_value("set_seed"),
        repeats = integer_value("repeats")
    )
    if (mode == "smoke") {
        config[c("n_features", "n_samples", "n_sets", "set_size", "repeats")] <-
            list(256L, 3L, 4L, 8L, 1L)
    }
    config
}

sensitivity_profile_fixture <- function(config) {
    RNGkind("Mersenne-Twister", "Inversion", "Rejection")
    features <- sprintf("feature_%04d", seq_len(config$n_features))
    samples <- sprintf("sample_%02d", seq_len(config$n_samples))
    set.seed(config$matrix_seed)
    mat <- matrix(
        stats::rexp(as.double(config$n_features) * config$n_samples, rate = 1),
        nrow = config$n_features,
        ncol = config$n_samples,
        dimnames = list(features, samples)
    )
    set.seed(config$set_seed)
    gene_sets <- lapply(seq_len(config$n_sets), function(index) {
        sample(features, config$set_size, replace = FALSE)
    })
    names(gene_sets) <- sprintf("set_%02d", seq_len(config$n_sets))
    if (anyDuplicated(features) || any(lengths(gene_sets) != config$set_size) ||
        any(vapply(gene_sets, anyDuplicated, integer(1L)) > 0L)) {
        stop("Sensitivity profile fixture is invalid.", call. = FALSE)
    }
    list(mat = mat, gene_sets = gene_sets)
}

sensitivity_profile_digest <- function(value) {
    path <- tempfile("genefunnel-sensitivity-profile-", fileext = ".rds")
    on.exit(unlink(path), add = TRUE)
    saveRDS(value, path, compress = FALSE, version = 3L)
    unname(tools::md5sum(path))
}

sensitivity_profile_rprof <- function(path) {
    lines <- readLines(path, warn = FALSE)
    records <- lines[nzchar(lines) & !startsWith(lines, "sample.interval=")]
    exact <- grepl("\\.gf_", records)
    total <- length(records)
    exact_count <- sum(exact)
    data.frame(
        sample_interval_sec = 0.01,
        total_samples = as.integer(total),
        exact_samples = as.integer(exact_count),
        exact_sample_share = if (total == 0L) NA_real_ else exact_count / total,
        stringsAsFactors = FALSE,
        check.names = FALSE
    )
}

sensitivity_profile_decision <- function(runs, profile) {
    median_elapsed <- stats::median(runs$elapsed_sec)
    exact_share <- profile$exact_sample_share[[1L]]
    elapsed_trigger <- is.finite(median_elapsed) && median_elapsed > 60
    exact_trigger <- is.finite(exact_share) && exact_share > 0.5
    data.frame(
        protocol_version = SENSITIVITY_PROFILE_VERSION,
        median_elapsed_sec = median_elapsed,
        exact_sample_share = exact_share,
        elapsed_trigger = elapsed_trigger,
        exact_trigger = exact_trigger,
        optimization_eligible = elapsed_trigger || exact_trigger,
        performance_claim = FALSE,
        stringsAsFactors = FALSE,
        check.names = FALSE
    )
}
