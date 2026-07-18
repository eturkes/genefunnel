# Assisted-by: OpenAI Codex.

BENCHMARK_PROTOCOL_INDEX_VERSION <- "F-I-1.0.0"
BENCHMARK_PROTOCOL_INDEX_SHA256 <-
    "3d34bf1bf40cadc0f64fdbc65fbe0e7ac6fff3e888a6cb74a6c069fc2d67577b"

benchmark_protocol_sha256 <- function(paths) {
    namespace <- asNamespace("tools")
    if (!exists("sha256sum", envir = namespace, inherits = FALSE)) {
        stop("The protocol harness requires tools::sha256sum().", call. = FALSE)
    }
    unname(get("sha256sum", envir = namespace, inherits = FALSE)(paths))
}

benchmark_protocol_read_index <- function(root = ".") {
    path <- file.path(root, "benchmark", "protocol-index.tsv")
    observed <- benchmark_protocol_sha256(path)
    if (!identical(observed, BENCHMARK_PROTOCOL_INDEX_SHA256)) {
        stop("Benchmark protocol index identity is invalid.", call. = FALSE)
    }
    utils::read.delim(
        path, stringsAsFactors = FALSE, check.names = FALSE,
        colClasses = "character", na.strings = character()
    )
}

benchmark_protocol_validate_rows <- function(index) {
    fields <- c(
        "index_version", "protocol_version", "scope", "role", "order_code",
        "path", "sha256"
    )
    scalar_valid <- identical(names(index), fields) && nrow(index) > 0L &&
        !anyNA(index) && all(nzchar(unlist(index))) &&
        all(index$index_version == BENCHMARK_PROTOCOL_INDEX_VERSION) &&
        all(grepl("^[0-9]+[.][0-9]+[.][0-9]+$", index$protocol_version)) &&
        all(grepl("^(shared|[a-z][a-z0-9_]*)$", index$scope)) &&
        all(grepl("^[a-z][a-z0-9_]*$", index$role)) &&
        identical(index$order_code, sprintf("%03d", seq_len(nrow(index)))) &&
        all(grepl("^[0-9a-f]{64}$", index$sha256))
    safe_paths <- all(startsWith(index$path, "benchmark/")) &&
        !any(grepl("(^|/)[.][.](/|$)|^[A-Za-z]:|\\\\", index$path))
    keys <- paste(index$protocol_version, index$scope, index$role, sep = "::")
    file_keys <- paste(index$protocol_version, index$path, sep = "::")
    if (!scalar_valid || !safe_paths || anyDuplicated(keys) ||
        anyDuplicated(file_keys)) {
        stop("Benchmark protocol index rows are invalid.", call. = FALSE)
    }
    invisible(index)
}

benchmark_protocol_resolve_paths <- function(index, root) {
    root <- normalizePath(root, mustWork = TRUE)
    paths <- file.path(root, index$path)
    resolved <- normalizePath(paths, mustWork = TRUE)
    prefix <- paste0(root, .Platform$file.sep)
    if (!all(startsWith(resolved, prefix))) {
        stop("Benchmark protocol file escapes the repository.", call. = FALSE)
    }
    resolved
}

benchmark_protocol_validate_manifests <- function(index, paths) {
    versions <- unique(index$protocol_version)
    for (version in versions) {
        rows <- index$protocol_version == version &
            index$scope == "shared" & index$role == "protocol_manifest"
        if (sum(rows) != 1L) {
            stop("Protocol version lacks one manifest: ", version, call. = FALSE)
        }
        manifest <- utils::read.delim(
            paths[rows], stringsAsFactors = FALSE, check.names = FALSE,
            colClasses = "character", na.strings = character()
        )
        if (!"protocol_version" %in% names(manifest) ||
            !identical(unique(manifest$protocol_version), version)) {
            stop("Protocol manifest version is invalid: ", version, call. = FALSE)
        }
    }
    invisible(index)
}

benchmark_protocol_validate_runners <- function(index) {
    available <- unique(index[index$scope != "shared", c(
        "protocol_version", "scope"
    ), drop = FALSE])
    valid <- nrow(available) > 0L && all(vapply(seq_len(nrow(available)), function(i) {
        rows <- index$protocol_version == available$protocol_version[[i]] &
            index$scope == available$scope[[i]] & index$role == "runner"
        sum(rows) == 1L
    }, logical(1L)))
    if (!valid) stop("Protocol suite runner mapping is invalid.", call. = FALSE)
    available
}

benchmark_protocol_validate_index <- function(root = ".") {
    index <- benchmark_protocol_read_index(root)
    benchmark_protocol_validate_rows(index)
    paths <- benchmark_protocol_resolve_paths(index, root)
    observed <- benchmark_protocol_sha256(paths)
    invalid <- is.na(observed) | observed != index$sha256
    if (any(invalid)) {
        stop(
            "Frozen protocol file identity is invalid: ",
            paste(index$path[invalid], collapse = ", "), call. = FALSE
        )
    }
    benchmark_protocol_validate_manifests(index, paths)
    available <- benchmark_protocol_validate_runners(index)
    list(
        index = index, index_version = BENCHMARK_PROTOCOL_INDEX_VERSION,
        protocol_versions = unique(index$protocol_version),
        suites = available, files = nrow(index)
    )
}

benchmark_protocol_select <- function(index, version, suite) {
    if (length(version) != 1L || !nzchar(version) ||
        length(suite) != 1L || !nzchar(suite) || identical(suite, "shared")) {
        stop("Protocol version and executable suite must be scalar.", call. = FALSE)
    }
    available <- index$protocol_version == version & index$scope == suite
    if (!any(available)) {
        stop("Unknown benchmark protocol/suite: ", version, "/", suite,
            call. = FALSE)
    }
    rows <- index$protocol_version == version &
        index$scope %in% c("shared", suite)
    selected <- index[rows, , drop = FALSE]
    selected[order(selected$order_code), , drop = FALSE]
}

benchmark_protocol_runner <- function(selected, root = ".") {
    rows <- selected$scope != "shared" & selected$role == "runner"
    if (sum(rows) != 1L) {
        stop("Selected protocol has no unique runner.", call. = FALSE)
    }
    benchmark_protocol_resolve_paths(selected[rows, , drop = FALSE], root)[[1L]]
}
