# Assisted-by: OpenAI Codex.

GENEFUNNEL_BENCHMARK_PROTOCOL <- "1.0.0"

benchmark_read_protocol <- function(path) {
    protocol <- utils::read.delim(
        path,
        stringsAsFactors = FALSE,
        check.names = FALSE,
        colClasses = "character",
        na.strings = character()
    )
    required <- c(
        "protocol_version", "suite", "scenario_id", "runner",
        "matrix_seed", "set_seed", "dimensions", "matrix_construction",
        "set_construction", "preprocessing", "methods", "assertion",
        "environment_artifacts"
    )
    missing <- setdiff(required, names(protocol))
    if (length(missing) > 0L) {
        stop(
            "Benchmark protocol lacks fields: ",
            paste(missing, collapse = ", "),
            call. = FALSE
        )
    }
    empty <- vapply(protocol[required], function(column) {
        any(is.na(column) | !nzchar(column))
    }, logical(1))
    if (any(empty)) {
        stop(
            "Benchmark protocol has empty fields: ",
            paste(names(empty)[empty], collapse = ", "),
            call. = FALSE
        )
    }
    versions <- unique(protocol$protocol_version)
    if (!identical(versions, GENEFUNNEL_BENCHMARK_PROTOCOL)) {
        stop(
            "Benchmark protocol version must be ",
            GENEFUNNEL_BENCHMARK_PROTOCOL,
            ".",
            call. = FALSE
        )
    }
    keys <- paste(protocol$suite, protocol$scenario_id, sep = "::")
    if (anyDuplicated(keys)) {
        stop("Benchmark protocol suite/scenario keys must be unique.", call. = FALSE)
    }
    protocol
}

benchmark_protocol_subset <- function(protocol, suite, scenarios = NULL) {
    selected <- protocol[protocol$suite == suite, , drop = FALSE]
    if (!is.null(scenarios)) {
        missing <- setdiff(scenarios, selected$scenario_id)
        extra <- setdiff(selected$scenario_id, scenarios)
        if (length(missing) > 0L || length(extra) > 0L) {
            stop(
                "Benchmark protocol scenarios disagree with the runner: ",
                paste(c(missing, extra), collapse = ", "),
                call. = FALSE
            )
        }
        selected <- selected[match(scenarios, selected$scenario_id), , drop = FALSE]
    }
    if (nrow(selected) == 0L) {
        stop("Benchmark protocol suite is empty: ", suite, call. = FALSE)
    }
    rownames(selected) <- NULL
    selected
}

benchmark_write_tsv <- function(value, path) {
    utils::write.table(
        value,
        path,
        sep = "\t",
        row.names = FALSE,
        quote = TRUE,
        na = "NA"
    )
}
