# Assisted-by: OpenAI Codex.

if (!nzchar(Sys.getenv("TZ", unset = ""))) {
    Sys.setenv(TZ = "UTC")
}

runner_file <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(runner_file) != 1L) {
    stop("Cannot locate benchmark/run-controlled.R.", call. = FALSE)
}
benchmark_dir <- dirname(normalizePath(sub("^--file=", "", runner_file)))
repo_root <- dirname(benchmark_dir)
source(file.path(benchmark_dir, "protocol.R"))
source(file.path(benchmark_dir, "provenance.R"))
source(file.path(benchmark_dir, "fixtures.R"))
source(file.path(benchmark_dir, "controlled.R"))
source(file.path(benchmark_dir, "report.R"))

benchmark_controlled_usage <- function() {
    cat(paste0(
        "Usage: Rscript benchmark/run-controlled.R [options]\n",
        "  --output=DIR  Result directory\n",
        "  --help        Show this text\n"
    ))
}

benchmark_controlled_options <- function(args) {
    output <- NULL
    for (arg in args) {
        if (identical(arg, "--help")) {
            benchmark_controlled_usage()
            quit(save = "no", status = 0L)
        } else if (grepl("^--output=", arg)) {
            output <- sub("^--output=", "", arg)
        } else {
            stop("Unknown controlled benchmark option: ", arg, call. = FALSE)
        }
    }
    list(output = output)
}

options <- benchmark_controlled_options(commandArgs(trailingOnly = TRUE))
if (is.null(options$output)) {
    stamp <- format(Sys.time(), "%Y%m%d-%H%M%S", tz = "UTC")
    options$output <- file.path(
        benchmark_dir,
        "results",
        paste(stamp, "controlled", sep = "-")
    )
}
dir.create(options$output, recursive = TRUE, showWarnings = FALSE)
output_dir <- normalizePath(options$output, mustWork = TRUE)

packages <- c("genefunnel", "Matrix", "BiocParallel", "Rcpp", "RcppArmadillo")
missing_packages <- packages[!vapply(
    packages,
    requireNamespace,
    quietly = TRUE,
    FUN.VALUE = logical(1)
)]
if (length(missing_packages) > 0L) {
    stop(
        "Controlled benchmark dependencies are unavailable: ",
        paste(missing_packages, collapse = ", "),
        call. = FALSE
    )
}

protocol <- benchmark_read_protocol(file.path(benchmark_dir, "protocol.tsv"))
manifest <- benchmark_protocol_subset(
    protocol,
    "controlled",
    benchmark_controlled_scenarios
)
benchmark_write_tsv(protocol, file.path(output_dir, "protocol.tsv"))
benchmark_write_tsv(manifest, file.path(output_dir, "manifest.tsv"))

metadata <- benchmark_metadata(
    repo_root = repo_root,
    runner = "benchmark/run-controlled.R",
    preset = "controlled",
    repeats = 1L,
    workers = 1L
)
benchmark_write_tsv(metadata, file.path(output_dir, "metadata.tsv"))
benchmark_write_session_info(
    file.path(output_dir, "session-info.txt"),
    packages
)

results <- benchmark_controlled_results()
summary <- data.frame(
    protocol_version = GENEFUNNEL_BENCHMARK_PROTOCOL,
    total_scenarios = length(unique(results$scenario_id)),
    total_assertions = nrow(results),
    passed_assertions = sum(results$passed),
    failed_assertions = sum(!results$passed),
    all_passed = all(results$passed),
    stringsAsFactors = FALSE,
    check.names = FALSE
)
benchmark_write_tsv(results, file.path(output_dir, "results.tsv"))
benchmark_write_tsv(summary, file.path(output_dir, "summary.tsv"))
benchmark_write_controlled_report(
    file.path(output_dir, "report.md"),
    manifest,
    results,
    summary,
    metadata
)

if (!all(results$passed)) {
    failed <- paste0(
        results$scenario_id[!results$passed],
        "/",
        results$assertion_id[!results$passed]
    )
    stop(
        "Controlled benchmark assertions failed: ",
        paste(failed, collapse = ", "),
        call. = FALSE
    )
}
cat("Controlled benchmark complete: ", output_dir, "\n", sep = "")
