# Assisted-by: OpenAI Codex.

smoke_file <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(smoke_file) != 1L) {
    stop("Cannot locate benchmark/ci-smoke.R.", call. = FALSE)
}
benchmark_dir <- dirname(normalizePath(sub("^--file=", "", smoke_file)))
output_root <- Sys.getenv(
    "GENEFUNNEL_BENCHMARK_OUTPUT",
    file.path(benchmark_dir, "results", "ci-smoke")
)
unlink(output_root, recursive = TRUE, force = TRUE)
performance_output <- file.path(output_root, "performance")
controlled_output <- file.path(output_root, "controlled")

status <- system2(
    Sys.which("Rscript"),
    c(
        "--vanilla",
        shQuote(file.path(benchmark_dir, "run.R")),
        "--preset=smoke",
        "--repeats=1",
        "--workers=2",
        shQuote(paste0("--output=", performance_output))
    )
)
if (!identical(status, 0L)) {
    stop("Benchmark smoke runner failed.", call. = FALSE)
}

summary <- utils::read.delim(
    file.path(performance_output, "summary.tsv"),
    stringsAsFactors = FALSE,
    check.names = FALSE,
    na.strings = "NA"
)
manifest <- utils::read.delim(
    file.path(performance_output, "manifest.tsv"),
    stringsAsFactors = FALSE,
    check.names = FALSE,
    na.strings = "NA"
)
stopifnot(
    nrow(summary) == 8L,
    all(summary$repeats == 1L),
    all(summary$elapsed_median_sec > 0),
    all(summary$score_cells_per_sec_median > 0),
    all(nzchar(summary$output_md5)),
    all(file.exists(file.path(
        performance_output,
        c(
            "protocol.tsv", "manifest.tsv", "metadata.tsv", "runs.tsv",
            "summary.tsv", "session-info.txt", "report.md"
        )
    ))),
    identical(unique(manifest$protocol_version), "1.0.0")
)
sparse <- summary$storage == "sparse"
fixture_digests <- split(summary$output_md5, summary$fixture_id, drop = TRUE)
stopifnot(
    any(sparse),
    all(summary$input_bytes[sparse] < summary$logical_dense_bytes[sparse]),
    all(vapply(
        fixture_digests,
        function(digests) length(unique(digests)) == 1L,
        logical(1)
    ))
)

controlled_status <- system2(
    Sys.which("Rscript"),
    c(
        "--vanilla",
        shQuote(file.path(benchmark_dir, "run-controlled.R")),
        shQuote(paste0("--output=", controlled_output))
    )
)
if (!identical(controlled_status, 0L)) {
    stop("Controlled benchmark runner failed.", call. = FALSE)
}
controlled_summary <- utils::read.delim(
    file.path(controlled_output, "summary.tsv"),
    stringsAsFactors = FALSE,
    check.names = FALSE,
    na.strings = "NA"
)
controlled_results <- utils::read.delim(
    file.path(controlled_output, "results.tsv"),
    stringsAsFactors = FALSE,
    check.names = FALSE,
    na.strings = "NA"
)
controlled_manifest <- utils::read.delim(
    file.path(controlled_output, "manifest.tsv"),
    stringsAsFactors = FALSE,
    check.names = FALSE,
    na.strings = "NA"
)
stopifnot(
    nrow(controlled_summary) == 1L,
    identical(controlled_summary$protocol_version, "1.0.0"),
    isTRUE(controlled_summary$all_passed),
    controlled_summary$failed_assertions == 0L,
    controlled_summary$total_scenarios == nrow(controlled_manifest),
    all(controlled_results$passed),
    setequal(controlled_results$scenario_id, controlled_manifest$scenario_id),
    all(file.exists(file.path(
        controlled_output,
        c(
            "protocol.tsv", "manifest.tsv", "metadata.tsv", "results.tsv",
            "summary.tsv", "session-info.txt", "report.md"
        )
    )))
)
controlled_report <- readLines(
    file.path(controlled_output, "report.md"),
    warn = FALSE
)
stopifnot(any(grepl("no claim that a thesis analysis was reproduced", controlled_report)))

source(file.path(benchmark_dir, "aggregation-protocol.R"), local = TRUE)
aggregation_protocol <- aggregation_validate_protocol(dirname(benchmark_dir))
stopifnot(
    identical(aggregation_protocol$registry_rows, 60L),
    identical(aggregation_protocol$data_files, 8L),
    identical(aggregation_protocol$synthetic_rows, 124416L),
    identical(aggregation_protocol$latent_scenarios, 62208)
)

cat(
    "Benchmark and aggregation-protocol smoke checks passed: ",
    normalizePath(output_root),
    "\n",
    sep = ""
)
