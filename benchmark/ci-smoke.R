# Assisted-by: OpenAI Codex.

smoke_file <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(smoke_file) != 1L) {
    stop("Cannot locate benchmark/ci-smoke.R.", call. = FALSE)
}
benchmark_dir <- dirname(normalizePath(sub("^--file=", "", smoke_file)))
output <- Sys.getenv(
    "GENEFUNNEL_BENCHMARK_OUTPUT",
    file.path(benchmark_dir, "results", "ci-smoke")
)
unlink(output, recursive = TRUE, force = TRUE)

status <- system2(
    Sys.which("Rscript"),
    c(
        "--vanilla",
        shQuote(file.path(benchmark_dir, "run.R")),
        "--preset=smoke",
        "--repeats=1",
        "--workers=2",
        shQuote(paste0("--output=", output))
    )
)
if (!identical(status, 0L)) {
    stop("Benchmark smoke runner failed.", call. = FALSE)
}

summary <- utils::read.delim(
    file.path(output, "summary.tsv"),
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
        output,
        c("manifest.tsv", "metadata.tsv", "runs.tsv", "session-info.txt")
    )))
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
cat("Benchmark smoke checks passed: ", normalizePath(output), "\n", sep = "")
