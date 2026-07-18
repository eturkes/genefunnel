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
source(file.path(benchmark_dir, "aggregation-runner.R"), local = TRUE)
aggregation_protocol <- aggregation_validate_protocol(dirname(benchmark_dir))
stopifnot(
    identical(aggregation_protocol$registry_rows, 116L),
    identical(aggregation_protocol$data_files, 8L),
    identical(aggregation_protocol$synthetic_rows, 124416L),
    identical(aggregation_protocol$latent_scenarios, 62208)
)
stopifnot(identical(
    names(aggregation_external_file_map(".")),
    c(
        "cellbench_celseq2_counts", "cellbench_celseq2_metadata",
        "cellbench_sortseq_counts", "cellbench_sortseq_metadata",
        "kang_batch2_counts", "kang_batch2_genes", "kang_batch2_metadata",
        "reactome_v97"
    )
))

source(file.path(benchmark_dir, "aggregation-synthetic.R"), local = TRUE)
source(
    file.path(benchmark_dir, "aggregation-synthetic-summary.R"),
    local = TRUE
)
aggregation_registry <- aggregation_read_registry(file.path(
    benchmark_dir,
    "aggregation-protocol.tsv"
))
aggregation_design <- aggregation_synthetic_design(aggregation_registry)
aggregation_smoke <- aggregation_smoke_design(aggregation_design)
aggregation_audit <- getFromNamespace(".aggregation_audit", "genefunnel")
aggregation_scale <- as.numeric(aggregation_registry_value(
    aggregation_registry,
    "synthetic",
    "common_scale"
))
aggregation_first <- aggregation_simulate_design(
    aggregation_smoke,
    aggregation_audit,
    aggregation_scale,
    workers = 1L,
    chunk_size = 2L
)
aggregation_second <- aggregation_simulate_design(
    aggregation_smoke,
    aggregation_audit,
    aggregation_scale,
    workers = 1L,
    chunk_size = 3L
)
aggregation_validate_smoke(aggregation_first, aggregation_smoke)
stopifnot(identical(aggregation_first, aggregation_second))
aggregation_validate_summary_smoke(
    aggregation_design,
    aggregation_registry
)

source(file.path(benchmark_dir, "aggregation-cellbench.R"), local = TRUE)
aggregation_validate_cellbench_smoke(
    aggregation_registry,
    aggregation_audit,
    getExportedValue("genefunnel", "genefunnel")
)

source(file.path(benchmark_dir, "aggregation-kang.R"), local = TRUE)
source(file.path(benchmark_dir, "aggregation-kang-summary.R"), local = TRUE)
aggregation_validate_kang_smoke(aggregation_registry, aggregation_audit)

source(file.path(benchmark_dir, "sensitivity-protocol.R"), local = TRUE)
sensitivity_protocol <- sensitivity_validate_protocol(dirname(benchmark_dir))
stopifnot(
    identical(sensitivity_protocol$registry_rows, 90L),
    identical(sensitivity_protocol$scenarios, 5760L),
    identical(sensitivity_protocol$feature_rows, 345600L),
    identical(sensitivity_protocol$technical_rows, 5760L),
    identical(sensitivity_protocol$scenarios_per_fold, 576L)
)
source(file.path(benchmark_dir, "sensitivity-profile.R"), local = TRUE)
sensitivity_profile_protocol <- sensitivity_profile_validate(
    dirname(benchmark_dir)
)
sensitivity_registry <- sensitivity_read_registry(file.path(
    benchmark_dir,
    "sensitivity-protocol.tsv"
))
sensitivity_smoke_config <- sensitivity_profile_config(
    sensitivity_registry,
    "smoke"
)
sensitivity_smoke_fixture <- sensitivity_profile_fixture(
    sensitivity_smoke_config
)
stopifnot(
    identical(sensitivity_profile_protocol$protocol_rows, 26L),
    identical(sensitivity_profile_protocol$parent_protocol, "E-1.0.0"),
    identical(dim(sensitivity_smoke_fixture$mat), c(256L, 3L)),
    identical(length(sensitivity_smoke_fixture$gene_sets), 4L),
    identical(
        sensitivity_profile_digest(sensitivity_smoke_fixture),
        sensitivity_profile_digest(sensitivity_profile_fixture(
            sensitivity_smoke_config
        ))
    )
)
sensitivity_profile_result <- sensitivity_profile_validate_result(
    dirname(benchmark_dir)
)
stopifnot(
    identical(sensitivity_profile_result$optimization_eligible, "TRUE"),
    identical(sensitivity_profile_result$performance_claim, "FALSE")
)
sensitivity_optimization_result <- sensitivity_profile_validate_optimization(
    dirname(benchmark_dir)
)
stopifnot(
    identical(sensitivity_optimization_result$fixed_workload_identity, "TRUE"),
    identical(sensitivity_optimization_result$performance_claim, "FALSE")
)
source(
    file.path(benchmark_dir, "sensitivity-controlled-protocol.R"),
    local = TRUE
)
sensitivity_controlled_protocol <- sensitivity_controlled_validate_protocol(
    dirname(benchmark_dir)
)
stopifnot(
    identical(sensitivity_controlled_protocol$protocol_version, "E-C-1.0.0"),
    identical(sensitivity_controlled_protocol$protocol_rows, 45L),
    identical(sensitivity_controlled_protocol$parent_protocol, "E-1.0.0")
)

cat(
    "Benchmark, aggregation, and sensitivity protocol smokes passed: ",
    normalizePath(output_root),
    "\n",
    sep = ""
)
