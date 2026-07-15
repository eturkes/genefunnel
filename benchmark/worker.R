# Assisted-by: OpenAI Codex.

if (!nzchar(Sys.getenv("TZ", unset = ""))) {
    Sys.setenv(TZ = "UTC")
}

worker_file <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(worker_file) != 1L) {
    stop("Cannot locate benchmark/worker.R.", call. = FALSE)
}
benchmark_dir <- dirname(normalizePath(sub("^--file=", "", worker_file)))
source(file.path(benchmark_dir, "fixtures.R"))
source(file.path(benchmark_dir, "measure.R"))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
    stop(
        "Usage: worker.R SCENARIO_RDS RESULT_TSV ALLOCATION_LOG",
        call. = FALSE
    )
}
scenario_path <- normalizePath(args[[1L]], mustWork = TRUE)
result_path <- args[[2L]]
allocation_path <- args[[3L]]

for (package in c("genefunnel", "Matrix", "BiocParallel")) {
    if (!requireNamespace(package, quietly = TRUE)) {
        stop("Benchmark dependency is not installed: ", package, call. = FALSE)
    }
}

scenario <- readRDS(scenario_path)
fixture <- benchmark_fixture(scenario)
fixture_stats <- benchmark_fixture_stats(fixture)
backend <- switch(
    scenario$backend,
    serial = BiocParallel::SerialParam(progressbar = FALSE),
    snow = BiocParallel::SnowParam(
        workers = scenario$backend_workers,
        type = "SOCK",
        progressbar = FALSE
    ),
    stop("Unknown benchmark backend.", call. = FALSE)
)
on.exit({
    if (BiocParallel::bpisup(backend)) {
        BiocParallel::bpstop(backend)
    }
}, add = TRUE)

invisible(gc())
monitor <- benchmark_start_memory_monitor()
monitor_stopped <- FALSE
on.exit({
    if (!monitor_stopped) {
        benchmark_stop_memory_monitor(monitor)
    }
}, add = TRUE)

timing <- system.time({
    scores <- genefunnel::genefunnel(
        fixture$mat,
        fixture$gene_sets,
        BPPARAM = backend
    )
})
memory_stats <- benchmark_stop_memory_monitor(monitor)
monitor_stopped <- TRUE

profiling_allocations <- isTRUE(capabilities("profmem"))
if (profiling_allocations) {
    utils::Rprofmem(allocation_path, threshold = 0)
    on.exit(utils::Rprofmem(NULL), add = TRUE)
    allocation_scores <- genefunnel::genefunnel(
        fixture$mat,
        fixture$gene_sets,
        BPPARAM = backend
    )
    utils::Rprofmem(NULL)
    if (!identical(allocation_scores, scores)) {
        stop("Timed and allocation-pass outputs differ.", call. = FALSE)
    }
}
allocation_stats <- if (profiling_allocations) {
    benchmark_allocation_stats(allocation_path)
} else {
    c(manager_r_alloc_bytes = NA_real_, manager_r_alloc_events = NA_real_)
}

expected_dimensions <- c(scenario$n_sets, scenario$n_samples)
if (!identical(dim(scores), expected_dimensions)) {
    stop("Benchmark output dimensions do not match the scenario.", call. = FALSE)
}

elapsed <- unname(timing[["elapsed"]])
score_cells <- as.double(scenario$n_sets) * scenario$n_samples
member_cells <- score_cells * scenario$set_size
metrics <- list(
    elapsed_sec = elapsed,
    user_sec = unname(timing[["user.self"]]),
    system_sec = unname(timing[["sys.self"]]),
    score_cells = score_cells,
    score_cells_per_sec = score_cells / elapsed,
    member_cells = member_cells,
    member_cells_per_sec = member_cells / elapsed,
    output_bytes = as.numeric(object.size(scores)),
    output_na_cells = sum(is.na(scores)),
    output_sum = sum(scores, na.rm = TRUE),
    output_md5 = benchmark_output_md5(scores),
    allocation_pass = if (profiling_allocations) "separate" else "unavailable",
    memory_scope = monitor$scope,
    r_version = paste(R.version$major, R.version$minor, sep = "."),
    genefunnel_version = as.character(utils::packageVersion("genefunnel")),
    Matrix_version = as.character(utils::packageVersion("Matrix")),
    BiocParallel_version = as.character(utils::packageVersion("BiocParallel")),
    Rcpp_version = as.character(utils::packageVersion("Rcpp")),
    RcppArmadillo_version = as.character(utils::packageVersion("RcppArmadillo"))
)

row <- as.data.frame(
    c(
        scenario,
        as.list(fixture_stats[1L, , drop = TRUE]),
        as.list(allocation_stats),
        as.list(memory_stats),
        metrics
    ),
    stringsAsFactors = FALSE,
    check.names = FALSE,
    optional = TRUE
)
utils::write.table(
    row,
    result_path,
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE,
    quote = TRUE,
    na = "NA"
)
cat(
    sprintf(
        "%s repeat %d: %.6f s, digest %s\n",
        scenario$scenario_id,
        scenario$repeat_id,
        elapsed,
        metrics$output_md5
    )
)
