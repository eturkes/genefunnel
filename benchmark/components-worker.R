# Assisted-by: OpenAI Codex.

if (!nzchar(Sys.getenv("TZ", unset = ""))) {
    Sys.setenv(TZ = "UTC")
}

worker_file <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(worker_file) != 1L) {
    stop("Cannot locate benchmark/components-worker.R.", call. = FALSE)
}
benchmark_dir <- dirname(normalizePath(sub("^--file=", "", worker_file)))
source(file.path(benchmark_dir, "fixtures.R"))
source(file.path(benchmark_dir, "measure.R"))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) {
    stop(
        paste(
            "Usage: components-worker.R SCENARIO_RDS RESULT_TSV",
            "ALLOCATION_LOG PACKAGE_LIBRARY"
        ),
        call. = FALSE
    )
}
scenario_path <- normalizePath(args[[1L]], mustWork = TRUE)
result_path <- args[[2L]]
allocation_path <- args[[3L]]
package_library <- normalizePath(args[[4L]], mustWork = TRUE)
.libPaths(c(package_library, .libPaths()))

for (package in c("genefunnel", "Matrix", "BiocParallel")) {
    if (!requireNamespace(package, quietly = TRUE)) {
        stop("Benchmark dependency is not installed: ", package, call. = FALSE)
    }
}
loaded_package <- normalizePath(find.package("genefunnel"), mustWork = TRUE)
expected_package <- normalizePath(
    file.path(package_library, "genefunnel"),
    mustWork = TRUE
)
if (!identical(loaded_package, expected_package)) {
    stop("Benchmark loaded genefunnel from the wrong library.", call. = FALSE)
}

scenario <- readRDS(scenario_path)
fixture <- benchmark_fixture(scenario)
fixture_stats <- benchmark_fixture_stats(fixture)
backend <- BiocParallel::SerialParam(progressbar = FALSE)

invisible(gc())
cold_timing <- system.time({
    cold_scores <- genefunnel::genefunnel(
        fixture$mat,
        fixture$gene_sets,
        BPPARAM = backend
    )
})
expected_dimensions <- c(scenario$n_sets, scenario$n_samples)
if (!identical(dim(cold_scores), expected_dimensions)) {
    stop("Benchmark output dimensions do not match the scenario.", call. = FALSE)
}
output_md5 <- benchmark_output_md5(cold_scores)
output_bytes <- as.numeric(object.size(cold_scores))
output_na_cells <- sum(is.na(cold_scores))
output_sum <- sum(cold_scores, na.rm = TRUE)

rm(cold_scores)
invisible(gc())
timed_calls <- scenario$timed_calls
if (!is.numeric(timed_calls) || length(timed_calls) != 1L ||
    !is.finite(timed_calls) || timed_calls < 1 ||
    timed_calls != floor(timed_calls)) {
    stop("Component benchmark timed-call count is invalid.", call. = FALSE)
}
timed_calls <- as.integer(timed_calls)
warm_timing <- system.time({
    for (call_index in seq_len(timed_calls)) {
        warm_scores <- genefunnel::genefunnel(
            fixture$mat,
            fixture$gene_sets,
            BPPARAM = backend
        )
    }
})
if (!identical(benchmark_output_md5(warm_scores), output_md5)) {
    stop("Cold and warm-pass outputs differ.", call. = FALSE)
}
rm(warm_scores)

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
    if (!identical(benchmark_output_md5(allocation_scores), output_md5)) {
        stop("Timed and allocation-pass outputs differ.", call. = FALSE)
    }
    rm(allocation_scores)
}
allocation_stats <- if (profiling_allocations) {
    benchmark_allocation_stats(allocation_path)
} else {
    c(manager_r_alloc_bytes = NA_real_, manager_r_alloc_events = NA_real_)
}

metrics <- list(
    cold_elapsed_sec = unname(cold_timing[["elapsed"]]),
    cold_user_sec = unname(cold_timing[["user.self"]]),
    cold_system_sec = unname(cold_timing[["sys.self"]]),
    elapsed_sec = unname(warm_timing[["elapsed"]]) / timed_calls,
    timing_batch_elapsed_sec = unname(warm_timing[["elapsed"]]),
    user_sec = unname(warm_timing[["user.self"]]) / timed_calls,
    system_sec = unname(warm_timing[["sys.self"]]) / timed_calls,
    output_bytes = output_bytes,
    output_na_cells = output_na_cells,
    output_sum = output_sum,
    output_md5 = output_md5,
    measurement_passes = "cold; five-call warm batch; allocation; identity checked",
    package_library = package_library,
    package_path = loaded_package,
    genefunnel_version = as.character(utils::packageVersion("genefunnel")),
    Matrix_version = as.character(utils::packageVersion("Matrix")),
    Matrix_path = normalizePath(find.package("Matrix"), mustWork = TRUE),
    BiocParallel_version = as.character(utils::packageVersion("BiocParallel")),
    BiocParallel_path = normalizePath(
        find.package("BiocParallel"),
        mustWork = TRUE
    ),
    Rcpp_version = as.character(utils::packageVersion("Rcpp")),
    Rcpp_path = normalizePath(find.package("Rcpp"), mustWork = TRUE),
    RcppArmadillo_version = as.character(utils::packageVersion("RcppArmadillo")),
    RcppArmadillo_path = normalizePath(
        find.package("RcppArmadillo"),
        mustWork = TRUE
    ),
    r_version = paste(R.version$major, R.version$minor, sep = ".")
)
row <- as.data.frame(
    c(
        scenario,
        as.list(fixture_stats[1L, , drop = TRUE]),
        as.list(allocation_stats),
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
cat(sprintf(
    "%s %s repeat %d: %.6f s, digest %s\n",
    scenario$scenario_id,
    scenario$method,
    scenario$repeat_id,
    metrics$elapsed_sec,
    metrics$output_md5
))
