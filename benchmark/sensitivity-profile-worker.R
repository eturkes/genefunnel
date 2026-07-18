# Assisted-by: OpenAI Codex.

if (!nzchar(Sys.getenv("TZ", unset = ""))) {
    Sys.setenv(TZ = "UTC")
}

worker_file <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(worker_file) != 1L) {
    stop("Cannot locate sensitivity-profile-worker.R.", call. = FALSE)
}
benchmark_dir <- dirname(normalizePath(sub("^--file=", "", worker_file)))
source(file.path(benchmark_dir, "measure.R"))
source(file.path(benchmark_dir, "sensitivity-profile.R"))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
    stop(
        "Usage: sensitivity-profile-worker.R CONFIG_RDS OUTPUT_DIR LIBRARY",
        call. = FALSE
    )
}
config <- readRDS(normalizePath(args[[1L]], mustWork = TRUE))
output_dir <- normalizePath(args[[2L]], mustWork = TRUE)
package_library <- normalizePath(args[[3L]], mustWork = TRUE)
.libPaths(c(package_library, .libPaths()))

for (package in c("genefunnel", "Matrix", "BiocParallel")) {
    if (!requireNamespace(package, quietly = TRUE)) {
        stop("Sensitivity profile dependency is unavailable: ", package)
    }
}
loaded_package <- normalizePath(find.package("genefunnel"), mustWork = TRUE)
expected_package <- normalizePath(
    file.path(package_library, "genefunnel"),
    mustWork = TRUE
)
if (!identical(loaded_package, expected_package)) {
    stop("Sensitivity profile loaded genefunnel from the wrong library.")
}

sensitivity <- getFromNamespace(".gene_set_sensitivity", "genefunnel")
validate_result <- getFromNamespace(".validate_sensitivity_chunk", "genefunnel")
backend <- BiocParallel::SerialParam(progressbar = FALSE)
fixture <- sensitivity_profile_fixture(config)
fixture_digest <- sensitivity_profile_digest(fixture)

profile_call <- function() {
    result <- sensitivity(fixture$mat, fixture$gene_sets, backend)
    validate_result(result, c(config$n_sets, config$n_samples))
    result
}

profile_check_digest <- function(value, expected = NULL) {
    observed <- sensitivity_profile_digest(value)
    if (!is.null(expected) && !identical(observed, expected)) {
        stop("Sensitivity profile output identity changed.", call. = FALSE)
    }
    observed
}

runs <- vector("list", config$repeats)
output_digest <- NULL
output_bytes <- NA_real_
for (repeat_id in seq_len(config$repeats)) {
    invisible(gc())
    timing <- system.time(result <- profile_call(), gcFirst = FALSE)
    digest <- profile_check_digest(result, output_digest)
    if (is.null(output_digest)) {
        output_digest <- digest
        output_bytes <- as.numeric(object.size(result))
    }
    runs[[repeat_id]] <- data.frame(
        repeat_id = repeat_id,
        elapsed_sec = unname(timing[["elapsed"]]),
        user_sec = unname(timing[["user.self"]]),
        system_sec = unname(timing[["sys.self"]]),
        output_md5 = digest,
        stringsAsFactors = FALSE,
        check.names = FALSE
    )
    rm(result)
}
runs <- do.call(rbind, runs)

rprof_path <- file.path(output_dir, "Rprof.out")
invisible(gc())
utils::Rprof(
    rprof_path,
    interval = 0.01,
    memory.profiling = FALSE,
    line.profiling = FALSE
)
on.exit(utils::Rprof(NULL), add = TRUE)
rprof_timing <- system.time(rprof_result <- profile_call(), gcFirst = FALSE)
utils::Rprof(NULL)
rprof_digest <- profile_check_digest(rprof_result, output_digest)
rm(rprof_result)
rprof <- sensitivity_profile_rprof(rprof_path)
rprof$elapsed_sec <- unname(rprof_timing[["elapsed"]])
rprof$output_md5 <- rprof_digest

allocation_path <- file.path(output_dir, "Rprofmem.out")
if (!isTRUE(capabilities("profmem"))) {
    stop("Sensitivity profile requires R allocation profiling.", call. = FALSE)
}
invisible(gc())
utils::Rprofmem(allocation_path, threshold = 0)
on.exit(utils::Rprofmem(NULL), add = TRUE)
allocation_timing <- system.time(
    allocation_result <- profile_call(),
    gcFirst = FALSE
)
utils::Rprofmem(NULL)
allocation_digest <- profile_check_digest(allocation_result, output_digest)
rm(allocation_result)
allocation <- as.list(benchmark_allocation_stats(allocation_path))
allocation$elapsed_sec <- unname(allocation_timing[["elapsed"]])
allocation$output_md5 <- allocation_digest
allocation <- as.data.frame(
    allocation,
    stringsAsFactors = FALSE,
    check.names = FALSE
)
if (any(!is.finite(runs$elapsed_sec)) || any(runs$elapsed_sec <= 0) ||
    rprof$total_samples[[1L]] < 1L ||
    !is.finite(rprof$exact_sample_share[[1L]]) ||
    allocation$manager_r_alloc_bytes[[1L]] <= 0 ||
    allocation$manager_r_alloc_events[[1L]] <= 0) {
    stop("Sensitivity profile measurement is incomplete.", call. = FALSE)
}

manifest <- data.frame(
    protocol_version = config$protocol_version,
    parent_protocol = config$parent_protocol,
    mode = config$mode,
    n_features = config$n_features,
    n_samples = config$n_samples,
    n_sets = config$n_sets,
    set_size = config$set_size,
    matrix_seed = config$matrix_seed,
    set_seed = config$set_seed,
    repeats = config$repeats,
    fixture_md5 = fixture_digest,
    matrix_bytes = as.numeric(object.size(fixture$mat)),
    gene_sets_bytes = as.numeric(object.size(fixture$gene_sets)),
    output_bytes = output_bytes,
    output_md5 = output_digest,
    stringsAsFactors = FALSE,
    check.names = FALSE
)
environment <- data.frame(
    key = c(
        "r_version", "r_platform", "genefunnel_version", "genefunnel_path",
        "Matrix_version", "Matrix_path", "BiocParallel_version",
        "BiocParallel_path"
    ),
    value = c(
        paste(R.version$major, R.version$minor, sep = "."),
        R.version$platform,
        as.character(utils::packageVersion("genefunnel")),
        loaded_package,
        as.character(utils::packageVersion("Matrix")),
        normalizePath(find.package("Matrix")),
        as.character(utils::packageVersion("BiocParallel")),
        normalizePath(find.package("BiocParallel"))
    ),
    stringsAsFactors = FALSE,
    check.names = FALSE
)

profile_write <- function(value, filename) {
    utils::write.table(
        value,
        file.path(output_dir, filename),
        sep = "\t",
        row.names = FALSE,
        col.names = TRUE,
        quote = TRUE,
        na = "NA"
    )
}
profile_write(runs, "runs.tsv")
profile_write(rprof, "profile.tsv")
profile_write(allocation, "allocation.tsv")
profile_write(manifest, "manifest.tsv")
profile_write(environment, "worker-environment.tsv")
profile_write(sensitivity_profile_decision(runs, rprof), "decision.tsv")
writeLines(capture.output(sessionInfo()), file.path(output_dir, "session-info.txt"))
cat(
    sprintf(
        "Sensitivity %s profile: median %.6f s, exact share %.6f, digest %s\n",
        config$mode,
        stats::median(runs$elapsed_sec),
        rprof$exact_sample_share,
        output_digest
    )
)
