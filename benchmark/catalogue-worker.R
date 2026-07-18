# Assisted-by: OpenAI Codex.

if (!nzchar(Sys.getenv("TZ", unset = ""))) {
    Sys.setenv(TZ = "UTC")
}

worker_file <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(worker_file) != 1L) {
    stop("Cannot locate benchmark/catalogue-worker.R.", call. = FALSE)
}
benchmark_dir <- dirname(normalizePath(sub("^--file=", "", worker_file)))
source(file.path(benchmark_dir, "fixtures.R"))
source(file.path(benchmark_dir, "measure.R"))
source(file.path(benchmark_dir, "catalogue-fixtures.R"))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 6L) {
    stop(paste(
        "Usage: catalogue-worker.R SCENARIO_RDS RESULT_TSV",
        "ALLOCATION_PREFIX PACKAGE_LIBRARY METHOD PHASE"
    ), call. = FALSE)
}
scenario_path <- normalizePath(args[[1L]], mustWork = TRUE)
result_path <- args[[2L]]
allocation_prefix <- args[[3L]]
package_library <- normalizePath(args[[4L]], mustWork = TRUE)
method <- args[[5L]]
phase <- args[[6L]]
if (!method %in% c("list", "compiled")) {
    stop("Catalogue worker method must be list or compiled.", call. = FALSE)
}
if (!phase %in% c("timing", "allocation", "rss")) {
    stop("Catalogue worker phase must be timing, allocation, or rss.", call. = FALSE)
}
.libPaths(c(package_library, .libPaths()))
for (package in c(
    "genefunnel", "Matrix", "BiocParallel", "Rcpp", "RcppArmadillo"
)) {
    if (!requireNamespace(package, quietly = TRUE)) {
        stop("Catalogue benchmark dependency is unavailable: ", package, call. = FALSE)
    }
}
loaded_package <- normalizePath(find.package("genefunnel"), mustWork = TRUE)
expected_package <- normalizePath(
    file.path(package_library, "genefunnel"),
    mustWork = TRUE
)
if (!identical(loaded_package, expected_package)) {
    stop("Catalogue worker loaded genefunnel from the wrong library.", call. = FALSE)
}

scenario <- readRDS(scenario_path)
fixture <- catalogue_fixture(scenario)
backend <- BiocParallel::SerialParam(progressbar = FALSE)
compile_catalogue <- if (method == "compiled") {
    getFromNamespace(".compile_gene_sets", "genefunnel")
} else {
    NULL
}
score_compiled <- if (method == "compiled") {
    getFromNamespace(".genefunnel_compiled", "genefunnel")
} else {
    NULL
}

elapsed_now <- function() unname(proc.time()[["elapsed"]])

catalogue_score_batch <- function(batch, compiled = NULL) {
    warnings <- character()
    score <- withCallingHandlers(
        if (method == "compiled") {
            score_compiled(fixture$batches[[batch]], compiled, backend)
        } else {
            genefunnel::genefunnel(
                fixture$batches[[batch]],
                fixture$gene_sets,
                BPPARAM = backend
            )
        },
        warning = function(condition) {
            warnings <<- c(warnings, conditionMessage(condition))
            invokeRestart("muffleWarning")
        }
    )
    if (length(warnings) > 0L) {
        stop(
            "Catalogue fixture emitted unexpected warnings: ",
            paste(unique(warnings), collapse = " | "),
            call. = FALSE
        )
    }
    score
}

catalogue_score_batches <- function(compiled = NULL) {
    lapply(
        seq_len(scenario$batches),
        catalogue_score_batch,
        compiled = compiled
    )
}

catalogue_output_facts <- function(scores) {
    list(
        output_sha256 = catalogue_score_digest(scores),
        output_bytes = sum(vapply(scores, function(score) {
            as.numeric(object.size(score))
        }, numeric(1))),
        output_na_cells = sum(vapply(scores, function(score) {
            sum(is.na(score))
        }, integer(1))),
        output_sum = sum(vapply(scores, sum, numeric(1), na.rm = TRUE))
    )
}

catalogue_timing_phase <- function() {
    origin <- elapsed_now()
    compiled <- NULL
    compile_elapsed <- 0
    if (method == "compiled") {
        before_compile <- elapsed_now()
        compiled <- compile_catalogue(
            fixture$gene_sets,
            rownames(fixture$batches[[1L]])
        )
        compile_elapsed <- elapsed_now() - before_compile
    }
    scores <- vector("list", scenario$batches)
    call_elapsed <- numeric(scenario$batches)
    cumulative_elapsed <- numeric(scenario$batches)
    for (batch in seq_len(scenario$batches)) {
        before_call <- elapsed_now()
        scores[[batch]] <- catalogue_score_batch(batch, compiled)
        call_elapsed[[batch]] <- elapsed_now() - before_call
        cumulative_elapsed[[batch]] <- elapsed_now() - origin
    }
    call_fields <- setNames(
        as.list(rep(NA_real_, 5L)),
        paste0("call_", seq_len(5L), "_sec")
    )
    cumulative_fields <- setNames(
        as.list(rep(NA_real_, 5L)),
        paste0("cumulative_", seq_len(5L), "_sec")
    )
    call_fields[seq_len(scenario$batches)] <- as.list(call_elapsed)
    cumulative_fields[seq_len(scenario$batches)] <- as.list(cumulative_elapsed)
    list(
        scores = scores,
        metrics = c(
            list(
                compile_elapsed_sec = compile_elapsed,
                total_elapsed_sec = cumulative_elapsed[[scenario$batches]]
            ),
            call_fields,
            cumulative_fields
        )
    )
}

catalogue_profile <- function(path, expression) {
    if (!isTRUE(capabilities("profmem"))) {
        value <- force(expression)
        return(list(
            value = value,
            stats = c(
                manager_r_alloc_bytes = NA_real_,
                manager_r_alloc_events = NA_real_
            )
        ))
    }
    utils::Rprofmem(path, threshold = 0)
    on.exit(utils::Rprofmem(NULL), add = TRUE)
    value <- force(expression)
    utils::Rprofmem(NULL)
    list(value = value, stats = benchmark_allocation_stats(path))
}

catalogue_allocation_phase <- function() {
    operation <- catalogue_profile(
        paste0(allocation_prefix, "-operation.log"),
        {
            compiled <- if (method == "compiled") {
                compile_catalogue(
                    fixture$gene_sets,
                    rownames(fixture$batches[[1L]])
                )
            } else {
                NULL
            }
            list(compiled = compiled, scores = catalogue_score_batches(compiled))
        }
    )
    prepared <- if (method == "compiled") {
        compile_catalogue(
            fixture$gene_sets,
            rownames(fixture$batches[[1L]])
        )
    } else {
        NULL
    }
    scoring <- catalogue_profile(
        paste0(allocation_prefix, "-scoring.log"),
        catalogue_score_batches(prepared)
    )
    if (!identical(
        catalogue_score_digest(scoring$value),
        catalogue_score_digest(operation$value$scores)
    )) {
        stop("Catalogue allocation-pass outputs differ.", call. = FALSE)
    }
    compile_stats <- c(
        manager_r_alloc_bytes = NA_real_,
        manager_r_alloc_events = NA_real_
    )
    object_bytes <- wire_bytes <- validation_sec <- canonical_encoding_sec <-
        serialization_sec <- deserialization_sec <- NA_real_
    catalogue_sha256 <- feature_sha256 <- content_sha256 <- wire_sha256 <-
        NA_character_
    if (method == "compiled") {
        compilation <- catalogue_profile(
            paste0(allocation_prefix, "-compile.log"),
            compile_catalogue(
                fixture$gene_sets,
                rownames(fixture$batches[[1L]])
            )
        )
        compile_stats <- compilation$stats
        compiled <- operation$value$compiled
        before_validation <- elapsed_now()
        getFromNamespace(".validate_compiled_catalogue", "genefunnel")(compiled)
        validation_sec <- elapsed_now() - before_validation
        state <- unclass(compiled)
        before_encoding <- elapsed_now()
        canonical_bytes <- getFromNamespace(
            ".catalogue_content_body",
            "genefunnel"
        )(state)
        canonical_encoding_sec <- elapsed_now() - before_encoding
        if (!identical(
            unname(tools::sha256sum(bytes = canonical_bytes)),
            state$fingerprints[["content"]]
        )) {
            stop("Catalogue canonical encoding fingerprint differs.", call. = FALSE)
        }
        before_serialization <- elapsed_now()
        wire <- getFromNamespace(
            ".serialize_compiled_catalogue",
            "genefunnel"
        )(compiled)
        serialization_sec <- elapsed_now() - before_serialization
        before_deserialization <- elapsed_now()
        restored <- getFromNamespace(
            ".unserialize_compiled_catalogue",
            "genefunnel"
        )(wire)
        deserialization_sec <- elapsed_now() - before_deserialization
        if (!identical(restored, compiled)) {
            stop("Catalogue wire round trip changed the object.", call. = FALSE)
        }
        object_bytes <- as.numeric(object.size(compiled))
        wire_bytes <- length(wire)
        catalogue_sha256 <- state$fingerprints[["catalogue"]]
        feature_sha256 <- state$fingerprints[["features"]]
        content_sha256 <- state$fingerprints[["content"]]
        wire_sha256 <- unname(tools::sha256sum(bytes = wire))
    }
    canonical_memberships <- scenario$sets * scenario$canonical_set_size
    list(
        scores = operation$value$scores,
        metrics = list(
            validation_sec = validation_sec,
            canonical_encoding_sec = canonical_encoding_sec,
            serialization_sec = serialization_sec,
            deserialization_sec = deserialization_sec,
            operation_alloc_bytes = operation$stats[["manager_r_alloc_bytes"]],
            scoring_alloc_bytes = scoring$stats[["manager_r_alloc_bytes"]],
            compile_alloc_bytes = compile_stats[["manager_r_alloc_bytes"]],
            canonical_memberships = canonical_memberships,
            object_bytes = object_bytes,
            object_bytes_per_membership = object_bytes / canonical_memberships,
            wire_bytes = wire_bytes,
            wire_bytes_per_membership = wire_bytes / canonical_memberships,
            catalogue_sha256 = catalogue_sha256,
            feature_sha256 = feature_sha256,
            content_sha256 = content_sha256,
            wire_sha256 = wire_sha256
        )
    )
}

catalogue_monitored <- function(label, expression) {
    invisible(gc())
    monitor <- benchmark_start_memory_monitor()
    stopped <- FALSE
    on.exit({
        if (!stopped) {
            benchmark_stop_memory_monitor(monitor)
        }
    }, add = TRUE)
    value <- force(expression)
    stats <- benchmark_stop_memory_monitor(monitor)
    stopped <- TRUE
    names(stats) <- sub("^score_", paste0(label, "_"), names(stats))
    names(stats)[names(stats) == "memory_samples"] <- paste0(
        label,
        "_memory_samples"
    )
    list(
        value = value,
        stats = stats,
        scope = setNames(monitor$scope, paste0(label, "_memory_scope"))
    )
}

catalogue_rss_phase <- function() {
    if (method == "compiled") {
        compilation <- catalogue_monitored(
            "compile",
            compile_catalogue(
                fixture$gene_sets,
                rownames(fixture$batches[[1L]])
            )
        )
        compiled <- compilation$value
        compile_stats <- c(as.list(compilation$stats), compilation$scope)
    } else {
        compiled <- NULL
        compile_stats <- list(
            compile_baseline_rss_kib = NA_real_,
            compile_peak_tree_rss_kib = NA_real_,
            compile_peak_increment_kib = NA_real_,
            compile_peak_processes = NA_real_,
            compile_memory_samples = NA_real_,
            compile_memory_scope = "not applicable"
        )
    }
    scoring <- catalogue_monitored(
        "scoring",
        catalogue_score_batches(compiled)
    )
    list(
        scores = scoring$value,
        metrics = c(compile_stats, as.list(scoring$stats), scoring$scope)
    )
}

phase_result <- switch(
    phase,
    timing = catalogue_timing_phase(),
    allocation = catalogue_allocation_phase(),
    rss = catalogue_rss_phase()
)
facts <- catalogue_output_facts(phase_result$scores)
environment <- list(
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
    RcppArmadillo_version = as.character(
        utils::packageVersion("RcppArmadillo")
    ),
    RcppArmadillo_path = normalizePath(
        find.package("RcppArmadillo"),
        mustWork = TRUE
    ),
    R_version = paste(R.version$major, R.version$minor, sep = "."),
    R_platform = R.version$platform
)
row <- as.data.frame(
    c(scenario, list(phase = phase), facts, phase_result$metrics, environment),
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
    "%s %s repeat %d %s: digest %s\n",
    scenario$scenario_id,
    method,
    scenario$repeat_id,
    phase,
    facts$output_sha256
))
