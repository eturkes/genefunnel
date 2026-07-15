# Assisted-by: OpenAI Codex.

if (!nzchar(Sys.getenv("TZ", unset = ""))) {
    Sys.setenv(TZ = "UTC")
}

runner_file <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(runner_file) != 1L) {
    stop("Cannot locate benchmark/run.R.", call. = FALSE)
}
benchmark_dir <- dirname(normalizePath(sub("^--file=", "", runner_file)))
repo_root <- dirname(benchmark_dir)
source(file.path(benchmark_dir, "protocol.R"))
source(file.path(benchmark_dir, "provenance.R"))
source(file.path(benchmark_dir, "fixtures.R"))
source(file.path(benchmark_dir, "report.R"))

benchmark_usage <- function() {
    cat(paste0(
        "Usage: Rscript benchmark/run.R [options]\n",
        "  --preset=full|smoke|hotspot  Scenario family (default: full)\n",
        "  --repeats=N                  Isolated repeats (default: 1)\n",
        "  --workers=N                  SOCK workers (default: 2)\n",
        "  --output=DIR                 Result directory\n",
        "  --help                       Show this text\n"
    ))
}

benchmark_parse_options <- function(args) {
    options <- list(preset = "full", repeats = 1L, workers = 2L, output = NULL)
    for (arg in args) {
        if (identical(arg, "--help")) {
            benchmark_usage()
            quit(save = "no", status = 0L)
        } else if (grepl("^--preset=", arg)) {
            options$preset <- sub("^--preset=", "", arg)
        } else if (grepl("^--repeats=", arg)) {
            options$repeats <- as.numeric(sub("^--repeats=", "", arg))
        } else if (grepl("^--workers=", arg)) {
            options$workers <- as.numeric(sub("^--workers=", "", arg))
        } else if (grepl("^--output=", arg)) {
            options$output <- sub("^--output=", "", arg)
        } else {
            stop("Unknown benchmark option: ", arg, call. = FALSE)
        }
    }
    options$repeats <- benchmark_assert_count(options$repeats, "`repeats`")
    options$workers <- benchmark_assert_count(options$workers, "`workers`")
    if (!options$preset %in% c("smoke", "full", "hotspot")) {
        stop("`preset` must be smoke, full, or hotspot.", call. = FALSE)
    }
    options
}

benchmark_summary <- function(results) {
    groups <- split(results, results$scenario_id, drop = TRUE)
    rows <- lapply(groups, function(group) {
        data.frame(
            scenario_id = group$scenario_id[[1L]],
            fixture_id = group$fixture_id[[1L]],
            storage = group$storage[[1L]],
            overlap = group$overlap[[1L]],
            backend = group$backend[[1L]],
            backend_workers = group$backend_workers[[1L]],
            n_features = group$n_features[[1L]],
            n_samples = group$n_samples[[1L]],
            n_sets = group$n_sets[[1L]],
            set_size = group$set_size[[1L]],
            density = group$density[[1L]],
            repeats = nrow(group),
            elapsed_min_sec = min(group$elapsed_sec),
            elapsed_median_sec = stats::median(group$elapsed_sec),
            score_cells_per_sec_median = stats::median(group$score_cells_per_sec),
            member_cells_per_sec_median = stats::median(group$member_cells_per_sec),
            manager_r_alloc_bytes_median = stats::median(group$manager_r_alloc_bytes),
            score_peak_tree_rss_kib_max = if (all(is.na(group$score_peak_tree_rss_kib))) {
                NA_real_
            } else {
                max(group$score_peak_tree_rss_kib, na.rm = TRUE)
            },
            score_peak_increment_kib_max = if (all(is.na(group$score_peak_increment_kib))) {
                NA_real_
            } else {
                max(group$score_peak_increment_kib, na.rm = TRUE)
            },
            process_max_rss_kib_max = if (all(is.na(group$process_max_rss_kib))) {
                NA_real_
            } else {
                max(group$process_max_rss_kib, na.rm = TRUE)
            },
            input_bytes = group$input_bytes[[1L]],
            logical_dense_bytes = group$logical_dense_bytes[[1L]],
            output_bytes = group$output_bytes[[1L]],
            output_md5 = paste(unique(group$output_md5), collapse = ","),
            stringsAsFactors = FALSE,
            check.names = FALSE
        )
    })
    result <- do.call(rbind, rows)
    rownames(result) <- NULL
    result
}

benchmark_verify_digests <- function(results) {
    repeated <- split(results$output_md5, results$scenario_id, drop = TRUE)
    unstable <- names(Filter(function(digests) length(unique(digests)) != 1L, repeated))
    paired_key <- paste(results$fixture_id, results$repeat_id, sep = "::")
    paired <- split(results$output_md5, paired_key, drop = TRUE)
    divergent <- names(Filter(function(digests) length(unique(digests)) != 1L, paired))
    if (length(unstable) > 0L || length(divergent) > 0L) {
        stop(
            "Output digest mismatch: ",
            paste(unique(c(unstable, divergent)), collapse = ", "),
            call. = FALSE
        )
    }
    invisible(TRUE)
}

options <- benchmark_parse_options(commandArgs(trailingOnly = TRUE))
if (is.null(options$output)) {
    stamp <- format(Sys.time(), "%Y%m%d-%H%M%S", tz = "UTC")
    options$output <- file.path(benchmark_dir, "results", paste(stamp, options$preset, sep = "-"))
}
dir.create(options$output, recursive = TRUE, showWarnings = FALSE)
output_dir <- normalizePath(options$output, mustWork = TRUE)
run_dir <- file.path(output_dir, "runs")
dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)

packages <- c("genefunnel", "Matrix", "BiocParallel", "Rcpp", "RcppArmadillo")
for (package in packages) {
    if (!requireNamespace(package, quietly = TRUE)) {
        stop("Benchmark dependency is not installed: ", package, call. = FALSE)
    }
}

protocol <- benchmark_read_protocol(file.path(benchmark_dir, "protocol.tsv"))
performance_protocol <- benchmark_protocol_subset(protocol, "performance")
if (!options$preset %in% performance_protocol$scenario_id) {
    stop("Benchmark preset is absent from protocol.tsv: ", options$preset, call. = FALSE)
}
scenarios <- benchmark_scenarios(options$preset, options$workers)
scenarios <- cbind(
    protocol_version = GENEFUNNEL_BENCHMARK_PROTOCOL,
    scenarios
)
benchmark_write_tsv(protocol, file.path(output_dir, "protocol.tsv"))
benchmark_write_tsv(scenarios, file.path(output_dir, "manifest.tsv"))
time_candidate <- if (file.exists("/usr/bin/time")) "/usr/bin/time" else ""
time_version <- if (nzchar(time_candidate)) {
    benchmark_command_output(time_candidate, "--version")
} else {
    character()
}
time_bin <- if (any(grepl("GNU Time", time_version, fixed = TRUE))) {
    time_candidate
} else {
    ""
}
metadata <- benchmark_metadata(
    repo_root = repo_root,
    runner = "benchmark/run.R",
    preset = options$preset,
    repeats = options$repeats,
    workers = options$workers,
    time_bin = time_bin
)
benchmark_write_tsv(metadata, file.path(output_dir, "metadata.tsv"))
benchmark_write_session_info(file.path(output_dir, "session-info.txt"), packages)

rscript <- Sys.which("Rscript")
worker <- file.path(benchmark_dir, "worker.R")
results <- list()
run_number <- 0L
for (scenario_index in seq_len(nrow(scenarios))) {
    for (repeat_id in seq_len(options$repeats)) {
        run_number <- run_number + 1L
        scenario <- as.list(scenarios[scenario_index, , drop = FALSE])
        scenario$repeat_id <- repeat_id
        prefix <- sprintf("%03d-%s-r%02d", run_number, scenario$scenario_id, repeat_id)
        scenario_path <- file.path(run_dir, paste0(prefix, "-scenario.rds"))
        result_path <- file.path(run_dir, paste0(prefix, "-result.tsv"))
        allocation_path <- file.path(run_dir, paste0(prefix, "-allocations.log"))
        stdout_path <- file.path(run_dir, paste0(prefix, "-stdout.log"))
        stderr_path <- file.path(run_dir, paste0(prefix, "-stderr.log"))
        time_path <- file.path(run_dir, paste0(prefix, "-time.txt"))
        saveRDS(scenario, scenario_path, version = 3L)

        worker_args <- c(
            "--vanilla",
            shQuote(worker),
            shQuote(scenario_path),
            shQuote(result_path),
            shQuote(allocation_path)
        )
        status <- if (nzchar(time_bin)) {
            system2(
                time_bin,
                c(
                    "--format=%M",
                    shQuote(paste0("--output=", time_path)),
                    "--",
                    shQuote(rscript),
                    worker_args
                ),
                stdout = stdout_path,
                stderr = stderr_path,
                env = paste0("GENEFUNNEL_BENCHMARK_TOKEN=", prefix)
            )
        } else {
            system2(
                rscript,
                worker_args,
                stdout = stdout_path,
                stderr = stderr_path,
                env = paste0("GENEFUNNEL_BENCHMARK_TOKEN=", prefix)
            )
        }
        if (!identical(status, 0L) || !file.exists(result_path)) {
            logs <- c(
                readLines(stdout_path, warn = FALSE),
                readLines(stderr_path, warn = FALSE)
            )
            stop(
                "Benchmark worker failed for ", scenario$scenario_id, ":\n",
                paste(utils::tail(logs, 30L), collapse = "\n"),
                call. = FALSE
            )
        }

        row <- utils::read.delim(
            result_path,
            stringsAsFactors = FALSE,
            check.names = FALSE,
            na.strings = "NA"
        )
        row$process_max_rss_kib <- if (file.exists(time_path)) {
            suppressWarnings(as.numeric(readLines(time_path, warn = FALSE)[[1L]]))
        } else {
            NA_real_
        }
        results[[length(results) + 1L]] <- row
        combined <- do.call(rbind, results)
        benchmark_write_tsv(combined, file.path(output_dir, "runs.tsv"))
        cat(readLines(stdout_path, warn = FALSE), sep = "\n")
    }
}

results <- do.call(rbind, results)
rownames(results) <- NULL
benchmark_verify_digests(results)
summary <- benchmark_summary(results)
benchmark_write_tsv(summary, file.path(output_dir, "summary.tsv"))
benchmark_write_performance_report(
    file.path(output_dir, "report.md"),
    protocol,
    summary,
    metadata,
    options
)
cat("Benchmark complete: ", output_dir, "\n", sep = "")
