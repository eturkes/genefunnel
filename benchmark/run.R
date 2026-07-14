if (!nzchar(Sys.getenv("TZ", unset = ""))) {
    Sys.setenv(TZ = "UTC")
}

runner_file <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(runner_file) != 1L) {
    stop("Cannot locate benchmark/run.R.", call. = FALSE)
}
benchmark_dir <- dirname(normalizePath(sub("^--file=", "", runner_file)))
repo_root <- dirname(benchmark_dir)
source(file.path(benchmark_dir, "fixtures.R"))

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

benchmark_command_output <- function(command, args = character(), wd = NULL) {
    old_wd <- getwd()
    if (!is.null(wd)) {
        setwd(wd)
        on.exit(setwd(old_wd), add = TRUE)
    }
    tryCatch(
        system2(command, args, stdout = TRUE, stderr = TRUE),
        error = function(error) character()
    )
}

benchmark_linux_field <- function(path, pattern) {
    lines <- tryCatch(readLines(path, warn = FALSE), error = function(error) character())
    value <- grep(pattern, lines, value = TRUE)
    if (length(value) == 0L) NA_character_ else trimws(sub(pattern, "", value[[1L]]))
}

benchmark_write_metadata <- function(path, options, time_bin) {
    git_head <- benchmark_command_output("git", c("rev-parse", "HEAD"), repo_root)
    git_status <- benchmark_command_output("git", c("status", "--porcelain"), repo_root)
    cpu_model <- benchmark_linux_field("/proc/cpuinfo", "^model name[[:space:]]*:[[:space:]]*")
    memory_kib <- benchmark_linux_field("/proc/meminfo", "^MemTotal:[[:space:]]*")
    time_version <- if (nzchar(time_bin)) {
        version <- benchmark_command_output(time_bin, "--version")
        if (length(version)) version[[1L]] else NA_character_
    } else {
        NA_character_
    }
    metadata <- data.frame(
        key = c(
            "generated_utc", "preset", "repeats", "requested_workers",
            "git_head", "git_dirty", "hostname", "sysname", "release",
            "machine", "cpu_model", "logical_cores", "memory_total",
            "rscript", "external_time"
        ),
        value = c(
            format(Sys.time(), tz = "UTC", usetz = TRUE),
            options$preset,
            options$repeats,
            options$workers,
            if (length(git_head)) git_head[[1L]] else NA_character_,
            length(git_status) > 0L,
            unname(Sys.info()[["nodename"]]),
            unname(Sys.info()[["sysname"]]),
            unname(Sys.info()[["release"]]),
            unname(Sys.info()[["machine"]]),
            cpu_model,
            parallel::detectCores(logical = TRUE),
            memory_kib,
            Sys.which("Rscript"),
            time_version
        ),
        stringsAsFactors = FALSE,
        check.names = FALSE
    )
    utils::write.table(
        metadata,
        path,
        sep = "\t",
        row.names = FALSE,
        quote = TRUE,
        na = "NA"
    )
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

if (!requireNamespace("genefunnel", quietly = TRUE)) {
    stop("Install GeneFunnel before running benchmarks.", call. = FALSE)
}
for (package in c("Matrix", "BiocParallel", "Rcpp", "RcppArmadillo")) {
    if (!requireNamespace(package, quietly = TRUE)) {
        stop("Benchmark dependency is not installed: ", package, call. = FALSE)
    }
}

scenarios <- benchmark_scenarios(options$preset, options$workers)
utils::write.table(
    scenarios,
    file.path(output_dir, "manifest.tsv"),
    sep = "\t",
    row.names = FALSE,
    quote = TRUE,
    na = "NA"
)
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
benchmark_write_metadata(file.path(output_dir, "metadata.tsv"), options, time_bin)
writeLines(
    capture.output({
        print(sessionInfo())
        cat("\nInstalled benchmark packages:\n")
        for (package in c(
            "genefunnel", "Matrix", "BiocParallel", "Rcpp", "RcppArmadillo"
        )) {
            cat(package, as.character(utils::packageVersion(package)), "\n")
        }
    }),
    file.path(output_dir, "session-info.txt")
)

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
        utils::write.table(
            combined,
            file.path(output_dir, "runs.tsv"),
            sep = "\t",
            row.names = FALSE,
            quote = TRUE,
            na = "NA"
        )
        cat(readLines(stdout_path, warn = FALSE), sep = "\n")
    }
}

results <- do.call(rbind, results)
rownames(results) <- NULL
benchmark_verify_digests(results)
summary <- benchmark_summary(results)
utils::write.table(
    summary,
    file.path(output_dir, "summary.tsv"),
    sep = "\t",
    row.names = FALSE,
    quote = TRUE,
    na = "NA"
)
cat("Benchmark complete: ", output_dir, "\n", sep = "")
