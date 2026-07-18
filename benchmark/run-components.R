# Assisted-by: OpenAI Codex.

if (!nzchar(Sys.getenv("TZ", unset = ""))) {
    Sys.setenv(TZ = "UTC")
}

runner_file <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(runner_file) != 1L) {
    stop("Cannot locate benchmark/run-components.R.", call. = FALSE)
}
benchmark_dir <- dirname(normalizePath(sub("^--file=", "", runner_file)))
repo_root <- dirname(benchmark_dir)
source(file.path(benchmark_dir, "fixtures.R"))
source(file.path(benchmark_dir, "protocol.R"))
source(file.path(benchmark_dir, "provenance.R"))
source(file.path(benchmark_dir, "report.R"))
source(file.path(benchmark_dir, "components-installation.R"))

COMPONENT_PROTOCOL_VERSION <- "A2-1.0.0"
COMPONENT_PROTOCOL_MD5 <- "974e3aeab67a8a5c21a01cd0f154a4ed"
COMPONENT_BASELINE_ID <- "9b60a3eb138e5fd267586624ccd8bf51907577e7"

components_usage <- function() {
    cat(paste0(
        "Usage: Rscript benchmark/run-components.R [options]\n",
        "  --baseline-library=DIR  Installed baseline-only library\n",
        "  --candidate-library=DIR Installed candidate-only library\n",
        "  --baseline-id=ID        Baseline Git SHA (required in gate mode)\n",
        "  --candidate-id=ID       Candidate Git SHA (required in gate mode)\n",
        "  --mode=smoke|gate       Four-repeat smoke or fixed gate (default: smoke)\n",
        "  --output=DIR            Result directory\n",
        "  --help                  Show this text\n"
    ))
}

components_parse_options <- function(args) {
    options <- list(
        baseline_library = NULL,
        candidate_library = NULL,
        baseline_id = "smoke",
        candidate_id = "smoke",
        mode = "smoke",
        output = NULL
    )
    for (arg in args) {
        if (identical(arg, "--help")) {
            components_usage()
            quit(save = "no", status = 0L)
        }
        matched <- FALSE
        for (key in names(options)) {
            flag <- paste0("--", gsub("_", "-", key), "=")
            if (startsWith(arg, flag)) {
                options[[key]] <- substring(arg, nchar(flag) + 1L)
                matched <- TRUE
                break
            }
        }
        if (!matched) {
            stop("Unknown benchmark option: ", arg, call. = FALSE)
        }
    }
    if (is.null(options$baseline_library) || is.null(options$candidate_library)) {
        stop("Both package libraries are required.", call. = FALSE)
    }
    if (!options$mode %in% c("smoke", "gate")) {
        stop("`--mode` must be smoke or gate.", call. = FALSE)
    }
    if (options$mode == "gate") {
        if (!identical(options$baseline_id, COMPONENT_BASELINE_ID)) {
            stop("Gate mode requires the locked baseline Git SHA.", call. = FALSE)
        }
        if (!grepl("^[0-9a-f]{40}$", options$candidate_id)) {
            stop("Gate mode requires a full candidate Git SHA.", call. = FALSE)
        }
    }
    options
}

components_read_protocol <- function(path, mode) {
    protocol_md5 <- unname(tools::md5sum(path))
    if (!identical(protocol_md5, COMPONENT_PROTOCOL_MD5)) {
        stop(
            "Component performance protocol bytes changed; version the protocol.",
            call. = FALSE
        )
    }
    protocol <- utils::read.delim(
        path,
        stringsAsFactors = FALSE,
        check.names = FALSE,
        na.strings = "NA"
    )
    required <- c(
        "protocol_version", "scenario_id", "storage", "overlap", "backend",
        "backend_workers", "n_features", "n_samples", "n_sets", "set_size",
        "matrix_seed", "set_seed", "density", "zero_fraction",
        "cell_missing_fraction", "stored_missing_fraction", "timed_calls",
        "order_code"
    )
    if (!identical(names(protocol), required)) {
        stop("Component performance protocol fields changed.", call. = FALSE)
    }
    if (!all(protocol$protocol_version == COMPONENT_PROTOCOL_VERSION) ||
        anyDuplicated(protocol$scenario_id)) {
        stop("Component performance protocol identity is invalid.", call. = FALSE)
    }
    if (any(protocol$timed_calls != 5L)) {
        stop("Component performance timed-call count is invalid.", call. = FALSE)
    }
    valid_order <- nchar(protocol$order_code) == 30L &
        grepl("^[BC]+$", protocol$order_code) &
        vapply(strsplit(protocol$order_code, ""), function(order) {
            sum(order == "B") == 15L && sum(order == "C") == 15L
        }, logical(1L))
    if (!all(valid_order)) {
        stop("Component performance order schedule is invalid.", call. = FALSE)
    }
    if (mode == "smoke") {
        dense <- protocol$storage == "dense"
        protocol$n_features <- 400L
        protocol$n_samples <- ifelse(dense, 16L, 32L)
        protocol$n_sets <- 20L
        protocol$set_size <- 8L
    }
    protocol
}

components_package_record <- function(library, label, expected_id, require_marker) {
    library <- normalizePath(library, mustWork = TRUE)
    package_path <- normalizePath(
        file.path(library, "genefunnel"),
        mustWork = TRUE
    )
    description <- utils::packageDescription("genefunnel", lib.loc = library)
    marker <- file.path(library, COMPONENT_INSTALLATION_MARKER)
    provenance <- if (file.exists(marker)) {
        components_read_installation(library)
    } else if (require_marker) {
        stop("Gate package library lacks installation provenance.", call. = FALSE)
    } else {
        NULL
    }
    if (require_marker) {
        value <- provenance
        source_archive <- file.path(
            dirname(library),
            paste0(label, "-source.tar")
        )
        source_digest <- if (file.exists(source_archive)) {
            unname(tools::md5sum(source_archive))
        } else {
            NA_character_
        }
        valid <- value$protocol_version[[1L]] == COMPONENT_PROTOCOL_VERSION &&
            value$method[[1L]] == label &&
            value$git_sha[[1L]] == expected_id &&
            value$source_archive_md5[[1L]] == source_digest &&
            value$source_archive_md5[[1L]] ==
                components_git_archive_digest(repo_root, expected_id) &&
            value$package_version[[1L]] == description[["Version"]] &&
            value$installed_manifest_md5[[1L]] ==
                components_manifest_digest(package_path)
        if (!isTRUE(valid)) {
            stop("Component package installation provenance is invalid.", call. = FALSE)
        }
    }
    data.frame(
        method = label,
        library = library,
        package_path = package_path,
        package_version = description[["Version"]],
        git_sha = if (is.null(provenance)) NA_character_ else provenance$git_sha,
        installed_manifest_md5 = if (is.null(provenance)) {
            NA_character_
        } else {
            provenance$installed_manifest_md5
        },
        stringsAsFactors = FALSE,
        check.names = FALSE
    )
}

components_load_state <- function() {
    load_line <- tryCatch(
        readLines("/proc/loadavg", n = 1L, warn = FALSE),
        error = function(error) character()
    )
    logical_cores <- parallel::detectCores(logical = TRUE)
    load_1_min <- if (length(load_line) == 1L) {
        suppressWarnings(as.numeric(strsplit(load_line, "[[:space:]]+")[[1L]][[1L]]))
    } else {
        NA_real_
    }
    normalized <- load_1_min / logical_cores
    data.frame(
        checked_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
        load_1_min = load_1_min,
        logical_cores = logical_cores,
        load_per_core = normalized,
        limit = 0.25,
        pass = is.finite(normalized) && normalized <= 0.25,
        stringsAsFactors = FALSE,
        check.names = FALSE
    )
}

components_summary <- function(results, mode) {
    baseline <- results[results$method == "baseline", , drop = FALSE]
    candidate <- results[results$method == "candidate", , drop = FALSE]
    keys <- c("scenario_id", "repeat_id")
    paired <- merge(
        baseline,
        candidate,
        by = keys,
        suffixes = c("_baseline", "_candidate"),
        sort = FALSE
    )
    if (nrow(paired) * 2L != nrow(results)) {
        stop("Component benchmark pairs are incomplete.", call. = FALSE)
    }
    paired$elapsed_ratio <- paired$elapsed_sec_candidate /
        paired$elapsed_sec_baseline
    paired$cold_elapsed_ratio <- paired$cold_elapsed_sec_candidate /
        paired$cold_elapsed_sec_baseline
    paired$digest_equal <- paired$output_md5_candidate == paired$output_md5_baseline
    paired$output_shape_equal <-
        paired$output_bytes_candidate == paired$output_bytes_baseline &
        paired$output_na_cells_candidate == paired$output_na_cells_baseline &
        paired$output_sum_candidate == paired$output_sum_baseline
    environment_fields <- c(
        "r_version", "Matrix_version", "Matrix_path",
        "BiocParallel_version", "BiocParallel_path", "Rcpp_version",
        "Rcpp_path", "RcppArmadillo_version", "RcppArmadillo_path"
    )
    paired$environment_equal <- Reduce(`&`, lapply(environment_fields, function(field) {
        paired[[paste0(field, "_baseline")]] ==
            paired[[paste0(field, "_candidate")]]
    }))

    groups <- split(paired, paired$scenario_id, drop = TRUE)
    summary <- do.call(rbind, lapply(groups, function(group) {
        log_ratios <- log(group$elapsed_ratio)
        cold_log_ratios <- log(group$cold_elapsed_ratio)
        estimate <- exp(mean(log_ratios))
        standard_error <- stats::sd(log_ratios) / sqrt(nrow(group))
        upper <- if (is.finite(standard_error)) {
            exp(mean(log_ratios) +
                stats::qt(0.95, df = nrow(group) - 1L) * standard_error)
        } else {
            estimate
        }
        cold_estimate <- exp(mean(cold_log_ratios))
        cold_standard_error <- stats::sd(cold_log_ratios) / sqrt(nrow(group))
        cold_upper <- if (is.finite(cold_standard_error)) {
            exp(mean(cold_log_ratios) +
                stats::qt(0.95, df = nrow(group) - 1L) * cold_standard_error)
        } else {
            cold_estimate
        }
        allocation_baseline <- stats::median(group$manager_r_alloc_bytes_baseline)
        allocation_candidate <- stats::median(group$manager_r_alloc_bytes_candidate)
        allocation_limit <- allocation_baseline + max(65536, 0.01 * allocation_baseline)
        rss_baseline <- stats::median(group$process_max_rss_kib_baseline)
        rss_candidate <- stats::median(group$process_max_rss_kib_candidate)
        rss_limit <- rss_baseline + max(4096, 0.10 * rss_baseline)
        data.frame(
            scenario_id = group$scenario_id[[1L]],
            storage = group$storage_baseline[[1L]],
            overlap = group$overlap_baseline[[1L]],
            repeats = nrow(group),
            elapsed_ratio_geomean = estimate,
            elapsed_ratio_upper_95 = upper,
            cold_elapsed_ratio_geomean = cold_estimate,
            cold_elapsed_ratio_upper_95 = cold_upper,
            baseline_elapsed_median_sec = stats::median(group$elapsed_sec_baseline),
            candidate_elapsed_median_sec = stats::median(group$elapsed_sec_candidate),
            baseline_cold_elapsed_median_sec = stats::median(
                group$cold_elapsed_sec_baseline
            ),
            candidate_cold_elapsed_median_sec = stats::median(
                group$cold_elapsed_sec_candidate
            ),
            baseline_allocation_median_bytes = allocation_baseline,
            candidate_allocation_median_bytes = allocation_candidate,
            allocation_limit_bytes = allocation_limit,
            baseline_process_max_median_kib = rss_baseline,
            candidate_process_max_median_kib = rss_candidate,
            process_max_limit_kib = rss_limit,
            digest_pass = isTRUE(all(group$digest_equal)),
            output_shape_pass = isTRUE(all(group$output_shape_equal)),
            environment_pass = isTRUE(all(group$environment_equal)),
            timing_pass = isTRUE(mode == "smoke") ||
                (is.finite(upper) && upper <= 1.05),
            allocation_pass = isTRUE(mode == "smoke") ||
                (is.finite(allocation_candidate) &&
                    allocation_candidate <= allocation_limit),
            rss_pass = isTRUE(mode == "smoke") ||
                (is.finite(rss_candidate) && rss_candidate <= rss_limit),
            stringsAsFactors = FALSE,
            check.names = FALSE
        )
    }))
    rownames(summary) <- NULL
    list(pairs = paired, summary = summary)
}

components_write_report <- function(path, summary, metadata, options) {
    view <- data.frame(
        Scenario = summary$scenario_id,
        Storage = summary$storage,
        Overlap = summary$overlap,
        Repeats = summary$repeats,
        `Elapsed ratio` = benchmark_format_number(summary$elapsed_ratio_geomean),
        `One-sided 95% upper` = benchmark_format_number(
            summary$elapsed_ratio_upper_95
        ),
        `Cold ratio context` = benchmark_format_number(
            summary$cold_elapsed_ratio_geomean
        ),
        `Allocation pass` = summary$allocation_pass,
        `RSS pass` = summary$rss_pass,
        `Identity pass` = summary$digest_pass & summary$output_shape_pass,
        `Environment pass` = summary$environment_pass,
        check.names = FALSE
    )
    lines <- c(
        "# GeneFunnel component default-path comparison",
        "",
        paste0(
            "Protocol `", COMPONENT_PROTOCOL_VERSION, "`; mode `", options$mode,
            "`; baseline `", options$baseline_id, "`; candidate `",
            options$candidate_id, "`."
        ),
        "",
        benchmark_markdown_table(view),
        "",
        "## Environment",
        "",
        benchmark_markdown_table(metadata),
        ""
    )
    writeLines(lines, path)
}

options <- components_parse_options(commandArgs(trailingOnly = TRUE))
repeats <- if (options$mode == "gate") 30L else 4L
git_head <- benchmark_command_output("git", c("rev-parse", "HEAD"), repo_root)
git_status <- benchmark_command_output("git", c("status", "--porcelain"), repo_root)
if (options$mode == "gate" &&
    (length(git_head) != 1L || !identical(git_head[[1L]], options$candidate_id) ||
        length(git_status) > 0L)) {
    stop(
        "Gate mode requires a clean repository at the candidate Git SHA.",
        call. = FALSE
    )
}
if (is.null(options$output)) {
    stamp <- format(Sys.time(), "%Y%m%d-%H%M%S", tz = "UTC")
    options$output <- file.path(
        benchmark_dir,
        "results",
        paste(stamp, "components", options$mode, sep = "-")
    )
}
dir.create(options$output, recursive = TRUE, showWarnings = FALSE)
output_dir <- normalizePath(options$output, mustWork = TRUE)
run_dir <- file.path(output_dir, "runs")
dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)

libraries <- rbind(
    components_package_record(
        options$baseline_library,
        "baseline",
        options$baseline_id,
        options$mode == "gate"
    ),
    components_package_record(
        options$candidate_library,
        "candidate",
        options$candidate_id,
        options$mode == "gate"
    )
)
if (options$mode == "gate" &&
    identical(libraries$library[[1L]], libraries$library[[2L]])) {
    stop("Gate mode requires separate package libraries.", call. = FALSE)
}
protocol_path <- file.path(benchmark_dir, "components-protocol.tsv")
protocol <- components_read_protocol(protocol_path, options$mode)
benchmark_write_tsv(protocol, file.path(output_dir, "protocol.tsv"))
benchmark_write_tsv(libraries, file.path(output_dir, "packages.tsv"))

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
if (options$mode == "gate" &&
    (!nzchar(time_bin) || !isTRUE(capabilities("profmem")))) {
    stop(
        "Gate mode requires compatible GNU time and R allocation profiling.",
        call. = FALSE
    )
}
metadata <- data.frame(
    key = c(
        "generated_utc", "protocol_version", "mode", "repeats",
        "baseline_id", "candidate_id", "git_head", "git_dirty",
        "r_version", "r_platform", "hostname", "external_time"
    ),
    value = c(
        format(Sys.time(), tz = "UTC", usetz = TRUE),
        COMPONENT_PROTOCOL_VERSION,
        options$mode,
        repeats,
        options$baseline_id,
        options$candidate_id,
        if (length(git_head)) git_head[[1L]] else NA_character_,
        length(git_status) > 0L,
        paste(R.version$major, R.version$minor, sep = "."),
        R.version$platform,
        unname(Sys.info()[["nodename"]]),
        if (length(time_version)) time_version[[1L]] else NA_character_
    ),
    stringsAsFactors = FALSE,
    check.names = FALSE
)
benchmark_write_tsv(metadata, file.path(output_dir, "metadata.tsv"))
benchmark_write_session_info(
    file.path(output_dir, "session-info.txt"),
    c("Matrix", "BiocParallel", "Rcpp", "RcppArmadillo")
)

load_checks <- list()
components_check_load <- function(label) {
    state <- components_load_state()
    state$position <- label
    load_checks[[length(load_checks) + 1L]] <<- state
    benchmark_write_tsv(
        do.call(rbind, load_checks),
        file.path(output_dir, "load.tsv")
    )
    if (options$mode == "gate" && !state$pass) {
        stop(
            "Gate environment is not quiescent: one-minute load per logical ",
            "CPU must be at most 0.25.",
            call. = FALSE
        )
    }
    invisible(state)
}
components_check_load("start")

rscript <- Sys.which("Rscript")
worker <- file.path(benchmark_dir, "components-worker.R")
results <- list()
run_number <- 0L
for (repeat_id in seq_len(repeats)) {
    scenario_order <- ((seq_len(nrow(protocol)) + repeat_id - 2L) %%
        nrow(protocol)) + 1L
    for (scenario_index in scenario_order) {
        components_check_load(paste0(
            "before-r",
            repeat_id,
            "-",
            protocol$scenario_id[[scenario_index]]
        ))
        candidate_first <- substring(
            protocol$order_code[[scenario_index]],
            repeat_id,
            repeat_id
        ) == "C"
        method_order <- if (candidate_first) {
            c("candidate", "baseline")
        } else {
            c("baseline", "candidate")
        }
        for (method in method_order) {
            run_number <- run_number + 1L
            scenario <- as.list(protocol[scenario_index, , drop = FALSE])
            scenario$fixture_id <- paste(
                options$mode,
                scenario$storage,
                scenario$overlap,
                sep = "-"
            )
            scenario$method <- method
            scenario$repeat_id <- repeat_id
            prefix <- sprintf(
                "%03d-%s-r%02d-%s",
                run_number,
                scenario$scenario_id,
                repeat_id,
                method
            )
            scenario_path <- file.path(run_dir, paste0(prefix, "-scenario.rds"))
            result_path <- file.path(run_dir, paste0(prefix, "-result.tsv"))
            allocation_path <- file.path(run_dir, paste0(prefix, "-allocations.log"))
            stdout_path <- file.path(run_dir, paste0(prefix, "-stdout.log"))
            stderr_path <- file.path(run_dir, paste0(prefix, "-stderr.log"))
            time_path <- file.path(run_dir, paste0(prefix, "-time.txt"))
            saveRDS(scenario, scenario_path, version = 3L)
            library <- libraries$library[libraries$method == method]
            worker_args <- c(
                "--vanilla",
                shQuote(worker),
                shQuote(scenario_path),
                shQuote(result_path),
                shQuote(allocation_path),
                shQuote(library)
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
                    "Component benchmark worker failed for ", prefix, ":\n",
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
}
components_check_load("complete")

results <- do.call(rbind, results)
rownames(results) <- NULL
comparison <- components_summary(results, options$mode)
benchmark_write_tsv(comparison$pairs, file.path(output_dir, "pairs.tsv"))
benchmark_write_tsv(comparison$summary, file.path(output_dir, "summary.tsv"))
components_write_report(
    file.path(output_dir, "report.md"),
    comparison$summary,
    metadata,
    options
)
all_pass <- isTRUE(with(
    comparison$summary,
    all(
        digest_pass,
        output_shape_pass,
        environment_pass,
        timing_pass,
        allocation_pass,
        rss_pass
    )
))
decision <- data.frame(
    protocol_version = COMPONENT_PROTOCOL_VERSION,
    mode = options$mode,
    all_pass = all_pass,
    stringsAsFactors = FALSE,
    check.names = FALSE
)
benchmark_write_tsv(decision, file.path(output_dir, "decision.tsv"))
cat("Component comparison complete: ", output_dir, "\n", sep = "")
if (!all_pass) {
    quit(save = "no", status = 1L)
}
