# Assisted-by: OpenAI Codex.

if (!nzchar(Sys.getenv("TZ", unset = ""))) {
    Sys.setenv(TZ = "UTC")
}

runner_file <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(runner_file) != 1L) {
    stop("Cannot locate benchmark/run-catalogue.R.", call. = FALSE)
}
benchmark_dir <- dirname(normalizePath(sub("^--file=", "", runner_file)))
repo_root <- dirname(benchmark_dir)
source(file.path(benchmark_dir, "fixtures.R"))
source(file.path(benchmark_dir, "protocol.R"))
source(file.path(benchmark_dir, "provenance.R"))
source(file.path(benchmark_dir, "report.R"))
source(file.path(benchmark_dir, "catalogue-fixtures.R"))
source(file.path(benchmark_dir, "catalogue-installation.R"))
source(file.path(benchmark_dir, "catalogue-summary.R"))

CATALOGUE_BASELINE_ID <- "a573c124909235e41bdbc3cfae950947465d8755"

catalogue_usage <- function() {
    cat(paste0(
        "Usage: Rscript benchmark/run-catalogue.R [options]\n",
        "  --baseline-library=DIR  Installed named-list control library\n",
        "  --candidate-library=DIR Installed compiled candidate library\n",
        "  --baseline-id=ID        Baseline Git SHA (required in gate mode)\n",
        "  --candidate-id=ID       Candidate Git SHA (required in gate mode)\n",
        "  --mode=smoke|gate       Four-repeat smoke or fixed gate\n",
        "  --output=DIR            Result directory\n",
        "  --help                  Show this text\n"
    ))
}

catalogue_parse_options <- function(args) {
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
            catalogue_usage()
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
            stop("Unknown catalogue benchmark option: ", arg, call. = FALSE)
        }
    }
    if (is.null(options$baseline_library) || is.null(options$candidate_library)) {
        stop("Both catalogue package libraries are required.", call. = FALSE)
    }
    if (!options$mode %in% c("smoke", "gate")) {
        stop("`--mode` must be smoke or gate.", call. = FALSE)
    }
    if (options$mode == "gate") {
        if (!identical(options$baseline_id, CATALOGUE_BASELINE_ID)) {
            stop("Gate mode requires the locked list baseline SHA.", call. = FALSE)
        }
        if (!grepl("^[0-9a-f]{40}$", options$candidate_id)) {
            stop("Gate mode requires a full candidate Git SHA.", call. = FALSE)
        }
    }
    options
}

catalogue_load_state <- function() {
    load_line <- tryCatch(
        readLines("/proc/loadavg", n = 1L, warn = FALSE),
        error = function(condition) character()
    )
    cores <- parallel::detectCores(logical = TRUE)
    load_1_min <- if (length(load_line) == 1L) {
        suppressWarnings(as.numeric(strsplit(
            load_line,
            "[[:space:]]+"
        )[[1L]][[1L]]))
    } else {
        NA_real_
    }
    normalized <- load_1_min / cores
    data.frame(
        checked_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
        load_1_min = load_1_min,
        logical_cores = cores,
        load_per_core = normalized,
        limit = 0.25,
        pass = is.finite(normalized) && normalized <= 0.25,
        stringsAsFactors = FALSE,
        check.names = FALSE
    )
}

catalogue_worker_env <- function(library, token) {
    paths <- unique(c(normalizePath(library, mustWork = TRUE), .libPaths()))
    c(
        paste0("R_LIBS=", paste(paths, collapse = .Platform$path.sep)),
        paste0("GENEFUNNEL_BENCHMARK_TOKEN=", token)
    )
}

options <- catalogue_parse_options(commandArgs(trailingOnly = TRUE))
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
        paste(stamp, "catalogue", options$mode, sep = "-")
    )
}
dir.create(options$output, recursive = TRUE, showWarnings = FALSE)
output_dir <- normalizePath(options$output, mustWork = TRUE)
run_dir <- file.path(output_dir, "runs")
dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)

libraries <- rbind(
    catalogue_package_record(
        repo_root,
        options$baseline_library,
        "list",
        options$baseline_id,
        CATALOGUE_PROTOCOL_VERSION,
        options$mode == "gate"
    ),
    catalogue_package_record(
        repo_root,
        options$candidate_library,
        "compiled",
        options$candidate_id,
        CATALOGUE_PROTOCOL_VERSION,
        options$mode == "gate"
    )
)
if (options$mode == "gate" &&
    identical(libraries$library[[1L]], libraries$library[[2L]])) {
    stop("Gate mode requires separate catalogue libraries.", call. = FALSE)
}
protocol_path <- file.path(benchmark_dir, "catalogue-protocol.tsv")
protocol <- catalogue_read_protocol(protocol_path, options$mode)
repeats <- unique(protocol$repeats)
batches <- unique(protocol$batches)
if (length(repeats) != 1L || length(batches) != 1L) {
    stop("Catalogue protocol repeat/batch counts are inconsistent.", call. = FALSE)
}
repeats <- as.integer(repeats)
batches <- as.integer(batches)
benchmark_write_tsv(
    utils::read.delim(
        protocol_path,
        stringsAsFactors = FALSE,
        check.names = FALSE,
        na.strings = "NA"
    ),
    file.path(output_dir, "protocol.tsv")
)
benchmark_write_tsv(protocol, file.path(output_dir, "manifest.tsv"))
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
    (.Platform$OS.type != "unix" || !dir.exists("/proc/self") ||
        !nzchar(time_bin) || !isTRUE(capabilities("profmem")))) {
    stop(
        "Gate mode requires Linux /proc, GNU time, and R allocation profiling.",
        call. = FALSE
    )
}
catalogue_package_path <- function(package) {
    tryCatch(
        normalizePath(find.package(package), mustWork = TRUE),
        error = function(condition) NA_character_
    )
}
system <- Sys.info()
metadata <- data.frame(
    key = c(
        "generated_utc", "protocol_version", "mode", "repeats", "batches",
        "baseline_id", "candidate_id", "git_head", "git_dirty", "hostname",
        "sysname", "release", "machine", "cpu_model", "logical_cores",
        "memory_total", "r_version", "r_platform", "rscript",
        "Matrix_version", "Matrix_path", "BiocParallel_version",
        "BiocParallel_path", "Rcpp_version", "Rcpp_path",
        "RcppArmadillo_version", "RcppArmadillo_path", "external_time_path",
        "external_time"
    ),
    value = c(
        format(Sys.time(), tz = "UTC", usetz = TRUE),
        CATALOGUE_PROTOCOL_VERSION,
        options$mode,
        repeats,
        batches,
        options$baseline_id,
        options$candidate_id,
        if (length(git_head)) git_head[[1L]] else NA_character_,
        length(git_status) > 0L,
        unname(system[["nodename"]]),
        unname(system[["sysname"]]),
        unname(system[["release"]]),
        unname(system[["machine"]]),
        benchmark_linux_field(
            "/proc/cpuinfo",
            "^model name[[:space:]]*:[[:space:]]*"
        ),
        parallel::detectCores(logical = TRUE),
        benchmark_linux_field(
            "/proc/meminfo",
            "^MemTotal:[[:space:]]*"
        ),
        paste(R.version$major, R.version$minor, sep = "."),
        R.version$platform,
        Sys.which("Rscript"),
        benchmark_package_version("Matrix"),
        catalogue_package_path("Matrix"),
        benchmark_package_version("BiocParallel"),
        catalogue_package_path("BiocParallel"),
        benchmark_package_version("Rcpp"),
        catalogue_package_path("Rcpp"),
        benchmark_package_version("RcppArmadillo"),
        catalogue_package_path("RcppArmadillo"),
        time_bin,
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

catalogue_write_decision <- function(all_pass, performance_claim, reason) {
    decision <- data.frame(
        protocol_version = CATALOGUE_PROTOCOL_VERSION,
        mode = options$mode,
        all_pass = isTRUE(all_pass),
        performance_claim = isTRUE(performance_claim),
        reason = reason,
        stringsAsFactors = FALSE,
        check.names = FALSE
    )
    benchmark_write_tsv(decision, file.path(output_dir, "decision.tsv"))
}

load_checks <- list()
catalogue_check_load <- function(label) {
    state <- catalogue_load_state()
    state$position <- label
    load_checks[[length(load_checks) + 1L]] <<- state
    benchmark_write_tsv(
        do.call(rbind, load_checks),
        file.path(output_dir, "load.tsv")
    )
    if (options$mode == "gate" && !state$pass) {
        catalogue_write_decision(
            FALSE,
            FALSE,
            paste0("quiescence_rejection:", label)
        )
        stop(
            "Gate environment is not quiescent: one-minute load per logical ",
            "CPU must be at most 0.25.",
            call. = FALSE
        )
    }
    invisible(state)
}
catalogue_check_load("start")

rscript <- Sys.which("Rscript")
candidate_library <- libraries$library[libraries$method == "compiled"]
correctness_path <- file.path(output_dir, "correctness.tsv")
correctness_stdout <- file.path(output_dir, "correctness-stdout.log")
correctness_stderr <- file.path(output_dir, "correctness-stderr.log")
correctness_status <- system2(
    rscript,
    c(
        "--vanilla",
        shQuote(file.path(benchmark_dir, "catalogue-correctness.R")),
        shQuote(candidate_library),
        shQuote(correctness_path)
    ),
    stdout = correctness_stdout,
    stderr = correctness_stderr,
    env = catalogue_worker_env(candidate_library, "catalogue-correctness")
)
if (!identical(correctness_status, 0L) || !file.exists(correctness_path)) {
    logs <- c(
        readLines(correctness_stdout, warn = FALSE),
        readLines(correctness_stderr, warn = FALSE)
    )
    stop(
        "Catalogue correctness worker failed:\n",
        paste(utils::tail(logs, 30L), collapse = "\n"),
        call. = FALSE
    )
}
correctness <- utils::read.delim(
    correctness_path,
    stringsAsFactors = FALSE,
    check.names = FALSE
)
correctness_pass <- isTRUE(all(correctness$passed))

catalogue_merge_phase_rows <- function(rows, scenario_fields) {
    expected_phases <- c("timing", "allocation", "rss")
    if (!identical(names(rows), expected_phases)) {
        stop("Catalogue worker phases are incomplete.", call. = FALSE)
    }
    shared <- c(
        scenario_fields,
        "output_sha256", "output_bytes", "output_na_cells", "output_sum",
        "package_library", "package_path", "genefunnel_version",
        "Matrix_version", "Matrix_path", "BiocParallel_version",
        "BiocParallel_path", "Rcpp_version", "Rcpp_path",
        "RcppArmadillo_version", "RcppArmadillo_path", "R_version",
        "R_platform"
    )
    metrics <- list(
        timing = c(
            "compile_elapsed_sec", "total_elapsed_sec",
            paste0("call_", seq_len(5L), "_sec"),
            paste0("cumulative_", seq_len(5L), "_sec")
        ),
        allocation = c(
            "validation_sec", "canonical_encoding_sec", "serialization_sec",
            "deserialization_sec", "operation_alloc_bytes", "scoring_alloc_bytes",
            "compile_alloc_bytes", "canonical_memberships", "object_bytes",
            "object_bytes_per_membership", "wire_bytes",
            "wire_bytes_per_membership", "catalogue_sha256",
            "feature_sha256", "content_sha256", "wire_sha256"
        ),
        rss = c(
            "compile_baseline_rss_kib", "compile_peak_tree_rss_kib",
            "compile_peak_increment_kib", "compile_peak_processes",
            "compile_memory_samples", "compile_memory_scope",
            "scoring_baseline_rss_kib", "scoring_peak_tree_rss_kib",
            "scoring_peak_increment_kib", "scoring_peak_processes",
            "scoring_memory_samples", "scoring_memory_scope"
        )
    )
    for (phase in expected_phases) {
        row <- rows[[phase]]
        expected <- c(shared, "phase", metrics[[phase]])
        if (nrow(row) != 1L || anyDuplicated(names(row)) ||
            !setequal(names(row), expected) || !identical(row$phase, phase)) {
            stop("Catalogue worker phase schema is invalid: ", phase, call. = FALSE)
        }
    }
    reference <- as.list(rows$timing[1L, shared, drop = FALSE])
    for (phase in expected_phases[-1L]) {
        observed <- as.list(rows[[phase]][1L, shared, drop = FALSE])
        if (!identical(observed, reference)) {
            stop(
                "Catalogue worker phases disagree on inputs, outputs, or environment.",
                call. = FALSE
            )
        }
    }
    merged <- rows$timing[1L, shared, drop = FALSE]
    for (phase in expected_phases) {
        merged <- cbind(
            merged,
            rows[[phase]][1L, metrics[[phase]], drop = FALSE]
        )
    }
    merged
}

worker <- file.path(benchmark_dir, "catalogue-worker.R")
catalogue_run_worker_phase <- function(context, phase) {
    phase_prefix <- paste(context$prefix, phase, sep = "-")
    result_path <- file.path(run_dir, paste0(phase_prefix, "-result.tsv"))
    allocation_prefix <- file.path(
        run_dir,
        paste0(phase_prefix, "-allocations")
    )
    stdout_path <- file.path(run_dir, paste0(phase_prefix, "-stdout.log"))
    stderr_path <- file.path(run_dir, paste0(phase_prefix, "-stderr.log"))
    worker_args <- c(
        "--vanilla",
        shQuote(worker),
        shQuote(context$scenario_path),
        shQuote(result_path),
        shQuote(allocation_prefix),
        shQuote(context$library),
        context$method,
        phase
    )
    environment <- catalogue_worker_env(context$library, phase_prefix)
    status <- if (phase == "timing" && nzchar(time_bin)) {
        system2(
            time_bin,
            c(
                "--format=%M",
                shQuote(paste0("--output=", context$time_path)),
                "--",
                shQuote(rscript),
                worker_args
            ),
            stdout = stdout_path,
            stderr = stderr_path,
            env = environment
        )
    } else {
        system2(
            rscript,
            worker_args,
            stdout = stdout_path,
            stderr = stderr_path,
            env = environment
        )
    }
    if (!identical(status, 0L) || !file.exists(result_path)) {
        logs <- c(
            readLines(stdout_path, warn = FALSE),
            readLines(stderr_path, warn = FALSE)
        )
        stop(
            "Catalogue worker failed for ", phase_prefix, ":\n",
            paste(utils::tail(logs, 30L), collapse = "\n"),
            call. = FALSE
        )
    }
    cat(readLines(stdout_path, warn = FALSE), sep = "\n")
    utils::read.delim(
        result_path,
        stringsAsFactors = FALSE,
        check.names = FALSE,
        na.strings = "NA"
    )
}

results <- list()
run_number <- 0L
for (repeat_id in seq_len(repeats)) {
    scenario_order <- ((seq_len(nrow(protocol)) + repeat_id - 2L) %%
        nrow(protocol)) + 1L
    for (scenario_index in scenario_order) {
        scenario_id <- protocol$scenario_id[[scenario_index]]
        catalogue_check_load(paste0("before-r", repeat_id, "-", scenario_id))
        compiled_first <- substring(
            protocol$order_code[[scenario_index]],
            repeat_id,
            repeat_id
        ) == "C"
        method_order <- if (compiled_first) {
            c("compiled", "list")
        } else {
            c("list", "compiled")
        }
        contexts <- list()
        for (method in method_order) {
            run_number <- run_number + 1L
            scenario <- as.list(protocol[scenario_index, , drop = FALSE])
            scenario$method <- method
            scenario$repeat_id <- repeat_id
            prefix <- sprintf(
                "%03d-%s-r%02d-%s",
                run_number,
                scenario_id,
                repeat_id,
                method
            )
            scenario_path <- file.path(run_dir, paste0(prefix, "-scenario.rds"))
            time_path <- file.path(run_dir, paste0(prefix, "-time.txt"))
            saveRDS(scenario, scenario_path, version = 3L)
            library <- libraries$library[libraries$method == method]
            contexts[[method]] <- list(
                prefix = prefix,
                scenario = scenario,
                scenario_path = scenario_path,
                time_path = time_path,
                library = library,
                method = method
            )
        }
        phase_rows <- setNames(lapply(method_order, function(method) {
            list()
        }), method_order)
        for (phase in c("timing", "allocation", "rss")) {
            for (method in method_order) {
                phase_rows[[method]][[phase]] <- catalogue_run_worker_phase(
                    contexts[[method]],
                    phase
                )
            }
        }
        for (method in method_order) {
            context <- contexts[[method]]
            row <- catalogue_merge_phase_rows(
                phase_rows[[method]],
                names(context$scenario)
            )
            row$process_max_rss_kib <- if (file.exists(context$time_path)) {
                values <- suppressWarnings(as.numeric(readLines(
                    context$time_path,
                    warn = FALSE
                )))
                if (length(values)) values[[1L]] else NA_real_
            } else {
                NA_real_
            }
            results[[length(results) + 1L]] <- row
            benchmark_write_tsv(
                do.call(rbind, results),
                file.path(output_dir, "runs.tsv")
            )
        }
    }
}
catalogue_check_load("complete")

results <- do.call(rbind, results)
rownames(results) <- NULL
paired <- catalogue_pair_results(results, options$mode, batches)
summary <- catalogue_scenario_summary(paired, options$mode, batches)
benchmark_write_tsv(paired, file.path(output_dir, "pairs.tsv"))
benchmark_write_tsv(summary, file.path(output_dir, "summary.tsv"))
catalogue_write_report(
    file.path(output_dir, "report.md"),
    summary,
    metadata,
    correctness,
    options
)
all_pass <- correctness_pass && catalogue_summary_pass(summary)
catalogue_write_decision(
    all_pass,
    identical(options$mode, "gate") && all_pass,
    if (all_pass) {
        if (identical(options$mode, "gate")) {
            "all_gates_pass"
        } else {
            "smoke_orchestration_pass"
        }
    } else {
        "correctness_or_threshold_failure"
    }
)
cat("Catalogue comparison complete: ", output_dir, "\n", sep = "")
if (!all_pass) {
    quit(save = "no", status = 1L)
}
