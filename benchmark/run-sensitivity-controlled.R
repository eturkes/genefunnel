# Assisted-by: OpenAI Codex.

if (!nzchar(Sys.getenv("TZ", unset = ""))) Sys.setenv(TZ = "UTC")

runner_file <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(runner_file) != 1L) {
    stop("Cannot locate run-sensitivity-controlled.R.", call. = FALSE)
}
benchmark_dir <- dirname(normalizePath(sub("^--file=", "", runner_file)))
repo_root <- dirname(benchmark_dir)
source(file.path(benchmark_dir, "protocol.R"))
source(file.path(benchmark_dir, "provenance.R"))
source(file.path(benchmark_dir, "report.R"))
source(file.path(benchmark_dir, "components-installation.R"))
source(file.path(benchmark_dir, "sensitivity-protocol.R"))
source(file.path(benchmark_dir, "sensitivity-controlled-protocol.R"))
source(file.path(benchmark_dir, "sensitivity-controlled.R"))
source(file.path(benchmark_dir, "sensitivity-controlled-summary.R"))
source(file.path(benchmark_dir, "sensitivity-controlled-runner.R"))

sensitivity_controlled_usage <- function() {
    cat(paste0(
        "Usage: Rscript benchmark/run-sensitivity-controlled.R [options]\n",
        "  --mode=smoke|gate       Orchestration smoke or frozen full gate\n",
        "  --candidate-id=SHA      Explicit clean current SHA required for gate\n",
        "  --workers=N             Fork workers on Unix (default: up to 4)\n",
        "  --chunk-size=N          Scenarios per checkpoint (default: 128)\n",
        "  --output=DIR            Result/checkpoint directory inside repository\n",
        "  --help                  Show this text\n"
    ))
}

sensitivity_controlled_positive_integer <- function(value, label) {
    parsed <- suppressWarnings(as.integer(value))
    if (length(parsed) != 1L || is.na(parsed) || parsed < 1L ||
        !identical(as.character(parsed), value)) {
        stop(label, " must be a positive integer.", call. = FALSE)
    }
    parsed
}

sensitivity_controlled_options <- function(args) {
    cores <- parallel::detectCores(logical = TRUE)
    workers <- if (is.finite(cores)) min(4L, cores) else 1L
    options <- list(
        mode = "smoke", candidate_id = NULL, workers = as.integer(workers),
        chunk_size = 128L, output = NULL
    )
    for (argument in args) {
        if (argument == "--help") {
            sensitivity_controlled_usage()
            quit(save = "no", status = 0L)
        } else if (startsWith(argument, "--mode=")) {
            options$mode <- sub("^--mode=", "", argument)
        } else if (startsWith(argument, "--candidate-id=")) {
            options$candidate_id <- sub("^--candidate-id=", "", argument)
        } else if (startsWith(argument, "--workers=")) {
            options$workers <- sensitivity_controlled_positive_integer(
                sub("^--workers=", "", argument), "`--workers`"
            )
        } else if (startsWith(argument, "--chunk-size=")) {
            options$chunk_size <- sensitivity_controlled_positive_integer(
                sub("^--chunk-size=", "", argument), "`--chunk-size`"
            )
        } else if (startsWith(argument, "--output=")) {
            options$output <- sub("^--output=", "", argument)
        } else stop("Unknown sensitivity controlled option: ", argument)
    }
    if (!options$mode %in% c("smoke", "gate")) {
        stop("`--mode` must be smoke or gate.", call. = FALSE)
    }
    if (.Platform$OS.type == "windows" && options$workers != 1L) {
        stop("Windows sensitivity execution requires `--workers=1`.")
    }
    options
}

options <- sensitivity_controlled_options(commandArgs(trailingOnly = TRUE))
parent <- sensitivity_validate_protocol(repo_root)
execution <- sensitivity_controlled_validate_protocol(repo_root)
registry <- sensitivity_read_registry(file.path(
    benchmark_dir, "sensitivity-protocol.tsv"
))
sensitivity_controlled_validate_endpoint_contract(registry)
git <- sensitivity_controlled_git_state(
    repo_root, options$candidate_id, options$mode
)
output_dir <- sensitivity_controlled_output_path(
    options$output, repo_root, options$mode
)

full_design <- sensitivity_controlled_design(registry)
design <- if (options$mode == "gate") full_design else
    sensitivity_controlled_smoke_design(full_design)
config <- sensitivity_controlled_run_config(
    git, options$mode, options$chunk_size, design
)
sensitivity_controlled_publish_tsv(
    config, file.path(output_dir, "run-config.tsv")
)
sensitivity_controlled_publish_tsv(
    design, file.path(output_dir, "manifest.tsv")
)
sensitivity_controlled_publish_file(
    file.path(benchmark_dir, "sensitivity-protocol.tsv"),
    file.path(output_dir, "parent-protocol.tsv"), SENSITIVITY_PROTOCOL_MD5
)
sensitivity_controlled_publish_file(
    file.path(benchmark_dir, "sensitivity-controlled-protocol.tsv"),
    file.path(output_dir, "protocol.tsv"), SENSITIVITY_CONTROLLED_MD5
)

package <- NULL
tryCatch({
    cat("Installing clean sensitivity candidate snapshot...\n")
    package <- sensitivity_controlled_install(
        repo_root, output_dir, git$candidate
    )
    backend <- BiocParallel::SerialParam(progressbar = FALSE)
    cat("Running four-stratum representation preflight...\n")
    sensitivity_controlled_preflight(
        full_design, registry, package$scorer, package$sensitivity, backend
    )

    cat(
        "Calculating ", nrow(design), " scenario(s) with ", options$workers,
        " worker(s)...\n", sep = ""
    )
    observed <- sensitivity_controlled_simulate_design(
        design, registry, package$scorer, package$sensitivity, backend,
        options$workers, options$chunk_size,
        file.path(output_dir, "checkpoints"),
        sensitivity_controlled_config_fingerprint(config)
    )
    if (options$mode == "gate") {
        sensitivity_controlled_gate_counts(observed, design, registry)
        frames <- observed
        repeats <- NULL
    } else {
        frames <- sensitivity_controlled_model_smoke(registry)
        repeats <- 8L
        benchmark_write_tsv(
            frames$feature, file.path(output_dir, "model-smoke-feature.tsv")
        )
        benchmark_write_tsv(
            frames$technical, file.path(output_dir, "model-smoke-technical.tsv")
        )
    }

    cat("Fitting frozen held-out models and fixed-prediction bootstrap...\n")
    analysis <- sensitivity_controlled_analyze(frames, registry, repeats)
    if (options$mode == "smoke") {
        sensitivity_controlled_validate_model_adversaries(
            frames, analysis$feature_cv, registry
        )
    }
    summary <- sensitivity_controlled_write_analysis(
        output_dir, observed, analysis, options$mode
    )
    metadata <- sensitivity_controlled_metadata(
        repo_root, git, package, options, observed, analysis
    )
    benchmark_write_tsv(metadata, file.path(output_dir, "metadata.tsv"))
    benchmark_write_session_info(
        file.path(output_dir, "session-info.txt"),
        c("genefunnel", "Matrix", "BiocParallel", "Rcpp", "RcppArmadillo")
    )
    sensitivity_controlled_report(
        file.path(output_dir, "report.md"), analysis, summary,
        metadata, options, git
    )
    artifacts <- c(
        "run-config.tsv", "manifest.tsv", "parent-protocol.tsv", "protocol.tsv",
        "source.tar", "installation.log", "feature-observations.tsv",
        "technical-observations.tsv", "fold-results.tsv",
        "feature-predictions.tsv", "technical-predictions.tsv",
        "model-coefficients.tsv", "model-scaling.tsv", "bootstrap.tsv",
        "endpoints.tsv", "strata.tsv", "summary.tsv", "metadata.tsv",
        "session-info.txt", "report.md"
    )
    if (options$mode == "smoke") artifacts <- c(
        artifacts, "model-smoke-feature.tsv", "model-smoke-technical.tsv"
    )
    benchmark_write_tsv(
        sensitivity_controlled_artifacts(output_dir, artifacts),
        file.path(output_dir, "artifacts.tsv")
    )
    outcome <- if (options$mode == "smoke") "SMOKE" else if (
        summary$controlled_gate_pass
    ) "PASS" else "FAIL"
    cat(
        "Sensitivity controlled execution complete (", outcome, "): ",
        output_dir, "\n", sep = ""
    )
}, finally = {
    sensitivity_controlled_uninstall(package)
})
