# Assisted-by: OpenAI Codex.

if (!nzchar(Sys.getenv("TZ", unset = ""))) {
    Sys.setenv(TZ = "UTC")
}

runner_file <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(runner_file) != 1L) {
    stop("Cannot locate run-aggregation-synthetic.R.", call. = FALSE)
}
benchmark_dir <- dirname(normalizePath(sub("^--file=", "", runner_file)))
repo_root <- dirname(benchmark_dir)
source(file.path(benchmark_dir, "protocol.R"))
source(file.path(benchmark_dir, "provenance.R"))
source(file.path(benchmark_dir, "report.R"))
source(file.path(benchmark_dir, "aggregation-protocol.R"))
source(file.path(benchmark_dir, "aggregation-synthetic.R"))
source(file.path(benchmark_dir, "aggregation-synthetic-summary.R"))

aggregation_synthetic_usage <- function() {
    cat(paste0(
        "Usage: Rscript benchmark/run-aggregation-synthetic.R [options]\n",
        "  --workers=N     Fork workers on Unix (default: up to 4)\n",
        "  --chunk-size=N  Latent pairs per checkpoint (default: 128)\n",
        "  --output=DIR    Result/checkpoint directory\n",
        "  --help          Show this text\n"
    ))
}

aggregation_positive_integer <- function(value, label) {
    parsed <- suppressWarnings(as.integer(value))
    if (length(parsed) != 1L || is.na(parsed) || parsed < 1L ||
        !identical(as.character(parsed), value)) {
        stop(label, " must be a positive integer.", call. = FALSE)
    }
    parsed
}

aggregation_synthetic_options <- function(args) {
    detected <- parallel::detectCores(logical = TRUE)
    default_workers <- if (is.finite(detected)) min(4L, detected) else 1L
    options <- list(
        workers = as.integer(default_workers),
        chunk_size = 128L,
        output = NULL
    )
    for (argument in args) {
        if (identical(argument, "--help")) {
            aggregation_synthetic_usage()
            quit(save = "no", status = 0L)
        } else if (startsWith(argument, "--workers=")) {
            options$workers <- aggregation_positive_integer(
                sub("^--workers=", "", argument),
                "`--workers`"
            )
        } else if (startsWith(argument, "--chunk-size=")) {
            options$chunk_size <- aggregation_positive_integer(
                sub("^--chunk-size=", "", argument),
                "`--chunk-size`"
            )
        } else if (startsWith(argument, "--output=")) {
            options$output <- sub("^--output=", "", argument)
        } else {
            stop("Unknown synthetic runner option: ", argument, call. = FALSE)
        }
    }
    if (.Platform$OS.type == "windows" && options$workers != 1L) {
        stop(
            "Windows synthetic execution requires `--workers=1`.",
            call. = FALSE
        )
    }
    options
}

aggregation_clean_git_head <- function(root) {
    head <- benchmark_command_output("git", c("rev-parse", "HEAD"), root)
    status <- benchmark_command_output(
        "git",
        c("status", "--porcelain", "--untracked-files=normal"),
        root
    )
    if (length(head) != 1L || !grepl("^[0-9a-f]{40}$", head[[1L]])) {
        stop("Cannot resolve the synthetic runner Git commit.", call. = FALSE)
    }
    if (length(status) > 0L) {
        stop(
            "Full synthetic execution requires a clean committed tree.",
            call. = FALSE
        )
    }
    head[[1L]]
}

aggregation_publish_locked_file <- function(source, target, expected_hash) {
    if (file.exists(target)) {
        if (!identical(unname(tools::sha256sum(target)), expected_hash)) {
            stop("Existing synthetic protocol artifact differs.", call. = FALSE)
        }
        return(invisible(target))
    }
    if (!file.copy(source, target, overwrite = FALSE) ||
        !identical(unname(tools::sha256sum(target)), expected_hash)) {
        stop("Cannot publish the synthetic protocol artifact.", call. = FALSE)
    }
    invisible(target)
}

aggregation_publish_locked_tsv <- function(value, target) {
    temporary <- tempfile("aggregation-manifest-", tmpdir = dirname(target))
    on.exit(unlink(temporary, force = TRUE), add = TRUE)
    benchmark_write_tsv(value, temporary)
    if (file.exists(target)) {
        if (!identical(
            unname(tools::sha256sum(target)),
            unname(tools::sha256sum(temporary))
        )) {
            stop("Existing synthetic manifest differs.", call. = FALSE)
        }
    } else if (!file.rename(temporary, target)) {
        stop("Cannot publish the synthetic manifest.", call. = FALSE)
    }
    invisible(target)
}

aggregation_install_package <- function(root, output_dir) {
    library <- tempfile("aggregation-library-")
    dir.create(library)
    published <- FALSE
    on.exit({
        if (!published) {
            unlink(library, recursive = TRUE, force = TRUE)
        }
    }, add = TRUE)
    log <- file.path(output_dir, "installation.log")
    status <- system2(
        file.path(R.home("bin"), "R"),
        c(
            "CMD", "INSTALL", "--preclean", "--clean",
            shQuote(paste0("--library=", library)),
            shQuote(root)
        ),
        stdout = log,
        stderr = log
    )
    if (!identical(status, 0L)) {
        stop("Committed package installation failed; see installation.log.",
            call. = FALSE)
    }
    .libPaths(c(library, .libPaths()))
    namespace <- loadNamespace("genefunnel", lib.loc = library)
    result <- list(
        library = library,
        namespace = namespace,
        audit = get(".aggregation_audit", envir = namespace, inherits = FALSE)
    )
    published <- TRUE
    result
}

aggregation_synthetic_metadata <- function(
    root,
    head,
    package,
    options,
    observations,
    cross_validation
) {
    metadata <- benchmark_metadata(
        repo_root = root,
        runner = "benchmark/run-aggregation-synthetic.R",
        preset = "aggregation_synthetic",
        repeats = 1L,
        workers = options$workers
    )
    metadata$value[metadata$key == "protocol_version"] <-
        AGGREGATION_PROTOCOL_VERSION
    additions <- data.frame(
        key = c(
            "aggregation_protocol_sha256", "aggregation_data_sha256",
            "committed_git_head", "chunk_size", "rng_kind",
            "measurement_rows", "latent_scenarios", "eligible_measurements",
            "model_excluded_pairs", "installed_package_path"
        ),
        value = c(
            AGGREGATION_PROTOCOL_SHA256,
            AGGREGATION_DATA_SHA256,
            head,
            options$chunk_size,
            "Mersenne-Twister/Inversion/Rejection",
            nrow(observations),
            length(unique(observations$latent_id)),
            sum(observations$eligible),
            cross_validation$excluded_pairs,
            normalizePath(file.path(package$library, "genefunnel"))
        ),
        stringsAsFactors = FALSE,
        check.names = FALSE
    )
    rbind(metadata, additions)
}

aggregation_artifact_manifest <- function(output_dir, files) {
    paths <- file.path(output_dir, files)
    if (any(!file.exists(paths))) {
        stop("Synthetic evidence artifacts are incomplete.", call. = FALSE)
    }
    data.frame(
        artifact = files,
        bytes = as.numeric(file.info(paths)$size),
        sha256 = unname(tools::sha256sum(paths)),
        stringsAsFactors = FALSE,
        check.names = FALSE
    )
}

options <- aggregation_synthetic_options(commandArgs(trailingOnly = TRUE))
if (is.null(options$output)) {
    stamp <- format(Sys.time(), "%Y%m%d-%H%M%S", tz = "UTC")
    options$output <- file.path(
        benchmark_dir,
        "results",
        paste(stamp, "aggregation-synthetic", sep = "-")
    )
}
dir.create(options$output, recursive = TRUE, showWarnings = FALSE)
output_dir <- normalizePath(options$output, mustWork = TRUE)
head <- aggregation_clean_git_head(repo_root)

registry_path <- file.path(benchmark_dir, "aggregation-protocol.tsv")
registry <- aggregation_read_registry(registry_path)
design <- aggregation_synthetic_design(registry)
aggregation_validate_endpoint_contract(registry)
aggregation_publish_locked_file(
    registry_path,
    file.path(output_dir, "protocol.tsv"),
    AGGREGATION_PROTOCOL_SHA256
)
aggregation_publish_locked_tsv(design, file.path(output_dir, "manifest.tsv"))

package <- NULL
tryCatch({
    cat("Installing committed package snapshot...\n")
    package <- aggregation_install_package(repo_root, output_dir)
    common_scale <- as.numeric(aggregation_registry_value(
        registry,
        "synthetic",
        "common_scale"
    ))
    smoke_design <- aggregation_smoke_design(design)
    smoke <- aggregation_simulate_design(
        smoke_design,
        package$audit,
        common_scale,
        workers = 1L,
        chunk_size = 2L
    )
    aggregation_validate_smoke(smoke, smoke_design)

    cat(
        "Simulating ", nrow(design), " measurements with ", options$workers,
        " worker(s)...\n",
        sep = ""
    )
    checkpoint_dir <- file.path(output_dir, "checkpoints")
    fingerprint <- paste(AGGREGATION_PROTOCOL_SHA256, head, sep = "::")
    observations <- aggregation_simulate_design(
        design,
        package$audit,
        common_scale,
        workers = options$workers,
        chunk_size = options$chunk_size,
        checkpoint_dir = checkpoint_dir,
        fingerprint = fingerprint
    )
    benchmark_write_tsv(observations, file.path(output_dir, "observations.tsv"))

    cat("Fitting frozen paired cross-validation models...\n")
    cross_validation <- aggregation_cross_validate(observations, registry)
    bootstrap <- aggregation_bootstrap_reductions(
        cross_validation$predictions,
        registry
    )
    endpoints <- aggregation_synthetic_endpoints(
        observations,
        cross_validation,
        bootstrap,
        registry
    )
    strata <- aggregation_stratified_errors(observations, registry)
    summary <- aggregation_synthetic_summary(
        observations,
        endpoints,
        cross_validation,
        bootstrap
    )

    metadata <- aggregation_synthetic_metadata(
        repo_root,
        head,
        package,
        options,
        observations,
        cross_validation
    )
    benchmark_write_tsv(endpoints, file.path(output_dir, "endpoints.tsv"))
    benchmark_write_tsv(
        cross_validation$folds,
        file.path(output_dir, "fold-results.tsv")
    )
    benchmark_write_tsv(
        cross_validation$predictions,
        file.path(output_dir, "predictions.tsv")
    )
    benchmark_write_tsv(
        cross_validation$coefficients,
        file.path(output_dir, "model-coefficients.tsv")
    )
    benchmark_write_tsv(bootstrap, file.path(output_dir, "bootstrap.tsv"))
    benchmark_write_tsv(strata, file.path(output_dir, "strata.tsv"))
    benchmark_write_tsv(summary, file.path(output_dir, "summary.tsv"))
    benchmark_write_tsv(metadata, file.path(output_dir, "metadata.tsv"))
    benchmark_write_session_info(
        file.path(output_dir, "session-info.txt"),
        c("genefunnel", "Matrix", "BiocParallel", "Rcpp", "RcppArmadillo")
    )
    aggregation_write_synthetic_report(
        file.path(output_dir, "report.md"),
        endpoints,
        summary,
        strata,
        metadata,
        options$workers
    )

    artifact_files <- c(
        "protocol.tsv", "manifest.tsv", "installation.log", "observations.tsv",
        "endpoints.tsv", "fold-results.tsv", "predictions.tsv",
        "model-coefficients.tsv", "bootstrap.tsv", "strata.tsv", "summary.tsv",
        "metadata.tsv", "session-info.txt", "report.md"
    )
    artifacts <- aggregation_artifact_manifest(output_dir, artifact_files)
    benchmark_write_tsv(artifacts, file.path(output_dir, "artifacts.tsv"))

    outcome <- if (summary$all_synthetic_gates_pass) "PASS" else "FAIL"
    cat(
        "Aggregation synthetic validation complete (", outcome, "): ",
        output_dir, "\n",
        sep = ""
    )
}, finally = {
    if (!is.null(package)) {
        unlink(package$library, recursive = TRUE, force = TRUE)
    }
})
