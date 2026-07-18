# Assisted-by: OpenAI Codex.

if (!nzchar(Sys.getenv("TZ", unset = ""))) {
    Sys.setenv(TZ = "UTC")
}

runner_file <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(runner_file) != 1L) {
    stop("Cannot locate run-aggregation-cellbench.R.", call. = FALSE)
}
benchmark_dir <- dirname(normalizePath(sub("^--file=", "", runner_file)))
repo_root <- dirname(benchmark_dir)
source(file.path(benchmark_dir, "protocol.R"))
source(file.path(benchmark_dir, "provenance.R"))
source(file.path(benchmark_dir, "report.R"))
source(file.path(benchmark_dir, "aggregation-protocol.R"))
source(file.path(benchmark_dir, "aggregation-runner.R"))
source(file.path(benchmark_dir, "aggregation-cellbench.R"))

aggregation_cellbench_usage <- function() {
    cat(paste0(
        "Usage: Rscript benchmark/run-aggregation-cellbench.R [options]\n",
        "  --data-dir=DIR  Directory containing four pinned CellBench files\n",
        "  --output=DIR    Result directory\n",
        "  --help          Show this text\n"
    ))
}

aggregation_cellbench_options <- function(args) {
    options <- list(data_dir = NULL, output = NULL)
    for (argument in args) {
        if (identical(argument, "--help")) {
            aggregation_cellbench_usage()
            quit(save = "no", status = 0L)
        } else if (startsWith(argument, "--data-dir=")) {
            options$data_dir <- sub("^--data-dir=", "", argument)
        } else if (startsWith(argument, "--output=")) {
            options$output <- sub("^--output=", "", argument)
        } else {
            stop("Unknown CellBench runner option: ", argument, call. = FALSE)
        }
    }
    if (is.null(options$data_dir) || !nzchar(options$data_dir) ||
        !dir.exists(options$data_dir)) {
        stop("`--data-dir` must name an existing directory.", call. = FALSE)
    }
    options$data_dir <- normalizePath(options$data_dir)
    options
}

aggregation_cellbench_verified_files <- function(data_dir, manifest) {
    ids <- c(
        "cellbench_celseq2_counts", "cellbench_celseq2_metadata",
        "cellbench_sortseq_counts", "cellbench_sortseq_metadata"
    )
    selected <- manifest[match(ids, manifest$data_id), , drop = FALSE]
    files <- aggregation_external_file_map(data_dir)[ids]
    aggregation_verify_data_files(selected, files)
}

aggregation_cellbench_read_inputs <- function(files, registry) {
    paths <- stats::setNames(files$path, files$data_id)
    labels <- c(
        aggregation_registry_value(registry, "cellbench", "training_platform"),
        aggregation_registry_value(registry, "cellbench", "heldout_platform")
    )
    counts <- c("cellbench_celseq2_counts", "cellbench_sortseq_counts")
    metadata <- c("cellbench_celseq2_metadata", "cellbench_sortseq_metadata")
    result <- lapply(seq_along(labels), function(index) {
        aggregation_cellbench_read_platform(
            paths[[counts[[index]]]], paths[[metadata[[index]]]],
            labels[[index]], registry
        )
    })
    stats::setNames(result, labels)
}

aggregation_cellbench_metadata <- function(
    root,
    head,
    package,
    analysis,
    files
) {
    metadata <- benchmark_metadata(
        repo_root = root, runner = "benchmark/run-aggregation-cellbench.R",
        preset = "aggregation_cellbench", repeats = 1L, workers = 1L
    )
    additions <- data.frame(
        key = c(
            "aggregation_protocol_sha256", "aggregation_data_sha256",
            "committed_git_head", "input_files", "input_bytes",
            "common_genes", "selected_sets", "mixed_library_set_rows",
            "installed_package_path"
        ),
        value = c(
            AGGREGATION_PROTOCOL_SHA256, AGGREGATION_DATA_SHA256, head,
            nrow(files), sum(files$bytes), analysis$summary$common_genes,
            analysis$summary$selected_sets,
            analysis$summary$mixed_library_set_rows,
            normalizePath(file.path(package$library, "genefunnel"))
        ),
        stringsAsFactors = FALSE, check.names = FALSE
    )
    metadata$value[metadata$key == "protocol_version"] <-
        AGGREGATION_PROTOCOL_VERSION
    rbind(metadata, additions)
}

aggregation_cellbench_report_curves <- function(curves) {
    selected <- curves[!curves$passed, , drop = FALSE]
    if (!nrow(selected)) {
        selected <- curves[order(
            -curves$quantile90_absolute_error,
            -curves$median_absolute_error
        ), , drop = FALSE]
    }
    fields <- c(
        "platform", "composition", "mRNA_amount", "set", "library_count",
        "median_absolute_error", "quantile90_absolute_error", "passed"
    )
    utils::head(selected[fields], 20L)
}

aggregation_write_cellbench_report <- function(
    path,
    analysis,
    metadata,
    data_dir
) {
    passed <- isTRUE(analysis$summary$all_cellbench_gates_pass)
    lines <- c(
        "# GeneFunnel CellBench mixture validation", "",
        paste0(
            "Protocol `", AGGREGATION_PROTOCOL_VERSION, "`. ",
            if (passed) "PASS" else "FAIL",
            " - frozen CellBench co-primary decision."
        ), "",
        paste0(
            "This process-control result uses known RNA mixtures. It is not ",
            "a biological effect or public-API validation."
        ), "", "## Endpoints", "",
        benchmark_markdown_table(analysis$endpoints), "",
        "## Failed groups or 20 worst passing groups", "",
        benchmark_markdown_table(aggregation_cellbench_report_curves(
            analysis$curve_groups
        )), "", "## Environment", "",
        benchmark_markdown_table(metadata), "", "## Reproduce", "",
        "```sh", paste0(
            "R_LIBS_USER=\"$PWD/.agent/R-library\" Rscript --vanilla ",
            "benchmark/run-aggregation-cellbench.R --data-dir=",
            shQuote(data_dir)
        ), "```", "",
        "All group rows, including failures, remain in `curve-groups.tsv`."
    )
    writeLines(lines, path, useBytes = TRUE)
}

options <- aggregation_cellbench_options(commandArgs(trailingOnly = TRUE))
if (is.null(options$output)) {
    stamp <- format(Sys.time(), "%Y%m%d-%H%M%S", tz = "UTC")
    options$output <- file.path(
        benchmark_dir,
        "results",
        paste(stamp, "aggregation-cellbench", sep = "-")
    )
}
dir.create(options$output, recursive = TRUE, showWarnings = FALSE)
output_dir <- normalizePath(options$output, mustWork = TRUE)
head <- aggregation_clean_git_head(repo_root)

registry_path <- file.path(benchmark_dir, "aggregation-protocol.tsv")
data_path <- file.path(benchmark_dir, "aggregation-data.tsv")
registry <- aggregation_read_registry(registry_path)
manifest <- aggregation_read_data_manifest(data_path)
files <- aggregation_cellbench_verified_files(options$data_dir, manifest)
aggregation_publish_locked_file(
    registry_path, file.path(output_dir, "protocol.tsv"),
    AGGREGATION_PROTOCOL_SHA256
)
aggregation_publish_locked_file(
    data_path, file.path(output_dir, "data-manifest.tsv"),
    AGGREGATION_DATA_SHA256
)
benchmark_write_tsv(files, file.path(output_dir, "data-files.tsv"))

package <- NULL
old_locale <- Sys.getlocale("LC_COLLATE")
tryCatch({
    if (is.na(suppressWarnings(Sys.setlocale("LC_COLLATE", "C")))) {
        stop("Cannot establish frozen C collation.", call. = FALSE)
    }
    cat("Installing committed package snapshot...\n")
    package <- aggregation_install_package(repo_root, output_dir)
    cat("Reading and validating pinned CellBench inputs...\n")
    platforms <- aggregation_cellbench_read_inputs(files, registry)
    cat("Selecting training-only sets and evaluating mixtures...\n")
    analysis <- aggregation_cellbench_analyze(
        platforms, registry, package$audit, package$scorer
    )
    metadata <- aggregation_cellbench_metadata(
        repo_root, head, package, analysis, files
    )

    benchmark_write_tsv(
        analysis$selection$manifest,
        file.path(output_dir, "set-manifest.tsv")
    )
    benchmark_write_tsv(
        analysis$selection$membership,
        file.path(output_dir, "set-membership.tsv")
    )
    benchmark_write_tsv(
        aggregation_cellbench_profile_table(analysis$profiles),
        file.path(output_dir, "pure-profiles.tsv")
    )
    benchmark_write_tsv(
        analysis$references,
        file.path(output_dir, "references.tsv")
    )
    benchmark_write_tsv(
        analysis$observations,
        file.path(output_dir, "observations.tsv")
    )
    benchmark_write_tsv(
        analysis$curve_groups,
        file.path(output_dir, "curve-groups.tsv")
    )
    benchmark_write_tsv(
        analysis$conditions,
        file.path(output_dir, "condition-medians.tsv")
    )
    benchmark_write_tsv(analysis$endpoints, file.path(output_dir, "endpoints.tsv"))
    benchmark_write_tsv(analysis$summary, file.path(output_dir, "summary.tsv"))
    benchmark_write_tsv(metadata, file.path(output_dir, "metadata.tsv"))
    benchmark_write_session_info(
        file.path(output_dir, "session-info.txt"),
        c("genefunnel", "Matrix", "BiocParallel", "Rcpp", "RcppArmadillo")
    )
    aggregation_write_cellbench_report(
        file.path(output_dir, "report.md"), analysis, metadata, options$data_dir
    )
    artifact_files <- c(
        "protocol.tsv", "data-manifest.tsv", "data-files.tsv",
        "installation.log", "set-manifest.tsv", "set-membership.tsv",
        "pure-profiles.tsv", "references.tsv", "observations.tsv",
        "curve-groups.tsv", "condition-medians.tsv", "endpoints.tsv",
        "summary.tsv", "metadata.tsv", "session-info.txt", "report.md"
    )
    artifacts <- aggregation_artifact_manifest(output_dir, artifact_files)
    benchmark_write_tsv(artifacts, file.path(output_dir, "artifacts.tsv"))
    outcome <- if (analysis$summary$all_cellbench_gates_pass) "PASS" else "FAIL"
    cat(
        "CellBench aggregation validation complete (", outcome, "): ",
        output_dir, "\n", sep = ""
    )
}, finally = {
    suppressWarnings(Sys.setlocale("LC_COLLATE", old_locale))
    if (!is.null(package)) {
        unlink(package$library, recursive = TRUE, force = TRUE)
    }
})
