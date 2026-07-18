# Assisted-by: OpenAI Codex.

if (!nzchar(Sys.getenv("TZ", unset = ""))) {
    Sys.setenv(TZ = "UTC")
}

runner_file <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(runner_file) != 1L) {
    stop("Cannot locate run-aggregation-kang.R.", call. = FALSE)
}
benchmark_dir <- dirname(normalizePath(sub("^--file=", "", runner_file)))
repo_root <- dirname(benchmark_dir)
source(file.path(benchmark_dir, "protocol.R"))
source(file.path(benchmark_dir, "provenance.R"))
source(file.path(benchmark_dir, "report.R"))
source(file.path(benchmark_dir, "aggregation-protocol.R"))
source(file.path(benchmark_dir, "aggregation-runner.R"))
source(file.path(benchmark_dir, "aggregation-kang.R"))
source(file.path(benchmark_dir, "aggregation-kang-summary.R"))

aggregation_kang_usage <- function() {
    cat(paste0(
        "Usage: Rscript benchmark/run-aggregation-kang.R [options]\n",
        "  --data-dir=DIR  Directory containing four pinned Kang/Reactome files\n",
        "  --output=DIR    Result directory\n",
        "  --help          Show this text\n"
    ))
}

aggregation_kang_options <- function(args) {
    options <- list(data_dir = NULL, output = NULL)
    for (argument in args) {
        if (identical(argument, "--help")) {
            aggregation_kang_usage()
            quit(save = "no", status = 0L)
        } else if (startsWith(argument, "--data-dir=")) {
            options$data_dir <- sub("^--data-dir=", "", argument)
        } else if (startsWith(argument, "--output=")) {
            options$output <- sub("^--output=", "", argument)
        } else {
            stop("Unknown Kang runner option: ", argument, call. = FALSE)
        }
    }
    if (is.null(options$data_dir) || !nzchar(options$data_dir) ||
        !dir.exists(options$data_dir)) {
        stop("`--data-dir` must name an existing directory.", call. = FALSE)
    }
    options$data_dir <- normalizePath(options$data_dir)
    options
}

aggregation_kang_verified_files <- function(data_dir, manifest) {
    ids <- c(
        "kang_batch2_counts", "kang_batch2_genes",
        "kang_batch2_metadata", "reactome_v97"
    )
    selected <- manifest[match(ids, manifest$data_id), , drop = FALSE]
    files <- aggregation_external_file_map(data_dir)[ids]
    aggregation_verify_data_files(selected, files)
}

aggregation_kang_preprocessing_facts <- function(inputs, analysis) {
    collapsed <- analysis$collapsed$facts
    data.frame(
        fact = c(
            "raw_barcodes", "raw_duplicate_barcodes", "metadata_cells",
            "gene_rows", "nonempty_gene_rows", "unique_symbols",
            "duplicate_symbol_rows", "aggregate_columns", "collapsed_nonzero",
            "catalogue_pathways", "retained_pathways"
        ),
        value = c(
            sum(inputs$barcodes$batch_sizes), inputs$barcodes$raw_duplicates,
            nrow(inputs$metadata), collapsed$gene_rows,
            collapsed$nonempty_gene_rows, collapsed$unique_symbols,
            collapsed$duplicate_symbol_rows, collapsed$aggregate_columns,
            collapsed$collapsed_nonzero, nrow(inputs$catalogue$manifest),
            length(inputs$catalogue$sets)
        ),
        stringsAsFactors = FALSE, check.names = FALSE
    )
}

aggregation_kang_assignment_table <- function(inputs, analysis) {
    value <- analysis$units$assignments
    data.frame(
        cell_id = rownames(inputs$metadata)[value$cell_index],
        unit_id = value$unit_id, half = value$half,
        stringsAsFactors = FALSE, check.names = FALSE
    )
}

aggregation_kang_metadata <- function(root, head, package, analysis, files) {
    metadata <- benchmark_metadata(
        repo_root = root, runner = "benchmark/run-aggregation-kang.R",
        preset = "aggregation_kang", repeats = 1L, workers = 1L
    )
    additions <- data.frame(
        key = c(
            "aggregation_protocol_sha256", "aggregation_data_sha256",
            "committed_git_head", "input_files", "input_bytes",
            "fixed_units", "eligible_units", "retained_cells",
            "retained_pathways", "audit_rows", "installed_package_path"
        ),
        value = c(
            AGGREGATION_PROTOCOL_SHA256, AGGREGATION_DATA_SHA256, head,
            nrow(files), sum(files$bytes), analysis$summary$fixed_units,
            analysis$summary$eligible_units, analysis$summary$retained_cells,
            analysis$summary$retained_pathways, analysis$summary$audit_rows,
            normalizePath(file.path(package$library, "genefunnel"))
        ),
        stringsAsFactors = FALSE, check.names = FALSE
    )
    metadata$value[metadata$key == "protocol_version"] <-
        AGGREGATION_PROTOCOL_VERSION
    rbind(metadata, additions)
}

aggregation_kang_report_stability <- function(stability) {
    stability[order(
        stability$defined, stability$split_spearman,
        stability$donor, stability$condition, na.last = FALSE,
        method = "radix"
    ), , drop = FALSE]
}

aggregation_write_kang_report <- function(
    path,
    analysis,
    metadata,
    data_dir
) {
    status <- function(value) if (isTRUE(value)) "PASS" else "FAIL"
    summary <- analysis$summary
    lines <- c(
        "# GeneFunnel Kang/Reactome characterization", "",
        paste0("Protocol `", AGGREGATION_PROTOCOL_VERSION, "`."), "",
        paste0(
            "Technical stability: ", status(summary$technical_stability_pass),
            "; held-out replication: ", status(summary$heldout_replication_pass),
            "; biological sign test: ", status(summary$biological_effect_pass), "."
        ), "",
        paste0(
            "These donor-level results characterize the frozen perturbation ",
            "design. They cannot override the separate known-mixture gate."
        ), "", "## Endpoints", "",
        benchmark_markdown_table(analysis$endpoints), "",
        "## Primary-pathway decisions", "",
        benchmark_markdown_table(analysis$decisions), "",
        "## Split stability", "",
        benchmark_markdown_table(aggregation_kang_report_stability(
            analysis$stability
        )), "", "## Environment", "",
        benchmark_markdown_table(metadata), "", "## Reproduce", "",
        "```sh", paste0(
            "R_LIBS_USER=\"$PWD/.agent/R-library\" Rscript --vanilla ",
            "benchmark/run-aggregation-kang.R --data-dir=", shQuote(data_dir)
        ), "```", "",
        "All pathway/view rows remain in `audit-summary.tsv`."
    )
    writeLines(lines, path, useBytes = TRUE)
}

aggregation_kang_write_analysis <- function(
    output_dir,
    inputs,
    analysis,
    metadata
) {
    benchmark_write_tsv(inputs$matrix_facts, file.path(output_dir, "matrix-facts.tsv"))
    benchmark_write_tsv(
        aggregation_kang_preprocessing_facts(inputs, analysis),
        file.path(output_dir, "preprocessing-facts.tsv")
    )
    benchmark_write_tsv(analysis$units$manifest, file.path(output_dir, "unit-manifest.tsv"))
    benchmark_write_tsv(
        aggregation_kang_assignment_table(inputs, analysis),
        file.path(output_dir, "cell-assignments.tsv")
    )
    benchmark_write_tsv(inputs$catalogue$manifest, file.path(output_dir, "pathway-manifest.tsv"))
    benchmark_write_tsv(inputs$catalogue$membership, file.path(output_dir, "pathway-membership.tsv"))
    benchmark_write_tsv(analysis$audits, file.path(output_dir, "audit-summary.tsv"))
    benchmark_write_tsv(analysis$stability, file.path(output_dir, "stability.tsv"))
    benchmark_write_tsv(analysis$contrasts, file.path(output_dir, "donor-contrasts.tsv"))
    benchmark_write_tsv(analysis$decisions, file.path(output_dir, "pathway-decisions.tsv"))
    benchmark_write_tsv(analysis$endpoints, file.path(output_dir, "endpoints.tsv"))
    benchmark_write_tsv(analysis$summary, file.path(output_dir, "summary.tsv"))
    benchmark_write_tsv(metadata, file.path(output_dir, "metadata.tsv"))
}

options <- aggregation_kang_options(commandArgs(trailingOnly = TRUE))
if (is.null(options$output)) {
    stamp <- format(Sys.time(), "%Y%m%d-%H%M%S", tz = "UTC")
    options$output <- file.path(
        benchmark_dir, "results", paste(stamp, "aggregation-kang", sep = "-")
    )
}
dir.create(options$output, recursive = TRUE, showWarnings = FALSE)
output_dir <- normalizePath(options$output, mustWork = TRUE)
head <- aggregation_clean_git_head(repo_root)

registry_path <- file.path(benchmark_dir, "aggregation-protocol.tsv")
data_path <- file.path(benchmark_dir, "aggregation-data.tsv")
registry <- aggregation_read_registry(registry_path)
manifest <- aggregation_read_data_manifest(data_path)
files <- aggregation_kang_verified_files(options$data_dir, manifest)
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
extraction_dir <- tempfile("aggregation-kang-extract-", tmpdir = output_dir)
dir.create(extraction_dir)
old_locale <- Sys.getlocale("LC_COLLATE")
tryCatch({
    if (is.na(suppressWarnings(Sys.setlocale("LC_COLLATE", "C")))) {
        stop("Cannot establish frozen C collation.", call. = FALSE)
    }
    cat("Installing committed package snapshot...\n")
    package <- aggregation_install_package(repo_root, output_dir)
    cat("Reading and validating pinned Kang/Reactome inputs...\n")
    inputs <- aggregation_kang_external_inputs(files, registry, extraction_dir)
    cat("Aggregating cells and evaluating frozen donor endpoints...\n")
    analysis <- aggregation_kang_analyze(inputs, registry, package$audit)
    metadata <- aggregation_kang_metadata(
        repo_root, head, package, analysis, files
    )
    aggregation_kang_write_analysis(output_dir, inputs, analysis, metadata)
    benchmark_write_session_info(
        file.path(output_dir, "session-info.txt"),
        c("genefunnel", "Matrix", "BiocParallel", "Rcpp", "RcppArmadillo")
    )
    aggregation_write_kang_report(
        file.path(output_dir, "report.md"), analysis, metadata, options$data_dir
    )
    artifact_files <- c(
        "protocol.tsv", "data-manifest.tsv", "data-files.tsv",
        "installation.log", "matrix-facts.tsv", "preprocessing-facts.tsv",
        "unit-manifest.tsv", "cell-assignments.tsv", "pathway-manifest.tsv",
        "pathway-membership.tsv", "audit-summary.tsv", "stability.tsv",
        "donor-contrasts.tsv", "pathway-decisions.tsv", "endpoints.tsv",
        "summary.tsv", "metadata.tsv", "session-info.txt", "report.md"
    )
    artifacts <- aggregation_artifact_manifest(output_dir, artifact_files)
    benchmark_write_tsv(artifacts, file.path(output_dir, "artifacts.tsv"))
    cat(
        "Kang aggregation characterization complete (technical ",
        if (analysis$summary$technical_stability_pass) "PASS" else "FAIL",
        ", held-out ",
        if (analysis$summary$heldout_replication_pass) "PASS" else "FAIL",
        ", biological ",
        if (analysis$summary$biological_effect_pass) "PASS" else "FAIL",
        "): ", output_dir, "\n", sep = ""
    )
}, finally = {
    suppressWarnings(Sys.setlocale("LC_COLLATE", old_locale))
    unlink(extraction_dir, recursive = TRUE, force = TRUE)
    if (!is.null(package)) {
        unlink(package$library, recursive = TRUE, force = TRUE)
    }
})
