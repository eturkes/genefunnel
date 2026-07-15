# Assisted-by: OpenAI Codex.

benchmark_markdown_text <- function(value) {
    value <- ifelse(is.na(value), "NA", as.character(value))
    value <- gsub("\\|", "\\\\|", value)
    gsub("[\r\n]+", " ", value)
}

benchmark_markdown_table <- function(value) {
    value[] <- lapply(value, benchmark_markdown_text)
    header <- paste0("| ", paste(names(value), collapse = " | "), " |")
    rule <- paste0("| ", paste(rep("---", ncol(value)), collapse = " | "), " |")
    rows <- apply(value, 1L, function(row) {
        paste0("| ", paste(row, collapse = " | "), " |")
    })
    c(header, rule, rows)
}

benchmark_metadata_rows <- function(metadata) {
    keep <- c(
        "generated_utc", "protocol_version", "git_head", "git_dirty",
        "sysname", "release", "machine", "cpu_model", "logical_cores",
        "memory_total", "r_version", "genefunnel_version", "Matrix_version",
        "BiocParallel_version"
    )
    selected <- metadata[match(keep, metadata$key), , drop = FALSE]
    names(selected) <- c("Environment", "Recorded value")
    selected
}

benchmark_format_number <- function(value, digits = 6L) {
    ifelse(
        is.na(value),
        "NA",
        formatC(value, digits = digits, format = "fg", flag = "#")
    )
}

benchmark_write_performance_report <- function(
    path,
    protocol,
    summary,
    metadata,
    options
) {
    view <- data.frame(
        Scenario = summary$scenario_id,
        Storage = summary$storage,
        Overlap = summary$overlap,
        Backend = summary$backend,
        Dimensions = paste0(
            summary$n_features, " x ", summary$n_samples,
            "; ", summary$n_sets, " sets x ", summary$set_size
        ),
        Repeats = summary$repeats,
        `Median seconds` = benchmark_format_number(summary$elapsed_median_sec),
        `Score cells/second` = benchmark_format_number(
            summary$score_cells_per_sec_median,
            digits = 7L
        ),
        `Peak scoring increment MiB` = benchmark_format_number(
            summary$score_peak_increment_kib_max / 1024
        ),
        Digest = summary$output_md5,
        check.names = FALSE
    )
    source_row <- benchmark_protocol_subset(protocol, "performance")
    source_row <- source_row[source_row$scenario_id == options$preset, , drop = FALSE]
    command <- paste0(
        "R_LIBS_USER=\"$PWD/.agent/R-library\" Rscript --vanilla ",
        "benchmark/run.R --preset=", options$preset,
        " --repeats=", options$repeats,
        " --workers=", options$workers
    )
    lines <- c(
        "# GeneFunnel performance benchmark report",
        "",
        paste0("Protocol: `", GENEFUNNEL_BENCHMARK_PROTOCOL, "`; preset: `", options$preset, "`."),
        "",
        paste0(
            "This report records synthetic computational evidence. It is not a ",
            "reproduction of thesis datasets or historical runtimes. Timings and ",
            "memory observations are specific to the recorded hardware and software ",
            "environment; no row is a universal performance threshold."
        ),
        "",
        "## Protocol",
        "",
        paste0("- Matrix construction: ", source_row$matrix_construction),
        paste0("- Gene-set construction: ", source_row$set_construction),
        paste0("- Preprocessing: ", source_row$preprocessing),
        paste0("- Methods: ", source_row$methods),
        paste0("- Assertion: ", source_row$assertion),
        "",
        paste0(
            "The expanded `manifest.tsv` records every scenario and seed. `runs.tsv` ",
            "retains isolated-run metrics and package versions; `metadata.tsv` and ",
            "`session-info.txt` retain the execution environment."
        ),
        "",
        "## Results",
        "",
        benchmark_markdown_table(view),
        "",
        "## Environment",
        "",
        benchmark_markdown_table(benchmark_metadata_rows(metadata)),
        "",
        "## Reproduce",
        "",
        "```sh",
        command,
        "```",
        ""
    )
    writeLines(lines, path)
}

benchmark_write_controlled_report <- function(
    path,
    manifest,
    results,
    summary,
    metadata
) {
    protocol_view <- data.frame(
        Scenario = manifest$scenario_id,
        Seeds = paste0(manifest$matrix_seed, "; ", manifest$set_seed),
        Dimensions = manifest$dimensions,
        `Matrix construction` = manifest$matrix_construction,
        `Set construction` = manifest$set_construction,
        Methods = manifest$methods,
        Assertion = manifest$assertion,
        check.names = FALSE
    )
    results_view <- data.frame(
        Scenario = results$scenario_id,
        Assertion = results$assertion_id,
        Expected = results$expected,
        Observed = results$observed,
        Tolerance = ifelse(
            is.na(results$tolerance),
            "NA",
            formatC(results$tolerance, digits = 3L, format = "e")
        ),
        Passed = ifelse(results$passed, "yes", "NO"),
        check.names = FALSE
    )
    lines <- c(
        "# GeneFunnel controlled scientific benchmark report",
        "",
        paste0("Protocol: `", GENEFUNNEL_BENCHMARK_PROTOCOL, "`."),
        "",
        paste0(
            "This is a deterministic synthetic protocol derived from the normative ",
            "GeneFunnel specification. It uses no thesis text, data, result, or ",
            "historical runtime and makes no claim that a thesis analysis was ",
            "reproduced."
        ),
        "",
        paste0(
            "Outcome: ", summary$passed_assertions, "/", summary$total_assertions,
            " assertions passed across ", summary$total_scenarios, " scenarios."
        ),
        "",
        "## Protocol",
        "",
        benchmark_markdown_table(protocol_view),
        "",
        paste0(
            "All values are constructed directly on a non-negative scale with a ",
            "meaningful zero; there is no normalization, imputation, identifier ",
            "conversion, or positive shifting. Analytic cases have no random seed. ",
            "The sparse case records its matrix and set seeds above."
        ),
        "",
        "## Results",
        "",
        benchmark_markdown_table(results_view),
        "",
        "## Environment",
        "",
        benchmark_markdown_table(benchmark_metadata_rows(metadata)),
        "",
        "## Reproduce",
        "",
        "```sh",
        "R_LIBS_USER=\"$PWD/.agent/R-library\" Rscript --vanilla benchmark/run-controlled.R",
        "```",
        ""
    )
    writeLines(lines, path)
}
