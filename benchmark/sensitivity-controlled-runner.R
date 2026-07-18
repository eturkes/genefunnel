# Assisted-by: OpenAI Codex.

sensitivity_controlled_git_state <- function(root, candidate_id, mode) {
    head <- benchmark_command_output("git", c("rev-parse", "HEAD"), root)
    status <- benchmark_command_output(
        "git", c("status", "--porcelain", "--untracked-files=normal"), root
    )
    if (length(head) != 1L || !grepl("^[0-9a-f]{40}$", head[[1L]])) {
        stop("Cannot resolve sensitivity controlled Git HEAD.", call. = FALSE)
    }
    candidate <- if (is.null(candidate_id)) head[[1L]] else candidate_id
    if (!grepl("^[0-9a-f]{40}$", candidate)) {
        stop("Sensitivity controlled candidate SHA is invalid.", call. = FALSE)
    }
    gate_valid <- !is.null(candidate_id) &&
        identical(candidate, head[[1L]]) && length(status) == 0L
    if (mode == "gate" && !gate_valid) {
        stop(
            "Gate mode requires the explicit clean current Git SHA.",
            call. = FALSE
        )
    }
    list(
        head = head[[1L]], candidate = candidate,
        dirty = length(status) > 0L
    )
}

sensitivity_controlled_output_path <- function(path, root, mode) {
    if (is.null(path)) {
        stamp <- format(Sys.time(), "%Y%m%d-%H%M%S", tz = "UTC")
        path <- file.path("benchmark", "results", paste(
            stamp, "sensitivity-controlled", mode, sep = "-"
        ))
    }
    root <- normalizePath(root, mustWork = TRUE)
    if (!startsWith(path, .Platform$file.sep)) path <- file.path(root, path)
    path <- normalizePath(path, mustWork = FALSE)
    if (!startsWith(path, paste0(root, .Platform$file.sep))) {
        stop("Sensitivity controlled output must be inside the repository.")
    }
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
    normalizePath(path, mustWork = TRUE)
}

sensitivity_controlled_publish_file <- function(source, target, expected_md5) {
    if (file.exists(target)) {
        if (!identical(unname(tools::md5sum(target)), expected_md5)) {
            stop("Existing sensitivity protocol artifact differs.", call. = FALSE)
        }
    } else if (!file.copy(source, target, overwrite = FALSE) ||
        !identical(unname(tools::md5sum(target)), expected_md5)) {
        stop("Cannot publish sensitivity protocol artifact.", call. = FALSE)
    }
    invisible(target)
}

sensitivity_controlled_publish_tsv <- function(value, target) {
    temporary <- tempfile("sensitivity-controlled-", tmpdir = dirname(target))
    on.exit(unlink(temporary, force = TRUE), add = TRUE)
    benchmark_write_tsv(value, temporary)
    if (file.exists(target)) {
        if (!identical(
            unname(tools::sha256sum(target)),
            unname(tools::sha256sum(temporary))
        )) stop("Existing sensitivity locked table differs.", call. = FALSE)
    } else if (!file.rename(temporary, target)) {
        stop("Cannot publish sensitivity locked table.", call. = FALSE)
    }
    invisible(target)
}

sensitivity_controlled_run_config <- function(
    git, mode, chunk_size, design
) {
    data.frame(
        key = c(
            "parent_protocol", "parent_registry_md5", "execution_protocol",
            "execution_protocol_md5", "candidate_id", "mode", "chunk_size",
            "scenario_count", "first_scenario_id", "last_scenario_id"
        ),
        value = c(
            SENSITIVITY_PROTOCOL_VERSION, SENSITIVITY_PROTOCOL_MD5,
            SENSITIVITY_CONTROLLED_VERSION, SENSITIVITY_CONTROLLED_MD5,
            git$candidate, mode, chunk_size, nrow(design),
            design$scenario_id[[1L]], design$scenario_id[[nrow(design)]]
        ),
        stringsAsFactors = FALSE, check.names = FALSE
    )
}

sensitivity_controlled_config_fingerprint <- function(config) {
    paste(config$key, config$value, sep = "=", collapse = "::")
}

sensitivity_controlled_git_archive <- function(root, output_dir, candidate) {
    target <- file.path(output_dir, "source.tar")
    temporary <- tempfile("sensitivity-source-", tmpdir = output_dir, fileext = ".tar")
    on.exit(unlink(temporary, force = TRUE), add = TRUE)
    status <- system2(
        "git",
        c(
            "-C", shQuote(root), "archive", "--format=tar",
            shQuote(paste0("--output=", temporary)), shQuote(candidate)
        ),
        stdout = FALSE, stderr = FALSE
    )
    if (!identical(status, 0L) || !file.exists(temporary)) {
        stop("Cannot archive sensitivity controlled candidate.", call. = FALSE)
    }
    digest <- unname(tools::sha256sum(temporary))
    if (file.exists(target)) {
        if (!identical(unname(tools::sha256sum(target)), digest)) {
            stop("Existing sensitivity source archive differs.", call. = FALSE)
        }
    } else if (!file.rename(temporary, target)) {
        stop("Cannot publish sensitivity source archive.", call. = FALSE)
    }
    list(
        path = target, sha256 = digest,
        md5 = unname(tools::md5sum(target))
    )
}

sensitivity_controlled_install <- function(root, output_dir, candidate) {
    archive <- sensitivity_controlled_git_archive(root, output_dir, candidate)
    source <- tempfile("sensitivity-source-", tmpdir = output_dir)
    library <- tempfile("sensitivity-library-", tmpdir = output_dir)
    dir.create(source)
    dir.create(library)
    success <- FALSE
    on.exit({
        unlink(source, recursive = TRUE, force = TRUE)
        if (!success) unlink(library, recursive = TRUE, force = TRUE)
    }, add = TRUE)
    utils::untar(archive$path, exdir = source)
    log <- file.path(output_dir, "installation.log")
    status <- system2(
        file.path(R.home("bin"), "R"),
        c(
            "CMD", "INSTALL", "--no-multiarch", "--with-keep.source",
            shQuote(paste0("--library=", library)), shQuote(source)
        ),
        stdout = log, stderr = log
    )
    package <- file.path(library, "genefunnel")
    if (!identical(status, 0L) || !dir.exists(package)) {
        stop("Cannot install sensitivity candidate; see installation.log.")
    }
    old_paths <- .libPaths()
    .libPaths(c(library, old_paths))
    namespace <- loadNamespace("genefunnel", lib.loc = library)
    success <- TRUE
    list(
        library = library, old_paths = old_paths, namespace = namespace,
        scorer = get("genefunnel", envir = namespace, inherits = FALSE),
        sensitivity = get(".sensitivity_cell", envir = namespace, inherits = FALSE),
        archive_sha256 = archive$sha256, archive_md5 = archive$md5,
        installed_manifest_md5 = components_manifest_digest(package)
    )
}

sensitivity_controlled_uninstall <- function(package) {
    if (is.null(package)) return(invisible(NULL))
    if ("genefunnel" %in% loadedNamespaces()) {
        try(unloadNamespace("genefunnel"), silent = TRUE)
    }
    .libPaths(package$old_paths)
    unlink(package$library, recursive = TRUE, force = TRUE)
    invisible(NULL)
}

sensitivity_controlled_chunk_path <- function(directory, chunk_id) {
    file.path(directory, sprintf("chunk-%04d.rds", chunk_id))
}

sensitivity_controlled_simulate_chunk <- function(
    chunk_id, indices, design, registry, scorer, sensitivity, backend,
    checkpoint_dir, fingerprint
) {
    expected_design <- design[indices, , drop = FALSE]
    expected_ids <- expected_design$scenario_id
    path <- sensitivity_controlled_chunk_path(checkpoint_dir, chunk_id)
    if (file.exists(path)) {
        cached <- readRDS(path)
        valid <- identical(attr(cached, "fingerprint"), fingerprint) &&
            identical(attr(cached, "scenario_ids"), expected_ids)
        if (!valid) stop("Sensitivity checkpoint identity mismatch.", call. = FALSE)
        sensitivity_controlled_validate_observations(
            cached, expected_design, registry
        )
        return(cached)
    }
    pieces <- lapply(indices, function(index) {
        sensitivity_controlled_observe_scenario(
            design[index, , drop = FALSE], registry, scorer, sensitivity, backend
        )
    })
    observed <- sensitivity_controlled_bind_observations(pieces)
    sensitivity_controlled_validate_observations(observed, expected_design, registry)
    attr(observed, "fingerprint") <- fingerprint
    attr(observed, "scenario_ids") <- expected_ids
    temporary <- paste0(path, ".tmp-", Sys.getpid())
    saveRDS(observed, temporary, compress = TRUE, version = 3L)
    if (!file.rename(temporary, path)) {
        unlink(temporary, force = TRUE)
        stop("Cannot publish sensitivity checkpoint.", call. = FALSE)
    }
    observed
}

sensitivity_controlled_execute_chunks <- function(sequence, run_chunk, workers) {
    if (workers == 1L || .Platform$OS.type == "windows") {
        return(lapply(sequence, run_chunk))
    }
    parallel::mclapply(
        sequence, run_chunk, mc.cores = workers,
        mc.preschedule = TRUE, mc.set.seed = FALSE
    )
}

sensitivity_controlled_simulate_design <- function(
    design, registry, scorer, sensitivity, backend, workers, chunk_size,
    checkpoint_dir, fingerprint
) {
    workers <- as.integer(workers)
    chunk_size <- as.integer(chunk_size)
    if (workers < 1L || chunk_size < 1L || nrow(design) < 1L) {
        stop("Sensitivity controlled execution shape is invalid.", call. = FALSE)
    }
    dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)
    chunks <- split(seq_len(nrow(design)), ceiling(seq_len(nrow(design)) / chunk_size))
    run_chunk <- function(chunk_id) {
        sensitivity_controlled_simulate_chunk(
            chunk_id, chunks[[chunk_id]], design, registry, scorer, sensitivity,
            backend, checkpoint_dir, fingerprint
        )
    }
    observed <- sensitivity_controlled_execute_chunks(
        seq_along(chunks), run_chunk, workers
    )
    if (any(vapply(observed, inherits, logical(1L), "try-error"))) {
        stop("One or more sensitivity chunks failed.", call. = FALSE)
    }
    result <- sensitivity_controlled_bind_observations(observed)
    sensitivity_controlled_validate_observations(result, design, registry)
    result
}

sensitivity_controlled_preflight <- function(
    full_design, registry, scorer, sensitivity, backend
) {
    design <- sensitivity_controlled_smoke_design(full_design)
    pieces <- lapply(seq_len(nrow(design)), function(index) {
        sensitivity_controlled_observe_scenario(
            design[index, , drop = FALSE], registry, scorer, sensitivity,
            backend, verify_encoding = TRUE
        )
    })
    observed <- sensitivity_controlled_bind_observations(pieces)
    sensitivity_controlled_validate_observations(observed, design, registry)
    invisible(observed)
}

sensitivity_controlled_gate_counts <- function(observed, design, registry) {
    expected_feature <- as.integer(sensitivity_registry_value(
        registry, "controlled", "feature_rows"
    ))
    expected_technical <- as.integer(sensitivity_registry_value(
        registry, "controlled", "technical_rows"
    ))
    valid <- nrow(design) == 5760L && nrow(observed$feature) == expected_feature &&
        nrow(observed$technical) == expected_technical &&
        length(unique(observed$technical$scenario_id)) == 5760L
    if (!valid) {
        stop("Sensitivity controlled gate is incomplete.", call. = FALSE)
    }
    invisible(observed)
}

sensitivity_controlled_metadata <- function(
    root, git, package, options, observed, analysis
) {
    metadata <- benchmark_metadata(
        root, "benchmark/run-sensitivity-controlled.R",
        paste0("sensitivity_controlled_", options$mode),
        1L, options$workers
    )
    metadata$value[metadata$key == "protocol_version"] <-
        SENSITIVITY_CONTROLLED_VERSION
    additions <- data.frame(
        key = c(
            "parent_protocol", "parent_registry_md5", "execution_protocol_md5",
            "candidate_id", "mode", "chunk_size", "source_archive_md5",
            "source_archive_sha256", "installed_manifest_md5", "feature_rows",
            "technical_rows", "latent_scenarios", "bootstrap_repeats_per_target",
            "scientific_decision"
        ),
        value = c(
            SENSITIVITY_PROTOCOL_VERSION, SENSITIVITY_PROTOCOL_MD5,
            SENSITIVITY_CONTROLLED_MD5, git$candidate, options$mode,
            options$chunk_size, package$archive_md5, package$archive_sha256,
            package$installed_manifest_md5, nrow(observed$feature),
            nrow(observed$technical), length(unique(observed$technical$scenario_id)),
            nrow(analysis$bootstrap) / 2L, options$mode == "gate"
        ),
        stringsAsFactors = FALSE, check.names = FALSE
    )
    rbind(metadata, additions)
}

sensitivity_controlled_write_analysis <- function(
    output_dir, observed, analysis, mode
) {
    benchmark_write_tsv(
        observed$feature, file.path(output_dir, "feature-observations.tsv")
    )
    benchmark_write_tsv(
        observed$technical, file.path(output_dir, "technical-observations.tsv")
    )
    folds <- rbind(analysis$feature_cv$folds, analysis$technical_cv$folds)
    coefficients <- rbind(
        analysis$feature_cv$coefficients, analysis$technical_cv$coefficients
    )
    scaling <- rbind(
        analysis$feature_cv$scaling, analysis$technical_cv$scaling
    )
    summary <- cbind(
        mode = mode, scientific_decision = mode == "gate", analysis$summary
    )
    benchmark_write_tsv(folds, file.path(output_dir, "fold-results.tsv"))
    benchmark_write_tsv(
        analysis$feature_cv$predictions,
        file.path(output_dir, "feature-predictions.tsv")
    )
    benchmark_write_tsv(
        analysis$technical_cv$predictions,
        file.path(output_dir, "technical-predictions.tsv")
    )
    benchmark_write_tsv(
        coefficients, file.path(output_dir, "model-coefficients.tsv")
    )
    benchmark_write_tsv(scaling, file.path(output_dir, "model-scaling.tsv"))
    benchmark_write_tsv(
        analysis$bootstrap, file.path(output_dir, "bootstrap.tsv")
    )
    benchmark_write_tsv(
        analysis$endpoints, file.path(output_dir, "endpoints.tsv")
    )
    benchmark_write_tsv(analysis$strata, file.path(output_dir, "strata.tsv"))
    benchmark_write_tsv(summary, file.path(output_dir, "summary.tsv"))
    summary
}

sensitivity_controlled_report <- function(
    path, analysis, summary, metadata, options, git
) {
    endpoints <- analysis$endpoints
    endpoints$estimate <- benchmark_format_number(endpoints$estimate)
    endpoints$threshold <- benchmark_format_number(endpoints$threshold)
    endpoints$passed <- ifelse(endpoints$passed, "yes", "NO")
    folds <- rbind(analysis$feature_cv$folds, analysis$technical_cv$folds)
    for (field in c("baseline_rmse", "augmented_rmse", "rmse_reduction")) {
        folds[[field]] <- benchmark_format_number(folds[[field]])
    }
    outcome <- if (options$mode == "smoke") {
        "SMOKE ONLY - planted model frames; no scientific decision."
    } else if (summary$controlled_gate_pass) {
        "PASS - both frozen controlled prediction gates passed."
    } else "FAIL - at least one frozen controlled prediction gate failed."
    command <- paste0(
        "candidate_id=$(git rev-parse HEAD)\nR_LIBS_USER=\"$PWD/.agent/R-library\" ",
        "Rscript --vanilla benchmark/run-sensitivity-controlled.R --mode=gate ",
        "--candidate-id=\"$candidate_id\" --workers=", options$workers,
        " --chunk-size=", options$chunk_size,
        " --output=benchmark/results/sensitivity-controlled"
    )
    lines <- c(
        "# GeneFunnel controlled sensitivity validation", "",
        paste0("Candidate `", git$candidate, "`; `", SENSITIVITY_PROTOCOL_VERSION,
            "` + `", SENSITIVITY_CONTROLLED_VERSION, "`. ", outcome), "",
        paste0(
            "This synthetic internal gate tests incremental prediction under the ",
            "frozen generator. It does not validate real technical or biological ",
            "replicates and never permits a public sensitivity API."
        ), "", "## Co-primary endpoints", "",
        benchmark_markdown_table(endpoints), "", "## Held-out folds", "",
        benchmark_markdown_table(folds), "", "## Environment", "",
        benchmark_markdown_table(benchmark_metadata_rows(metadata)), "",
        "## Reproduce", "", "```sh", command, "```", ""
    )
    writeLines(lines, path)
}

sensitivity_controlled_artifacts <- function(output_dir, files) {
    paths <- file.path(output_dir, files)
    if (any(!file.exists(paths))) {
        stop("Sensitivity controlled evidence is incomplete.", call. = FALSE)
    }
    data.frame(
        artifact = files, bytes = as.numeric(file.info(paths)$size),
        sha256 = unname(tools::sha256sum(paths)),
        stringsAsFactors = FALSE, check.names = FALSE
    )
}
