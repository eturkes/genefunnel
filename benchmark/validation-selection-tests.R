# Assisted-by: OpenAI Codex.

validation_test_write_tsv <- function(table, path) {
    utils::write.table(
        table, path, sep = "\t", quote = FALSE, row.names = FALSE,
        col.names = TRUE, na = "", eol = "\n", fileEncoding = "UTF-8"
    )
}

validation_test_identity <- function(rows) {
    data.frame(
        schema_version = rep(VALIDATION_SELECTION_SCHEMA, rows),
        protocol_version = rep(VALIDATION_AMENDMENT_VERSION, rows),
        protocol_sha256 = rep(VALIDATION_AMENDMENT_SHA256, rows),
        stringsAsFactors = FALSE, check.names = FALSE
    )
}

validation_test_fixture <- function(bundle_dir) {
    dir.create(bundle_dir, recursive = TRUE)
    assays <- c("bulk_rnaseq", "pseudobulk_rnaseq", "bulk_proteomics")
    prefixes <- c("bulk", "pseudo", "protein")
    sources <- do.call(rbind, lapply(seq_along(assays), function(a) {
        ids <- paste0(prefixes[[a]], c("_dev", "_h1", "_h2", "_h3"))
        data.frame(
            source_id = ids,
            assay = assays[[a]],
            access_role = c("development", rep("heldout", 3L)),
            source_record_id = toupper(ids),
            metadata_locator = paste0("https://metadata.example/", ids),
            license_id = "CC0-1.0",
            license_locator = "https://creativecommons.org/publicdomain/zero/1.0/",
            retrieval_permission = "scripted_retrieval_and_analysis",
            redistribution_permission = "allowed",
            independence_group = paste0(ids, "_collection"),
            stringsAsFactors = FALSE
        )
    }))
    sources$source_order <- as.character(seq_len(nrow(sources)))
    sources <- cbind(validation_test_identity(nrow(sources)), sources[c(
        "source_id", "source_order", "assay", "access_role",
        "source_record_id", "metadata_locator", "license_id",
        "license_locator", "retrieval_permission", "redistribution_permission",
        "independence_group"
    )])

    task_rows <- list()
    next_task <- 0L
    next_seed <- 260800000L
    append_task <- function(
        source_id, task_kind, contrast_id, population_id,
        reference_condition_id, intervention_condition_id, unit_frame_id,
        dependence_group_id, target_id, expected_direction, citation_id,
        citation_locator, null_assignment_mode, assignment_seed
    ) {
        next_task <<- next_task + 1L
        task_rows[[next_task]] <<- data.frame(
            task_id = paste0("task_", sprintf("%03d", next_task)),
            task_order = as.character(next_task), source_id = source_id,
            task_kind = task_kind, contrast_id = contrast_id,
            population_id = population_id,
            reference_condition_id = reference_condition_id,
            intervention_condition_id = intervention_condition_id,
            exposure_id = "registered_exposure",
            collection_time_id = "registered_time",
            unit_frame_id = unit_frame_id,
            dependence_group_id = dependence_group_id, target_id = target_id,
            expected_direction = expected_direction, citation_id = citation_id,
            citation_locator = citation_locator,
            null_assignment_mode = null_assignment_mode,
            label_seed = as.character(assignment_seed),
            permutation_seed = as.character(if (assignment_seed == 0L) {
                0L
            } else {
                assignment_seed + 1000000L
            }),
            stringsAsFactors = FALSE
        )
    }
    for (a in seq_along(assays)) {
        prefix <- prefixes[[a]]
        population <- if (assays[[a]] == "pseudobulk_rnaseq") {
            "monocyte"
        } else {
            "all"
        }
        append_task(
            paste0(prefix, "_dev"), "positive", "dev_c1", population,
            "control", "stimulus_a", "dev_frame", paste0(prefix, "_dev_dep"),
            "R-HSA-1059683", "+1", "PMID:100",
            "https://pubmed.ncbi.nlm.nih.gov/100/", "none", 0L
        )
        null_counts <- c(7L, 7L, 6L)
        for (study in seq_len(3L)) {
            source_id <- paste0(prefix, "_h", study)
            dependence <- paste0(source_id, "_dep")
            append_task(
                source_id, "positive", "contrast_a", population, "control",
                "stimulus_a", "frame_a", dependence, "R-HSA-1059683", "+1",
                "PMID:100", "https://pubmed.ncbi.nlm.nih.gov/100/", "none", 0L
            )
            append_task(
                source_id, "positive", "contrast_a", population, "control",
                "stimulus_a", "frame_a", dependence, "R-HSA-109581", "+1",
                "PMID:101", "https://pubmed.ncbi.nlm.nih.gov/101/", "none", 0L
            )
            append_task(
                source_id, "positive", "contrast_b", population, "control",
                "stimulus_b", "frame_b", dependence, "R-HSA-109606", "-1",
                "PMID:200", "https://pubmed.ncbi.nlm.nih.gov/200/", "none", 0L
            )
            for (draw in seq_len(null_counts[[study]])) {
                next_seed <- next_seed + 1L
                append_task(
                    source_id, "null", paste0("null_", study, "_", draw),
                    population, "control", "artificial_labels",
                    "null_frame", dependence, "all_eligible_pathways", "0",
                    "not_applicable", "not_applicable", "unpaired_balance",
                    next_seed
                )
            }
        }
    }
    tasks <- do.call(rbind, task_rows)
    tasks <- cbind(validation_test_identity(nrow(tasks)), tasks)

    design_rows <- do.call(rbind, lapply(seq_len(nrow(tasks)), function(i) {
        data.frame(
            task_id = tasks$task_id[[i]], field_order = "1",
            field_role = "unit", metadata_field = "biological_unit_id",
            field_type = "identifier",
            stringsAsFactors = FALSE
        )
    }))
    design <- cbind(validation_test_identity(nrow(design_rows)), design_rows)

    resampling_rows <- list()
    next_pool_row <- 0L
    for (i in seq_len(nrow(tasks))) {
        task <- tasks[i, , drop = FALSE]
        if (task$task_kind == "null") {
            rows <- list(c(
                resample_pool_id = paste0(task$source_id, "_null_pool"),
                pool_role = "joint", arm = "pool"
            ))
        } else {
            rows <- list(
                c(
                    resample_pool_id = paste0(
                        task$source_id, "_", task$unit_frame_id,
                        "_control_pool"
                    ),
                    pool_role = "separate", arm = "control"
                ),
                c(
                    resample_pool_id = paste0(
                        task$source_id, "_", task$unit_frame_id, "_",
                        task$intervention_condition_id, "_pool"
                    ),
                    pool_role = "separate", arm = "intervention"
                )
            )
        }
        for (j in seq_along(rows)) {
            next_pool_row <- next_pool_row + 1L
            resampling_rows[[next_pool_row]] <- data.frame(
                task_id = task$task_id, pool_order = as.character(j),
                resample_pool_id = rows[[j]][["resample_pool_id"]],
                pool_role = rows[[j]][["pool_role"]],
                arm = rows[[j]][["arm"]], stringsAsFactors = FALSE
            )
        }
    }
    resampling <- do.call(rbind, resampling_rows)
    resampling <- cbind(
        validation_test_identity(nrow(resampling)), resampling
    )

    filter_rows <- list()
    next_filter <- 0L
    for (i in seq_len(nrow(tasks))) {
        task <- tasks[i, , drop = FALSE]
        rows <- list()
        if (task$task_kind == "positive") {
            if (task$population_id != "all") {
                rows[[length(rows) + 1L]] <- c(
                    arm = "both", filter_role = "population",
                    metadata_scope = "cell",
                    metadata_field = "source_cell_type",
                    metadata_value = "CD14 monocyte"
                )
            }
            rows[[length(rows) + 1L]] <- c(
                arm = "both", filter_role = "exposure",
                metadata_scope = "sample",
                metadata_field = "source_exposure",
                metadata_value = "registered exposure"
            )
            rows[[length(rows) + 1L]] <- c(
                arm = "both", filter_role = "collection_time",
                metadata_scope = "sample",
                metadata_field = "source_time",
                metadata_value = "registered time"
            )
            rows[[length(rows) + 1L]] <- c(
                arm = "control", filter_role = "condition",
                metadata_scope = "sample",
                metadata_field = "source_arm", metadata_value = "control"
            )
            rows[[length(rows) + 1L]] <- c(
                arm = "intervention", filter_role = "condition",
                metadata_scope = "sample",
                metadata_field = "source_arm",
                metadata_value = task$intervention_condition_id
            )
        } else {
            if (task$population_id != "all") {
                rows[[length(rows) + 1L]] <- c(
                    arm = "pool", filter_role = "population",
                    metadata_scope = "cell",
                    metadata_field = "source_cell_type",
                    metadata_value = "CD14 monocyte"
                )
            }
            rows[[length(rows) + 1L]] <- c(
                arm = "pool", filter_role = "exposure",
                metadata_scope = "sample",
                metadata_field = "source_exposure",
                metadata_value = "registered exposure"
            )
            rows[[length(rows) + 1L]] <- c(
                arm = "pool", filter_role = "collection_time",
                metadata_scope = "sample",
                metadata_field = "source_time",
                metadata_value = "registered time"
            )
            rows[[length(rows) + 1L]] <- c(
                arm = "pool", filter_role = "condition",
                metadata_scope = "sample",
                metadata_field = "source_arm", metadata_value = "control"
            )
        }
        for (j in seq_along(rows)) {
            next_filter <- next_filter + 1L
            filter_rows[[next_filter]] <- data.frame(
                task_id = task$task_id, filter_order = as.character(j),
                arm = rows[[j]][["arm"]],
                filter_role = rows[[j]][["filter_role"]],
                metadata_scope = rows[[j]][["metadata_scope"]],
                metadata_field = rows[[j]][["metadata_field"]],
                metadata_value = rows[[j]][["metadata_value"]],
                stringsAsFactors = FALSE
            )
        }
    }
    filters <- do.call(rbind, filter_rows)
    filters <- cbind(validation_test_identity(nrow(filters)), filters)

    value_ids <- paste0("values_", sources$source_id)
    sample_ids <- paste0("samples_", sources$source_id)
    pseudo_source <- sources$assay == "pseudobulk_rnaseq"
    cell_ids <- paste0("cells_", sources$source_id[pseudo_source])
    mapping_id <- "mapping_protein_dev"
    value_kind <- rep("molecular_values", nrow(sources))
    value_kind[sources$source_id == "bulk_dev"] <- "combined_archive"
    value_suffix <- rep("values.tsv", nrow(sources))
    value_suffix[sources$source_id == "bulk_dev"] <- "values.zip"
    sample_kind <- rep("sample_metadata", nrow(sources))
    sample_kind[sources$source_id == "bulk_h2"] <- "metadata_archive"
    sample_suffix <- rep("samples.tsv", nrow(sources))
    sample_suffix[sources$source_id == "bulk_h2"] <- "samples.tar.gz"
    objects <- rbind(
        data.frame(
            object_id = value_ids, source_id = sources$source_id,
            object_kind = value_kind,
            access_class = "molecular_sealed",
            request_locator = paste0(
                "https://objects.example/", sources$source_id, "/", value_suffix
            ), stringsAsFactors = FALSE
        ),
        data.frame(
            object_id = sample_ids, source_id = sources$source_id,
            object_kind = sample_kind,
            access_class = "metadata_allowlisted",
            request_locator = paste0(
                "https://objects.example/", sources$source_id, "/", sample_suffix
            ), stringsAsFactors = FALSE
        ),
        data.frame(
            object_id = cell_ids,
            source_id = sources$source_id[pseudo_source],
            object_kind = "cell_metadata",
            access_class = "metadata_allowlisted",
            request_locator = paste0(
                "https://objects.example/", sources$source_id[pseudo_source],
                "/cells.tsv"
            ), stringsAsFactors = FALSE
        ),
        data.frame(
            object_id = mapping_id, source_id = "protein_dev",
            object_kind = "mapping_metadata",
            access_class = "mapping_allowlisted",
            request_locator =
                "https://objects.example/protein_dev/accession-map.tsv",
            stringsAsFactors = FALSE
        )
    )
    objects$object_order <- as.character(seq_len(nrow(objects)))
    objects <- objects[c(
        "object_id", "object_order", "source_id", "object_kind",
        "access_class", "request_locator"
    )]
    objects <- cbind(validation_test_identity(nrow(objects)), objects)

    rna <- sources$assay %in% c("bulk_rnaseq", "pseudobulk_rnaseq")
    protein <- sources$assay == "bulk_proteomics"
    cell_object <- rep("none", nrow(sources))
    cell_object[pseudo_source] <- cell_ids
    values_member <- rep("whole_object", nrow(sources))
    values_member[sources$source_id == "bulk_dev"] <- "matrix/counts.tsv"
    values_container <- rep("none", nrow(sources))
    values_container[sources$source_id == "bulk_dev"] <- "zip"
    values_container[sources$source_id == "bulk_h3"] <- "gzip"
    values_adapter <- rep("delimited_wide", nrow(sources))
    values_adapter[sources$source_id == "bulk_h1"] <- "delimited_long"
    values_adapter[sources$source_id == "pseudo_dev"] <- "tenx_h5"
    values_adapter[sources$source_id == "protein_dev"] <-
        "maxquant_protein_groups"
    values_adapter[sources$source_id == "protein_h1"] <-
        "fragpipe_combined_protein"
    long <- values_adapter == "delimited_long"
    matrix_orientation <- rep("features_by_observations", nrow(sources))
    matrix_orientation[long] <- "long_observation_feature"
    matrix_orientation[sources$source_id == "bulk_h3"] <-
        "observations_by_features"
    observation_id_field <- rep("column_names", nrow(sources))
    feature_id_field <- rep("row_names", nrow(sources))
    observation_id_field[long] <- "sample_id"
    feature_id_field[long] <- "ensembl_gene_id"
    observation_id_field[sources$source_id == "bulk_h3"] <- "row_names"
    feature_id_field[sources$source_id == "bulk_h3"] <- "column_names"
    observation_id_field[values_adapter == "tenx_h5"] <- "barcodes"
    feature_id_field[values_adapter == "tenx_h5"] <- "feature_id"
    observation_id_field[values_adapter == "maxquant_protein_groups"] <-
        "sample_columns"
    feature_id_field[values_adapter == "maxquant_protein_groups"] <-
        "majority_protein_ids"
    observation_id_field[values_adapter == "fragpipe_combined_protein"] <-
        "sample_columns"
    feature_id_field[values_adapter == "fragpipe_combined_protein"] <-
        "protein_id"
    quantity_field <- ifelse(long, "raw_count", "matrix_values")
    quantity_measure <- ifelse(rna, "raw_count", "abundance")
    quantity_measure[sources$source_id == "protein_dev"] <- "maxlfq"
    quantity_measure[sources$source_id == "protein_h1"] <- "intensity"
    quantity_selector <- ifelse(
        long, "exact_quantity_field", "exact_observation_columns"
    )
    quantity_selector[matrix_orientation == "observations_by_features"] <-
        "exact_observation_rows"
    quantity_selector[values_adapter == "tenx_h5"] <- "tenx_matrix"
    quantity_selector[values_adapter == "maxquant_protein_groups"] <-
        "maxquant_lfq_intensity_columns"
    quantity_selector[values_adapter == "fragpipe_combined_protein"] <-
        "fragpipe_intensity_columns"
    sample_member <- rep("whole_object", nrow(sources))
    sample_member[sources$source_id == "bulk_h2"] <- "metadata/samples.tsv"
    sample_container <- rep("none", nrow(sources))
    sample_container[sources$source_id == "bulk_h2"] <- "tar_gzip"
    technical_id <- rep("none", nrow(sources))
    technical_aggregation <- rep("none", nrow(sources))
    technical_id[sources$source_id == "bulk_h3"] <- "technical_run_id"
    technical_aggregation[sources$source_id == "bulk_h3"] <-
        "sum_raw_counts"
    technical_id[sources$source_id == "protein_h2"] <- "technical_run_id"
    technical_aggregation[sources$source_id == "protein_h2"] <-
        "arithmetic_mean_linear"
    mapping_present <- sources$source_id == "protein_dev"
    inputs <- data.frame(
        source_id = sources$source_id,
        values_object_id = value_ids, values_member = values_member,
        values_container = values_container, values_adapter = values_adapter,
        matrix_orientation = matrix_orientation,
        observation_id_field = observation_id_field,
        feature_id_field = feature_id_field, quantity_field = quantity_field,
        quantity_measure = quantity_measure,
        quantity_selector = quantity_selector,
        observation_selection_rule = ifelse(
            pseudo_source, "exact_cell_ids_from_cell_metadata",
            "exact_sample_ids_from_sample_metadata"
        ),
        sample_metadata_object_id = sample_ids,
        sample_metadata_member = sample_member,
        sample_metadata_container = sample_container,
        sample_metadata_adapter = "tsv_header", sample_id_field = "sample_id",
        technical_repeat_id_field = technical_id,
        technical_repeat_aggregation = technical_aggregation,
        cell_metadata_object_id = cell_object,
        cell_metadata_member = ifelse(pseudo_source, "whole_object", "none"),
        cell_metadata_container = "none",
        cell_metadata_adapter = ifelse(pseudo_source, "tsv_header", "none"),
        cell_id_field = ifelse(pseudo_source, "cell_id", "none"),
        cell_sample_id_field = ifelse(pseudo_source, "sample_id", "none"),
        quantity_level = ifelse(rna, "gene", "protein"),
        quantity_scale = ifelse(
            rna, "raw_nonnegative_integer_counts", "linear_nonnegative"
        ),
        normalization_state = ifelse(
            rna, "raw_counts", ifelse(
                sources$source_id == "protein_dev",
                "source_normalized_linear", "raw_linear"
            )
        ),
        imputation_state = "none",
        zero_semantics = "measured_zero",
        missing_semantics = ifelse(
            rna, "missing_forbidden", "explicit_tokens"
        ),
        missing_tokens = ifelse(rna, "none", "<blank>|NA"),
        assignment_or_inference = ifelse(
            rna, "direct_gene_identifier", "source_protein_inference"
        ),
        species_rule = ifelse(
            rna, "human_ensembl_only", "human_only_before_values"
        ),
        group_rule = ifelse(rna, "not_applicable", "exclude_ambiguous_groups"),
        duplicate_accession_operator = ifelse(
            rna, "not_applicable", ifelse(
                sources$source_id == "protein_dev",
                "source_nonoverlapping_linear_sum", "require_unique"
            )
        ),
        mapping_object_id = ifelse(mapping_present, mapping_id, "none"),
        mapping_member = ifelse(mapping_present, "whole_object", "none"),
        mapping_container = "none",
        mapping_adapter = ifelse(mapping_present, "tsv_header", "none"),
        mapping_source_field = ifelse(mapping_present, "accession", "none"),
        mapping_target_field = ifelse(
            mapping_present, "canonical_accession", "none"
        ),
        mapping_rule = ifelse(
            rna, "strip_ensembl_version", "canonical_uniprot_strip_isoform"
        ),
        evidence_locator = paste0(
            "https://metadata.example/", sources$source_id, "/input-contract"
        ), stringsAsFactors = FALSE
    )
    inputs <- cbind(validation_test_identity(nrow(inputs)), inputs)
    inputs <- inputs[VALIDATION_SELECTION_FIELDS$inputs]

    tables <- list(
        sources, tasks, design, resampling, filters, inputs, objects
    )
    for (i in seq_along(tables)) {
        validation_test_write_tsv(
            tables[[i]], file.path(bundle_dir, VALIDATION_SELECTION_FILES[[i]])
        )
    }
    invisible(bundle_dir)
}

validation_test_rejected <- function(expression) {
    inherits(try(force(expression), silent = TRUE), "try-error")
}

validation_test_production_access <- function(
    root, scratch, runtime_drift = FALSE, selection_tree_drift = FALSE,
    early_closure = FALSE, runtime_symlink_mode = FALSE
) {
    production_dir <- file.path(
        scratch, if (runtime_drift) {
            "production-access-runtime-drift"
        } else if (selection_tree_drift) {
            "production-access-selection-tree-drift"
        } else if (early_closure) {
            "production-access-early-closure"
        } else if (runtime_symlink_mode) {
            "production-access-runtime-symlink-mode"
        } else {
            "production-access"
        }
    )
    git_root <- file.path(production_dir, "repository")
    object_root <- file.path(production_dir, "objects")
    access_path <- file.path(production_dir, "access.tsv")
    dir.create(file.path(git_root, "benchmark"), recursive = TRUE)
    dir.create(object_root)
    root_files <- c(
        "validation-protocol.R", "validation-protocol.tsv",
        "validation-amendment.R", "validation-amendment.tsv",
        "validation-reactome97.tsv",
        "generate-validation-reactome97.R", "validation-selection.R"
    )
    copied <- file.copy(
        file.path(root, "benchmark", root_files),
        file.path(git_root, "benchmark", root_files)
    )
    if (!all(copied)) stop("Cannot construct production fixture.", call. = FALSE)
    if (runtime_drift) {
        drift_path <- file.path(git_root, "benchmark", "validation-amendment.R")
        connection <- file(drift_path, open = "ab")
        writeBin(charToRaw("# preselection semantic drift fixture\n"), connection)
        close(connection)
    }
    bundle <- file.path(git_root, "selection")
    validation_test_fixture(bundle)
    selection_extra <- file.path(bundle, "unregistered-at-selection.tsv")
    if (selection_tree_drift) {
        connection <- file(selection_extra, open = "wb")
        writeBin(charToRaw("unregistered\n"), connection)
        close(connection)
    }

    git <- function(arguments, stdout = FALSE) {
        suppressWarnings(system2(
            "git", c("-C", shQuote(git_root), arguments),
            stdout = stdout, stderr = FALSE, env = "LC_ALL=C"
        ))
    }
    marker_path <- file.path(git_root, "benchmark", "F-E-1.0.0.marker")
    closure_path <- file.path(
        git_root, "benchmark", "validation-execution-closure.R"
    )
    write_fixture_closure <- function() {
        marker_connection <- file(marker_path, open = "wb")
        writeBin(charToRaw("fixture execution closure\n"), marker_connection)
        close(marker_connection)
        marker_sha256 <- validation_contract_sha256(marker_path)
        closure_text <- paste0(
            "validation_execution_closure_validate <- function(root, ",
            "implementation_commit, validator_sha256) {\n",
            "    output <- tempfile(\"fixture-execution-marker-\")\n",
            "    on.exit(unlink(output), add = TRUE)\n",
            "    spec <- paste0(implementation_commit, ",
            "\":benchmark/F-E-1.0.0.marker\")\n",
            "    status <- suppressWarnings(system2(\"git\", c(\"-C\", ",
            "shQuote(root), \"show\", shQuote(spec)), stdout = output, ",
            "stderr = FALSE, env = \"LC_ALL=C\"))\n",
            "    observed <- unname(tools::sha256sum(output))\n",
            "    if (!identical(status, 0L) || !identical(observed, ",
            "\"", marker_sha256, "\")) stop(\"invalid fixture closure\")\n",
            "    list(schema_version = \"F-E-1.0.0\", ",
            "implementation_commit = implementation_commit, ",
            "validator_sha256 = validator_sha256, executable = TRUE)\n",
            "}\n"
        )
        closure_connection <- file(closure_path, open = "wb")
        writeBin(charToRaw(closure_text), closure_connection)
        close(closure_connection)
        validation_contract_sha256(closure_path)
    }
    if (early_closure) closure_sha256 <- write_fixture_closure()
    statuses <- c(
        git(c("init", "--object-format=sha1", "--quiet")),
        if (runtime_symlink_mode) {
            git(c("config", "core.symlinks", "false"))
        } else {
            0L
        },
        git(c("config", "user.name", shQuote("Fixture"))),
        git(c("config", "user.email", shQuote("fixture@example.invalid"))),
        git(c("config", "core.autocrlf", "false")),
        git(c("config", "commit.gpgsign", "false")),
        git(c("add", "--all")),
        git(c(
            "-c", "core.hooksPath=/dev/null", "commit", "--quiet",
            "--no-gpg-sign", "-m", shQuote("freeze selection")
        ))
    )
    selection_commit <- git(c("rev-parse", "HEAD"), stdout = TRUE)
    if (!all(statuses == 0L) || length(selection_commit) != 1L ||
        !grepl("^[0-9a-f]{40}$", selection_commit)) {
        stop("Cannot freeze production selection fixture.", call. = FALSE)
    }
    if (runtime_symlink_mode) {
        relative <- "benchmark/validation-amendment.R"
        blob <- git(c(
            "rev-parse", paste0("HEAD:", relative)
        ), stdout = TRUE)
        mode_statuses <- c(
            git(c(
                "update-index", "--cacheinfo", "120000", blob,
                shQuote(relative)
            )),
            git(c(
                "-c", "core.hooksPath=/dev/null", "commit", "--amend",
                "--quiet", "--no-gpg-sign", "--no-edit"
            ))
        )
        selection_commit <- git(c("rev-parse", "HEAD"), stdout = TRUE)
        restore_mode <- git(c(
            "update-index", "--cacheinfo", "100644", blob,
            shQuote(relative)
        ))
        if (length(blob) != 1L || !all(mode_statuses == 0L) ||
            length(selection_commit) != 1L || !identical(restore_mode, 0L)) {
            stop("Cannot construct runtime mode fixture.", call. = FALSE)
        }
    }
    if (runtime_drift && !file.copy(
            file.path(root, "benchmark", "validation-amendment.R"),
            file.path(git_root, "benchmark", "validation-amendment.R"),
            overwrite = TRUE
        )) {
        stop("Cannot restore production semantic runtime fixture.", call. = FALSE)
    }
    if (selection_tree_drift) unlink(selection_extra)
    selection <- validation_selection_validate(bundle, git_root)

    manifest_path <- file.path(git_root, "bytes.tsv")
    manifest <- data.frame(
        schema_version = VALIDATION_BYTE_SCHEMA,
        protocol_version = VALIDATION_AMENDMENT_VERSION,
        protocol_sha256 = VALIDATION_AMENDMENT_SHA256,
        selection_commit = selection_commit,
        selection_sha256 = selection$selection_sha256,
        object_id = selection$objects$object_id,
        status = "unavailable", resolved_locator = "none",
        terminal_http_status = "none",
        retrieved_utc = "2026-07-31T01:00:00Z", bytes = "0",
        sha256 = "none", failure_code = "transport_error",
        stringsAsFactors = FALSE, check.names = FALSE
    )
    validation_test_write_tsv(manifest, manifest_path)
    byte_add <- c("add", "--all")
    statuses <- c(
        git(byte_add),
        git(c(
            "-c", "core.hooksPath=/dev/null", "commit", "--quiet",
            "--no-gpg-sign", "-m", shQuote("freeze byte manifest")
        ))
    )
    byte_commit <- git(c("rev-parse", "HEAD"), stdout = TRUE)
    if (!all(statuses == 0L) || length(byte_commit) != 1L ||
        identical(byte_commit, selection_commit)) {
        stop("Cannot freeze production byte fixture.", call. = FALSE)
    }

    if (!early_closure) closure_sha256 <- write_fixture_closure()
    implementation_marker <- file.path(
        git_root, "benchmark", "implementation-freeze.marker"
    )
    marker_connection <- file(implementation_marker, open = "wb")
    writeBin(charToRaw("implementation freeze\n"), marker_connection)
    close(marker_connection)
    statuses <- c(
        git(c(
            "add", shQuote("benchmark/F-E-1.0.0.marker"),
            shQuote("benchmark/validation-execution-closure.R"),
            shQuote("benchmark/implementation-freeze.marker")
        )),
        git(c(
            "-c", "core.hooksPath=/dev/null", "commit", "--quiet",
            "--no-gpg-sign", "-m", shQuote("freeze implementation")
        ))
    )
    implementation_commit <- git(c("rev-parse", "HEAD"), stdout = TRUE)
    if (!all(statuses == 0L) || length(implementation_commit) != 1L ||
        identical(implementation_commit, byte_commit)) {
        stop("Cannot freeze production implementation fixture.", call. = FALSE)
    }

    expected <- validation_access_expected(selection, manifest)
    event_count <- length(expected$event_type)
    repository_head <- rep(selection_commit, event_count)
    repository_head[
        expected$byte_freeze:(expected$implementation_freeze - 1L)
    ] <- byte_commit
    repository_head[expected$implementation_freeze:event_count] <-
        implementation_commit
    byte_sha256 <- validation_contract_sha256(manifest_path)
    byte_manifest_sha256 <- rep("none", event_count)
    byte_manifest_sha256[expected$byte_freeze:event_count] <- byte_sha256
    implementation <- rep("none", event_count)
    implementation[expected$implementation_freeze:event_count] <-
        implementation_commit
    execution_closure <- rep("none", event_count)
    execution_closure[expected$implementation_freeze:event_count] <-
        closure_sha256
    event_utc <- rep("2026-07-31T01:00:00Z", event_count)
    event_utc[expected$stream_rows] <- manifest$retrieved_utc
    attempt_id <- paste0("attempt_", paste(rep("c", 32L), collapse = ""))
    access <- data.frame(
        schema_version = VALIDATION_ACCESS_SCHEMA,
        protocol_version = VALIDATION_AMENDMENT_VERSION,
        protocol_sha256 = VALIDATION_AMENDMENT_SHA256,
        event_order = as.character(seq_len(event_count)),
        event_type = expected$event_type,
        repository_head = repository_head, worktree_state = "clean",
        selection_sha256 = selection$selection_sha256,
        byte_manifest_sha256 = byte_manifest_sha256,
        implementation_commit = implementation,
        execution_closure_sha256 = execution_closure,
        execution_attempt_id = attempt_id,
        object_id = expected$object_id, event_utc = event_utc,
        bytes = expected$bytes, sha256 = expected$sha256,
        stringsAsFactors = FALSE, check.names = FALSE
    )
    validation_test_write_tsv(access, access_path)

    validator_environment <- environment(validation_access_git_artifacts)
    validator_name <- "validation_execution_closure_validate"
    had_validator <- exists(
        validator_name, envir = validator_environment, inherits = FALSE
    )
    if (had_validator) {
        old_validator <- get(
            validator_name, envir = validator_environment, inherits = FALSE
        )
    }
    on.exit({
        if (had_validator) {
            assign(validator_name, old_validator, envir = validator_environment)
        } else if (exists(
            validator_name, envir = validator_environment, inherits = FALSE
        )) {
            rm(list = validator_name, envir = validator_environment)
        }
    }, add = TRUE)
    assign(validator_name, function(...) {
        stop("ambient execution closure must not run", call. = FALSE)
    }, envir = validator_environment)

    production_result <- try(validation_access_log_validate(
        bundle, manifest_path, access_path, git_root, object_root,
        verify_git = TRUE
    ), silent = TRUE)
    if (runtime_drift || selection_tree_drift || early_closure ||
        runtime_symlink_mode) {
        expected_error <- if (selection_tree_drift) {
            "Selection commit bundle tree closure is invalid"
        } else if (early_closure) {
            "Final F-E-1.0.0 closure predates implementation freeze"
        } else {
            "Semantic runtime changed after selection freeze"
        }
        return(list(
            fault_rejected = inherits(production_result, "try-error") &&
                grepl(expected_error, as.character(production_result),
                    fixed = TRUE)
        ))
    }
    if (inherits(production_result, "try-error")) {
        stop("Production fixture unexpectedly rejected: ", production_result,
            call. = FALSE)
    }
    production <- production_result
    bad_closure_access <- access
    bad_closure_access$execution_closure_sha256[
        expected$implementation_freeze:event_count
    ] <- paste(rep("0", 64L), collapse = "")
    validation_test_write_tsv(bad_closure_access, access_path)
    closure_hash_rejected <- validation_test_rejected(
        validation_access_log_validate(
            bundle, manifest_path, access_path, git_root, object_root,
            verify_git = TRUE
        )
    )
    validation_test_write_tsv(access, access_path)
    replace_created <- identical(git(c(
        "replace", selection_commit, byte_commit
    )), 0L)
    replace_ref_rejected <- replace_created && validation_test_rejected(
        validation_access_log_validate(
            bundle, manifest_path, access_path, git_root, object_root,
            verify_git = TRUE
        )
    )
    replace_deleted <- replace_created && identical(git(c(
        "replace", "-d", selection_commit
    )), 0L)
    if (!replace_deleted) {
        stop("Cannot remove production replacement-ref fixture.", call. = FALSE)
    }
    marker_connection <- file(marker_path, open = "ab")
    writeBin(as.raw(1L), marker_connection)
    close(marker_connection)
    dirty_rejected <- validation_test_rejected(validation_access_log_validate(
        bundle, manifest_path, access_path, git_root, object_root,
        verify_git = TRUE
    ))
    marker_connection <- file(marker_path, open = "wb")
    writeBin(charToRaw("fixture execution closure\n"), marker_connection)
    close(marker_connection)
    descendant_path <- file.path(git_root, "benchmark", "descendant.marker")
    descendant_connection <- file(descendant_path, open = "wb")
    writeBin(charToRaw("clean descendant\n"), descendant_connection)
    close(descendant_connection)
    statuses <- c(
        git(c("add", shQuote("benchmark/descendant.marker"))),
        git(c(
            "-c", "core.hooksPath=/dev/null", "commit", "--quiet",
            "--no-gpg-sign", "-m", shQuote("descendant checkout")
        ))
    )
    descendant_rejected <- all(statuses == 0L) &&
        validation_test_rejected(validation_access_log_validate(
            bundle, manifest_path, access_path, git_root, object_root,
            verify_git = TRUE
        ))
    list(
        production_git_valid = isTRUE(production$executable) &&
            production$events == event_count &&
            production$objects == nrow(manifest) &&
            identical(production$attempt_id, attempt_id),
        dirty_worktree_rejected = dirty_rejected,
        descendant_head_rejected = descendant_rejected,
        ambient_closure_ignored = isTRUE(production$executable),
        closure_hash_rejected = closure_hash_rejected,
        replace_ref_rejected = replace_ref_rejected,
        events = event_count, objects = nrow(manifest)
    )
}

validation_selection_contract_self_test <- function(root, scratch) {
    bundle <- file.path(scratch, "selection")
    validation_test_fixture(bundle)
    valid <- validation_selection_validate(bundle, root)
    if (!identical(
            valid$selection_sha256,
            validation_selection_bundle_sha256(bundle)
        )) {
        stop("Selection fixture hash is unstable.", call. = FALSE)
    }

    fifo_path <- file.path(scratch, "contract-fifo.tsv")
    fifo_status <- suppressWarnings(system2("mkfifo", shQuote(fifo_path)))
    fifo_contract_rejected <- identical(fifo_status, 0L) &&
        validation_test_rejected(validation_contract_read_tsv(
            fifo_path, VALIDATION_SELECTION_FIELDS$sources
        ))
    unlink(fifo_path)

    git_root <- file.path(scratch, "git fixture")
    git_bundle <- file.path(git_root, "bundle space")
    validation_test_fixture(git_bundle)
    git_status <- c(
        system2("git", c(
            "-C", shQuote(git_root), "init", "--object-format=sha1", "--quiet"
        )),
        system2("git", c(
            "-C", shQuote(git_root), "config", "user.name", "Fixture"
        )),
        system2("git", c(
            "-C", shQuote(git_root), "config", "user.email",
            "fixture@example.invalid"
        )),
        system2("git", c("-C", shQuote(git_root), "add", "--all")),
        system2("git", c(
            "-C", shQuote(git_root), "commit", "--quiet", "--no-gpg-sign",
            "-m", shQuote("selection fixture")
        ))
    )
    git_commit <- system2(
        "git", c("-C", shQuote(git_root), "rev-parse", "HEAD"),
        stdout = TRUE
    )
    git_selection <- validation_selection_validate(git_bundle, root)
    git_spaced_path_valid <- all(git_status == 0L) &&
        length(git_commit) == 1L && !validation_test_rejected(
            validation_byte_git_selection(
                git_selection, git_bundle, git_root, git_commit
            )
        )
    prior_git_dir <- Sys.getenv("GIT_DIR", unset = NA_character_)
    restore_git_dir <- function() {
        if (is.na(prior_git_dir)) {
            Sys.unsetenv("GIT_DIR")
        } else {
            Sys.setenv(GIT_DIR = prior_git_dir)
        }
    }
    on.exit(restore_git_dir(), add = TRUE)
    Sys.setenv(GIT_DIR = file.path(scratch, "hostile-inherited-git-dir"))
    git_env_ignored <- !validation_test_rejected(
        validation_byte_git_selection(
            git_selection, git_bundle, git_root, git_commit
        )
    )
    restore_git_dir()
    hostile_git_name <- "GIT_X /usr/bin/false; printf FORGED_GIT_OUTPUT; #"
    prior_hostile_git <- Sys.getenv(hostile_git_name, unset = NA_character_)
    restore_hostile_git <- function() {
        if (is.na(prior_hostile_git)) {
            Sys.unsetenv(hostile_git_name)
        } else {
            do.call(Sys.setenv, as.list(stats::setNames(
                prior_hostile_git, hostile_git_name
            )))
        }
    }
    on.exit(restore_hostile_git(), add = TRUE)
    do.call(Sys.setenv, as.list(stats::setNames("1", hostile_git_name)))
    git_env_name_ignored <- !validation_test_rejected(
        validation_byte_git_selection(
            git_selection, git_bundle, git_root, git_commit
        )
    )
    restore_hostile_git()
    fixture_git <- function(arguments) suppressWarnings(system2(
        "git", c("-C", shQuote(git_root), arguments),
        stdout = FALSE, stderr = FALSE, env = "LC_ALL=C"
    ))
    include_path <- file.path(git_root, ".git", "config.include")
    include_status <- c(
        fixture_git(c(
            "config", "--file", shQuote(include_path),
            "filter.hostile.clean", "/bin/false"
        )),
        fixture_git(c("config", "include.path", "config.include"))
    )
    include_config_rejected <- all(include_status == 0L) &&
        validation_test_rejected(validation_byte_git_selection(
            git_selection, git_bundle, git_root, git_commit
        ))
    include_cleanup <- fixture_git(c(
        "config", "--unset-all", "include.path"
    ))
    unlink(include_path)
    if (!identical(include_cleanup, 0L)) {
        stop("Cannot clean included Git config fixture.", call. = FALSE)
    }
    worktree_status <- c(
        fixture_git(c("config", "extensions.worktreeConfig", "true")),
        fixture_git(c(
            "config", "--worktree", "filter.hostile.clean", "/bin/false"
        ))
    )
    worktree_config_rejected <- all(worktree_status == 0L) &&
        validation_test_rejected(validation_byte_git_selection(
            git_selection, git_bundle, git_root, git_commit
        ))
    worktree_cleanup <- c(
        fixture_git(c(
            "config", "--worktree", "--unset-all", "filter.hostile.clean"
        )),
        fixture_git(c(
            "config", "--unset-all", "extensions.worktreeConfig"
        ))
    )
    if (!all(worktree_cleanup == 0L)) {
        stop("Cannot clean worktree Git config fixture.", call. = FALSE)
    }

    mutation <- file.path(scratch, "mutation")
    copy_bundle <- function() {
        unlink(mutation, recursive = TRUE, force = TRUE)
        dir.create(mutation)
        copied <- file.copy(
            file.path(bundle, VALIDATION_SELECTION_FILES), mutation
        )
        if (!all(copied)) stop("Cannot copy selection fixture.", call. = FALSE)
    }

    copy_bundle()
    sources <- utils::read.delim(
        file.path(mutation, "sources.tsv"), colClasses = "character",
        check.names = FALSE
    )
    sources$access_role[sources$source_id == "bulk_h1"] <- "development"
    validation_test_write_tsv(sources, file.path(mutation, "sources.tsv"))
    role_rejected <- validation_test_rejected(
        validation_selection_validate(mutation, root)
    )

    copy_bundle()
    extra_path <- file.path(mutation, "unregistered.tsv")
    extra_connection <- file(extra_path, open = "wb")
    writeBin(charToRaw("unregistered\n"), extra_connection)
    close(extra_connection)
    bundle_extra_rejected <- validation_test_rejected(
        validation_selection_validate(mutation, root)
    )

    copy_bundle()
    sources <- utils::read.delim(
        file.path(mutation, "sources.tsv"), colClasses = "character",
        check.names = FALSE
    )
    sources$license_locator[[1L]] <-
        "https://licenses.example/contradictory-terms"
    validation_test_write_tsv(sources, file.path(mutation, "sources.tsv"))
    license_binding_rejected <- validation_test_rejected(
        validation_selection_validate(mutation, root)
    )

    copy_bundle()
    sources <- utils::read.delim(
        file.path(mutation, "sources.tsv"), colClasses = "character",
        check.names = FALSE
    )
    bulk_row <- which(sources$source_id == "bulk_h1")
    protein_row <- which(sources$source_id == "protein_h1")
    sources$source_record_id[[protein_row]] <-
        sources$source_record_id[[bulk_row]]
    sources$metadata_locator[[protein_row]] <-
        sources$metadata_locator[[bulk_row]]
    validation_test_write_tsv(sources, file.path(mutation, "sources.tsv"))
    source_identity_rejected <- validation_test_rejected(
        validation_selection_validate(mutation, root)
    )

    copy_bundle()
    sources <- utils::read.delim(
        file.path(mutation, "sources.tsv"), colClasses = "character",
        check.names = FALSE
    )
    sources$source_record_id[[1L]] <- "none"
    validation_test_write_tsv(sources, file.path(mutation, "sources.tsv"))
    source_sentinel_rejected <- validation_test_rejected(
        validation_selection_validate(mutation, root)
    )

    copy_bundle()
    sources <- utils::read.delim(
        file.path(mutation, "sources.tsv"), colClasses = "character",
        check.names = FALSE
    )
    sources$independence_group[[1L]] <- "none"
    validation_test_write_tsv(sources, file.path(mutation, "sources.tsv"))
    independence_sentinel_rejected <- validation_test_rejected(
        validation_selection_validate(mutation, root)
    )

    copy_bundle()
    sources <- utils::read.delim(
        file.path(mutation, "sources.tsv"), colClasses = "character",
        check.names = FALSE
    )
    sources$source_id[[1L]] <- "none"
    validation_test_write_tsv(sources, file.path(mutation, "sources.tsv"))
    source_id_sentinel_rejected <- validation_test_rejected(
        validation_selection_validate(mutation, root)
    )

    copy_bundle()
    tasks <- utils::read.delim(
        file.path(mutation, "tasks.tsv"), colClasses = "character",
        check.names = FALSE
    )
    filters <- utils::read.delim(
        file.path(mutation, "task-filters.tsv"), colClasses = "character",
        check.names = FALSE
    )
    pseudo_task <- tasks$task_id[tasks$source_id == "pseudo_h1"][[1L]]
    filters <- filters[!(filters$task_id == pseudo_task &
        filters$filter_role == "population"), , drop = FALSE]
    validation_test_write_tsv(
        filters, file.path(mutation, "task-filters.tsv")
    )
    population_rejected <- validation_test_rejected(
        validation_selection_validate(mutation, root)
    )

    copy_bundle()
    filters <- utils::read.delim(
        file.path(mutation, "task-filters.tsv"), colClasses = "character",
        check.names = FALSE
    )
    population_row <- which(filters$filter_role == "population")[[1L]]
    filters$metadata_scope[[population_row]] <- "sample"
    validation_test_write_tsv(
        filters, file.path(mutation, "task-filters.tsv")
    )
    metadata_scope_rejected <- validation_test_rejected(
        validation_selection_validate(mutation, root)
    )

    copy_bundle()
    tasks <- utils::read.delim(
        file.path(mutation, "tasks.tsv"), colClasses = "character",
        check.names = FALSE
    )
    filters <- utils::read.delim(
        file.path(mutation, "task-filters.tsv"), colClasses = "character",
        check.names = FALSE
    )
    shared <- tasks$task_id[
        tasks$source_id == "bulk_h1" & tasks$contrast_id == "contrast_a"
    ]
    row <- which(filters$task_id == shared[[2L]] &
        filters$arm == "intervention")[[1L]]
    filters$metadata_value[[row]] <- "changed_after_target_split"
    validation_test_write_tsv(
        filters, file.path(mutation, "task-filters.tsv")
    )
    contrast_rejected <- validation_test_rejected(
        validation_selection_validate(mutation, root)
    )

    copy_bundle()
    tasks <- utils::read.delim(
        file.path(mutation, "tasks.tsv"), colClasses = "character",
        check.names = FALSE
    )
    tasks$dependence_group_id[tasks$source_id == "bulk_h2"] <- "bulk_h1_dep"
    validation_test_write_tsv(tasks, file.path(mutation, "tasks.tsv"))
    dependence_study_rejected <- validation_test_rejected(
        validation_selection_validate(mutation, root)
    )

    copy_bundle()
    tasks <- utils::read.delim(
        file.path(mutation, "tasks.tsv"), colClasses = "character",
        check.names = FALSE
    )
    split_rows <- tasks$source_id == "bulk_h1" &
        tasks$unit_frame_id == "frame_b"
    tasks$dependence_group_id[split_rows] <- "bulk_h1_split_dep"
    validation_test_write_tsv(tasks, file.path(mutation, "tasks.tsv"))
    dependence_split_rejected <- validation_test_rejected(
        validation_selection_validate(mutation, root)
    )

    copy_bundle()
    tasks <- utils::read.delim(
        file.path(mutation, "tasks.tsv"), colClasses = "character",
        check.names = FALSE
    )
    null_rows <- which(tasks$task_kind == "null")
    tasks$label_seed[null_rows[[2L]]] <- tasks$label_seed[null_rows[[1L]]]
    validation_test_write_tsv(tasks, file.path(mutation, "tasks.tsv"))
    seed_rejected <- validation_test_rejected(
        validation_selection_validate(mutation, root)
    )

    copy_bundle()
    sources <- utils::read.delim(
        file.path(mutation, "sources.tsv"), colClasses = "character",
        check.names = FALSE
    )
    sources$bytes <- "0"
    validation_test_write_tsv(sources, file.path(mutation, "sources.tsv"))
    forbidden_rejected <- validation_test_rejected(
        validation_selection_validate(mutation, root)
    )

    copy_bundle()
    source_path <- file.path(mutation, "sources.tsv")
    source_lines <- readLines(source_path, warn = FALSE)
    source_lines[-1L] <- paste0("injected\t", source_lines[-1L])
    writeLines(source_lines, source_path, useBytes = TRUE)
    implicit_rownames_rejected <- validation_test_rejected(
        validation_selection_validate(mutation, root)
    )

    copy_bundle()
    tasks <- utils::read.delim(
        file.path(mutation, "tasks.tsv"), colClasses = "character",
        check.names = FALSE
    )
    tasks$target_id[which(tasks$task_kind == "positive")[[1L]]] <-
        "R-HSA-999999999999"
    validation_test_write_tsv(tasks, file.path(mutation, "tasks.tsv"))
    roster_rejected <- validation_test_rejected(
        validation_selection_validate(mutation, root)
    )

    copy_bundle()
    tasks <- utils::read.delim(
        file.path(mutation, "tasks.tsv"), colClasses = "character",
        check.names = FALSE
    )
    filters <- utils::read.delim(
        file.path(mutation, "task-filters.tsv"), colClasses = "character",
        check.names = FALSE
    )
    task_id <- tasks$task_id[
        tasks$source_id == "bulk_h1" & tasks$contrast_id == "contrast_a"
    ][[1L]]
    row <- which(filters$task_id == task_id & filters$arm == "control")[[1L]]
    extra <- filters[row, , drop = FALSE]
    extra$arm <- "both"
    extra$filter_order <- as.character(sum(filters$task_id == task_id) + 1L)
    filters <- rbind(filters, extra)
    validation_test_write_tsv(
        filters, file.path(mutation, "task-filters.tsv")
    )
    extra_condition_rejected <- validation_test_rejected(
        validation_selection_validate(mutation, root)
    )

    copy_bundle()
    design <- utils::read.delim(
        file.path(mutation, "task-design.tsv"), colClasses = "character",
        check.names = FALSE
    )
    design$metadata_field[[1L]] <- "source_exposure"
    validation_test_write_tsv(design, file.path(mutation, "task-design.tsv"))
    filtered_model_field_rejected <- validation_test_rejected(
        validation_selection_validate(mutation, root)
    )

    copy_bundle()
    design <- utils::read.delim(
        file.path(mutation, "task-design.tsv"), colClasses = "character",
        check.names = FALSE
    )
    design$metadata_field[[1L]] <- "none"
    validation_test_write_tsv(design, file.path(mutation, "task-design.tsv"))
    design_field_rejected <- validation_test_rejected(
        validation_selection_validate(mutation, root)
    )

    copy_bundle()
    design <- utils::read.delim(
        file.path(mutation, "task-design.tsv"), colClasses = "character",
        check.names = FALSE
    )
    design$field_type[[1L]] <- "categorical"
    validation_test_write_tsv(design, file.path(mutation, "task-design.tsv"))
    design_type_rejected <- validation_test_rejected(
        validation_selection_validate(mutation, root)
    )

    copy_bundle()
    design <- utils::read.delim(
        file.path(mutation, "task-design.tsv"), colClasses = "character",
        check.names = FALSE
    )
    extra <- design[1L, , drop = FALSE]
    extra$field_order <- "2"
    extra$field_role <- "covariate"
    extra$field_type <- "numeric"
    design <- rbind(design, extra)
    validation_test_write_tsv(design, file.path(mutation, "task-design.tsv"))
    design_role_alias_rejected <- validation_test_rejected(
        validation_selection_validate(mutation, root)
    )

    copy_bundle()
    design <- utils::read.delim(
        file.path(mutation, "task-design.tsv"), colClasses = "character",
        check.names = FALSE
    )
    extras <- design[rep(1L, 2L), , drop = FALSE]
    extras$field_order <- c("2", "3")
    extras$field_role <- "block"
    extras$metadata_field <- c("registered_block_a", "registered_block_b")
    extras$field_type <- "categorical"
    design <- rbind(design, extras)
    validation_test_write_tsv(design, file.path(mutation, "task-design.tsv"))
    multiple_block_rejected <- validation_test_rejected(
        validation_selection_validate(mutation, root)
    )

    copy_bundle()
    tasks <- utils::read.delim(
        file.path(mutation, "tasks.tsv"), colClasses = "character",
        check.names = FALSE
    )
    design <- utils::read.delim(
        file.path(mutation, "task-design.tsv"), colClasses = "character",
        check.names = FALSE
    )
    resampling <- utils::read.delim(
        file.path(mutation, "task-resampling.tsv"), colClasses = "character",
        check.names = FALSE
    )
    null_id <- tasks$task_id[
        tasks$source_id == "bulk_h1" & tasks$task_kind == "null"
    ][[2L]]
    block <- design[design$task_id == null_id, , drop = FALSE]
    block$field_order <- "2"
    block$field_role <- "block"
    block$metadata_field <- "null_assignment_block"
    block$field_type <- "categorical"
    design <- rbind(design, block)
    resampling$resample_pool_id[resampling$task_id == null_id] <-
        "bulk_h1_fragmented_null_pool"
    validation_test_write_tsv(design, file.path(mutation, "task-design.tsv"))
    validation_test_write_tsv(
        resampling, file.path(mutation, "task-resampling.tsv")
    )
    nuisance_pool_fragment_rejected <- validation_test_rejected(
        validation_selection_validate(mutation, root)
    )

    copy_bundle()
    tasks <- utils::read.delim(
        file.path(mutation, "tasks.tsv"), colClasses = "character",
        check.names = FALSE
    )
    design <- utils::read.delim(
        file.path(mutation, "task-design.tsv"), colClasses = "character",
        check.names = FALSE
    )
    resampling <- utils::read.delim(
        file.path(mutation, "task-resampling.tsv"), colClasses = "character",
        check.names = FALSE
    )
    frame_task <- tasks$task_id[
        tasks$source_id == "bulk_h1" & tasks$contrast_id == "contrast_b"
    ][[1L]]
    tasks$unit_frame_id[tasks$task_id == frame_task] <- "frame_a"
    design$metadata_field[
        design$task_id == frame_task & design$field_role == "unit"
    ] <- "alternate_biological_unit_id"
    resampling$resample_pool_id[resampling$task_id == frame_task] <- paste0(
        resampling$resample_pool_id[resampling$task_id == frame_task], "_alternate"
    )
    validation_test_write_tsv(tasks, file.path(mutation, "tasks.tsv"))
    validation_test_write_tsv(design, file.path(mutation, "task-design.tsv"))
    validation_test_write_tsv(
        resampling, file.path(mutation, "task-resampling.tsv")
    )
    frame_unit_binding_rejected <- validation_test_rejected(
        validation_selection_validate(mutation, root)
    )

    copy_bundle()
    tasks <- utils::read.delim(
        file.path(mutation, "tasks.tsv"), colClasses = "character",
        check.names = FALSE
    )
    resampling <- utils::read.delim(
        file.path(mutation, "task-resampling.tsv"), colClasses = "character",
        check.names = FALSE
    )
    null_id <- tasks$task_id[
        tasks$source_id == "bulk_h1" & tasks$task_kind == "null"
    ][[1L]]
    tasks$unit_frame_id[tasks$task_id == null_id] <- "frame_a"
    resampling$resample_pool_id[resampling$task_id == null_id] <-
        "bulk_h1_fresh_control_equivalent_pool"
    validation_test_write_tsv(tasks, file.path(mutation, "tasks.tsv"))
    validation_test_write_tsv(
        resampling, file.path(mutation, "task-resampling.tsv")
    )
    membership_pool_fragment_rejected <- validation_test_rejected(
        validation_selection_validate(mutation, root)
    )

    copy_bundle()
    filters <- utils::read.delim(
        file.path(mutation, "task-filters.tsv"), colClasses = "character",
        check.names = FALSE
    )
    filters$metadata_field[[1L]] <- "*"
    validation_test_write_tsv(
        filters, file.path(mutation, "task-filters.tsv")
    )
    filter_field_rejected <- validation_test_rejected(
        validation_selection_validate(mutation, root)
    )

    copy_bundle()
    filters <- utils::read.delim(
        file.path(mutation, "task-filters.tsv"), colClasses = "character",
        check.names = FALSE
    )
    task_id <- filters$task_id[[1L]]
    extra <- filters[1L, , drop = FALSE]
    extra$filter_order <- as.character(sum(filters$task_id == task_id) + 1L)
    extra$filter_role <- "subset"
    filters <- rbind(filters, extra)
    validation_test_write_tsv(
        filters, file.path(mutation, "task-filters.tsv")
    )
    predicate_role_alias_rejected <- validation_test_rejected(
        validation_selection_validate(mutation, root)
    )

    copy_bundle()
    filters <- utils::read.delim(
        file.path(mutation, "task-filters.tsv"), colClasses = "character",
        check.names = FALSE
    )
    tasks <- utils::read.delim(
        file.path(mutation, "tasks.tsv"), colClasses = "character",
        check.names = FALSE
    )
    bulk_tasks <- tasks$task_id[tasks$source_id == "bulk_h1"]
    contradictory <- filters$task_id %in% bulk_tasks &
        filters$filter_role == "collection_time"
    filters$metadata_field[contradictory] <- "source_exposure"
    validation_test_write_tsv(
        filters, file.path(mutation, "task-filters.tsv")
    )
    contradictory_filter_rejected <- validation_test_rejected(
        validation_selection_validate(mutation, root)
    )

    copy_bundle()
    objects <- utils::read.delim(
        file.path(mutation, "objects-planned.tsv"), colClasses = "character",
        check.names = FALSE
    )
    objects$object_kind[[2L]] <- "combined_archive"
    validation_test_write_tsv(
        objects, file.path(mutation, "objects-planned.tsv")
    )
    archive_member_rejected <- validation_test_rejected(
        validation_selection_validate(mutation, root)
    )

    copy_bundle()
    objects <- utils::read.delim(
        file.path(mutation, "objects-planned.tsv"), colClasses = "character",
        check.names = FALSE
    )
    inputs <- utils::read.delim(
        file.path(mutation, "assay-inputs.tsv"), colClasses = "character",
        check.names = FALSE
    )
    objects$object_kind[[1L]] <- "combined_archive"
    inputs$values_member[[1L]] <- "."
    validation_test_write_tsv(
        objects, file.path(mutation, "objects-planned.tsv")
    )
    validation_test_write_tsv(inputs, file.path(mutation, "assay-inputs.tsv"))
    unsafe_member_rejected <- validation_test_rejected(
        validation_selection_validate(mutation, root)
    )

    copy_bundle()
    objects <- utils::read.delim(
        file.path(mutation, "objects-planned.tsv"), colClasses = "character",
        check.names = FALSE
    )
    unused <- objects[1L, , drop = FALSE]
    unused$object_id <- "unused_metadata"
    unused$object_order <- as.character(nrow(objects) + 1L)
    unused$object_kind <- "sample_metadata"
    unused$access_class <- "metadata_allowlisted"
    unused$request_locator <- "https://objects.example/unused-metadata.tsv"
    objects <- rbind(objects, unused)
    validation_test_write_tsv(
        objects, file.path(mutation, "objects-planned.tsv")
    )
    unused_object_rejected <- validation_test_rejected(
        validation_selection_validate(mutation, root)
    )

    copy_bundle()
    objects <- utils::read.delim(
        file.path(mutation, "objects-planned.tsv"), colClasses = "character",
        check.names = FALSE
    )
    inputs <- utils::read.delim(
        file.path(mutation, "assay-inputs.tsv"), colClasses = "character",
        check.names = FALSE
    )
    old_id <- objects$object_id[[1L]]
    objects$object_id[[1L]] <- "none"
    inputs$values_object_id[inputs$values_object_id == old_id] <- "none"
    validation_test_write_tsv(
        objects, file.path(mutation, "objects-planned.tsv")
    )
    validation_test_write_tsv(inputs, file.path(mutation, "assay-inputs.tsv"))
    reserved_object_rejected <- validation_test_rejected(
        validation_selection_validate(mutation, root)
    )

    copy_bundle()
    objects <- utils::read.delim(
        file.path(mutation, "objects-planned.tsv"), colClasses = "character",
        check.names = FALSE
    )
    metadata_row <- which(objects$access_class == "metadata_allowlisted")[[1L]]
    objects$request_locator[[metadata_row]] <- paste0(
        objects$request_locator[[1L]], "#metadata-copy"
    )
    validation_test_write_tsv(
        objects, file.path(mutation, "objects-planned.tsv")
    )
    fragment_alias_rejected <- validation_test_rejected(
        validation_selection_validate(mutation, root)
    )

    copy_bundle()
    objects <- utils::read.delim(
        file.path(mutation, "objects-planned.tsv"), colClasses = "character",
        check.names = FALSE
    )
    objects$request_locator[[2L]] <- sub(
        "objects.example", "OBJECTS.EXAMPLE",
        objects$request_locator[[1L]], fixed = TRUE
    )
    validation_test_write_tsv(
        objects, file.path(mutation, "objects-planned.tsv")
    )
    canonical_locator_rejected <- validation_test_rejected(
        validation_selection_validate(mutation, root)
    )

    copy_bundle()
    sources <- utils::read.delim(
        file.path(mutation, "sources.tsv"), colClasses = "character",
        check.names = FALSE
    )
    sources$metadata_locator[[1L]] <- "https://metadata.example/a|b"
    validation_test_write_tsv(sources, file.path(mutation, "sources.tsv"))
    illegal_uri_character_rejected <- validation_test_rejected(
        validation_selection_validate(mutation, root)
    )

    copy_bundle()
    sources <- utils::read.delim(
        file.path(mutation, "sources.tsv"), colClasses = "character",
        check.names = FALSE
    )
    label <- paste(rep("a", 63L), collapse = "")
    long_host <- paste(rep(label, 5L), collapse = ".")
    sources$metadata_locator[[1L]] <- paste0(
        "https://", long_host, "/record"
    )
    validation_test_write_tsv(sources, file.path(mutation, "sources.tsv"))
    long_authority_rejected <- validation_test_rejected(
        validation_selection_validate(mutation, root)
    )

    copy_bundle()
    sources <- utils::read.delim(
        file.path(mutation, "sources.tsv"), colClasses = "character",
        check.names = FALSE
    )
    objects <- utils::read.delim(
        file.path(mutation, "objects-planned.tsv"), colClasses = "character",
        check.names = FALSE
    )
    sources$metadata_locator[[1L]] <- objects$request_locator[[1L]]
    validation_test_write_tsv(sources, file.path(mutation, "sources.tsv"))
    evidence_alias_rejected <- validation_test_rejected(
        validation_selection_validate(mutation, root)
    )

    copy_bundle()
    inputs <- utils::read.delim(
        file.path(mutation, "assay-inputs.tsv"), colClasses = "character",
        check.names = FALSE
    )
    objects <- utils::read.delim(
        file.path(mutation, "objects-planned.tsv"), colClasses = "character",
        check.names = FALSE
    )
    inputs$evidence_locator[[1L]] <- paste0(
        objects$request_locator[[1L]], "#evidence"
    )
    validation_test_write_tsv(inputs, file.path(mutation, "assay-inputs.tsv"))
    canonical_evidence_alias_rejected <- validation_test_rejected(
        validation_selection_validate(mutation, root)
    )

    copy_bundle()
    inputs <- utils::read.delim(
        file.path(mutation, "assay-inputs.tsv"), colClasses = "character",
        check.names = FALSE
    )
    protein <- which(inputs$quantity_level == "protein")[[1L]]
    inputs$missing_tokens[[protein]] <- "NA|NA"
    validation_test_write_tsv(inputs, file.path(mutation, "assay-inputs.tsv"))
    input_semantics_rejected <- validation_test_rejected(
        validation_selection_validate(mutation, root)
    )

    copy_bundle()
    inputs <- utils::read.delim(
        file.path(mutation, "assay-inputs.tsv"), colClasses = "character",
        check.names = FALSE
    )
    protein <- which(inputs$quantity_level == "protein")[[1L]]
    inputs$missing_tokens[[protein]] <- "0.0|NA"
    validation_test_write_tsv(inputs, file.path(mutation, "assay-inputs.tsv"))
    zero_missing_rejected <- validation_test_rejected(
        validation_selection_validate(mutation, root)
    )

    copy_bundle()
    inputs <- utils::read.delim(
        file.path(mutation, "assay-inputs.tsv"), colClasses = "character",
        check.names = FALSE
    )
    protein <- which(inputs$quantity_level == "protein")[[1L]]
    inputs$missing_tokens[[protein]] <- "1e-9999|NA"
    validation_test_write_tsv(inputs, file.path(mutation, "assay-inputs.tsv"))
    underflow_missing_rejected <- validation_test_rejected(
        validation_selection_validate(mutation, root)
    )

    copy_bundle()
    inputs <- utils::read.delim(
        file.path(mutation, "assay-inputs.tsv"), colClasses = "character",
        check.names = FALSE
    )
    protein <- which(inputs$quantity_level == "protein")[[1L]]
    inputs$missing_tokens[[protein]] <- "NA|"
    validation_test_write_tsv(inputs, file.path(mutation, "assay-inputs.tsv"))
    terminal_missing_token_rejected <- validation_test_rejected(
        validation_selection_validate(mutation, root)
    )

    copy_bundle()
    inputs <- utils::read.delim(
        file.path(mutation, "assay-inputs.tsv"), colClasses = "character",
        check.names = FALSE
    )
    inputs$values_adapter[[1L]] <- "maxquant_protein_groups"
    validation_test_write_tsv(inputs, file.path(mutation, "assay-inputs.tsv"))
    adapter_rejected <- validation_test_rejected(
        validation_selection_validate(mutation, root)
    )

    copy_bundle()
    inputs <- utils::read.delim(
        file.path(mutation, "assay-inputs.tsv"), colClasses = "character",
        check.names = FALSE
    )
    inputs$sample_metadata_adapter[[1L]] <- "unregistered_adapter"
    validation_test_write_tsv(inputs, file.path(mutation, "assay-inputs.tsv"))
    metadata_adapter_rejected <- validation_test_rejected(
        validation_selection_validate(mutation, root)
    )

    copy_bundle()
    inputs <- utils::read.delim(
        file.path(mutation, "assay-inputs.tsv"), colClasses = "character",
        check.names = FALSE
    )
    inputs$sample_id_field[[1L]] <- "none"
    validation_test_write_tsv(inputs, file.path(mutation, "assay-inputs.tsv"))
    selector_rejected <- validation_test_rejected(
        validation_selection_validate(mutation, root)
    )

    copy_bundle()
    inputs <- utils::read.delim(
        file.path(mutation, "assay-inputs.tsv"), colClasses = "character",
        check.names = FALSE
    )
    inputs$observation_id_field[[1L]] <- inputs$feature_id_field[[1L]]
    validation_test_write_tsv(inputs, file.path(mutation, "assay-inputs.tsv"))
    axis_alias_rejected <- validation_test_rejected(
        validation_selection_validate(mutation, root)
    )

    copy_bundle()
    inputs <- utils::read.delim(
        file.path(mutation, "assay-inputs.tsv"), colClasses = "character",
        check.names = FALSE
    )
    inputs$values_adapter[[1L]] <- "delimited_long"
    inputs$matrix_orientation[[1L]] <- "long_observation_feature"
    inputs$observation_id_field[[1L]] <- "sample_column"
    inputs$feature_id_field[[1L]] <- "feature_column"
    inputs$quantity_field[[1L]] <- "sample_column"
    inputs$quantity_selector[[1L]] <- "exact_quantity_field"
    validation_test_write_tsv(inputs, file.path(mutation, "assay-inputs.tsv"))
    quantity_alias_rejected <- validation_test_rejected(
        validation_selection_validate(mutation, root)
    )

    copy_bundle()
    inputs <- utils::read.delim(
        file.path(mutation, "assay-inputs.tsv"), colClasses = "character",
        check.names = FALSE
    )
    pseudo <- which(inputs$cell_metadata_object_id != "none")[[1L]]
    inputs$cell_id_field[[pseudo]] <- inputs$cell_sample_id_field[[pseudo]]
    validation_test_write_tsv(inputs, file.path(mutation, "assay-inputs.tsv"))
    cell_alias_rejected <- validation_test_rejected(
        validation_selection_validate(mutation, root)
    )

    copy_bundle()
    inputs <- utils::read.delim(
        file.path(mutation, "assay-inputs.tsv"), colClasses = "character",
        check.names = FALSE
    )
    archive_row <- which(inputs$values_container == "zip")[[1L]]
    inputs$values_container[[archive_row]] <- "none"
    validation_test_write_tsv(inputs, file.path(mutation, "assay-inputs.tsv"))
    container_rejected <- validation_test_rejected(
        validation_selection_validate(mutation, root)
    )

    copy_bundle()
    inputs <- utils::read.delim(
        file.path(mutation, "assay-inputs.tsv"), colClasses = "character",
        check.names = FALSE
    )
    maxquant_row <- which(
        inputs$values_adapter == "maxquant_protein_groups"
    )[[1L]]
    inputs$quantity_selector[[maxquant_row]] <-
        "fragpipe_maxlfq_intensity_columns"
    validation_test_write_tsv(inputs, file.path(mutation, "assay-inputs.tsv"))
    quantity_selector_rejected <- validation_test_rejected(
        validation_selection_validate(mutation, root)
    )

    copy_bundle()
    inputs <- utils::read.delim(
        file.path(mutation, "assay-inputs.tsv"), colClasses = "character",
        check.names = FALSE
    )
    row_oriented <- which(
        inputs$matrix_orientation == "observations_by_features"
    )[[1L]]
    inputs$quantity_selector[[row_oriented]] <- "exact_observation_columns"
    validation_test_write_tsv(inputs, file.path(mutation, "assay-inputs.tsv"))
    orientation_selector_rejected <- validation_test_rejected(
        validation_selection_validate(mutation, root)
    )

    copy_bundle()
    inputs <- utils::read.delim(
        file.path(mutation, "assay-inputs.tsv"), colClasses = "character",
        check.names = FALSE
    )
    maxlfq_row <- which(inputs$quantity_measure == "maxlfq")[[1L]]
    inputs$normalization_state[[maxlfq_row]] <- "raw_linear"
    validation_test_write_tsv(inputs, file.path(mutation, "assay-inputs.tsv"))
    maxlfq_state_rejected <- validation_test_rejected(
        validation_selection_validate(mutation, root)
    )

    copy_bundle()
    inputs <- utils::read.delim(
        file.path(mutation, "assay-inputs.tsv"), colClasses = "character",
        check.names = FALSE
    )
    rna_row <- which(inputs$quantity_level == "gene")[[1L]]
    inputs$duplicate_accession_operator[[rna_row]] <- "require_unique"
    validation_test_write_tsv(inputs, file.path(mutation, "assay-inputs.tsv"))
    duplicate_operator_rejected <- validation_test_rejected(
        validation_selection_validate(mutation, root)
    )

    copy_bundle()
    inputs <- utils::read.delim(
        file.path(mutation, "assay-inputs.tsv"), colClasses = "character",
        check.names = FALSE
    )
    technical_row <- which(inputs$technical_repeat_id_field != "none")[[1L]]
    inputs$technical_repeat_id_field[[technical_row]] <-
        inputs$sample_id_field[[technical_row]]
    validation_test_write_tsv(inputs, file.path(mutation, "assay-inputs.tsv"))
    technical_sample_alias_rejected <- validation_test_rejected(
        validation_selection_validate(mutation, root)
    )

    copy_bundle()
    inputs <- utils::read.delim(
        file.path(mutation, "assay-inputs.tsv"), colClasses = "character",
        check.names = FALSE
    )
    technical_row <- which(inputs$technical_repeat_id_field != "none")[[1L]]
    inputs$technical_repeat_id_field[[technical_row]] <- "biological_unit_id"
    validation_test_write_tsv(inputs, file.path(mutation, "assay-inputs.tsv"))
    technical_model_alias_rejected <- validation_test_rejected(
        validation_selection_validate(mutation, root)
    )

    copy_bundle()
    inputs <- utils::read.delim(
        file.path(mutation, "assay-inputs.tsv"), colClasses = "character",
        check.names = FALSE
    )
    technical_row <- which(inputs$technical_repeat_id_field != "none")[[1L]]
    inputs$technical_repeat_id_field[[technical_row]] <- "source_arm"
    validation_test_write_tsv(inputs, file.path(mutation, "assay-inputs.tsv"))
    technical_filter_alias_rejected <- validation_test_rejected(
        validation_selection_validate(mutation, root)
    )

    copy_bundle()
    inputs <- utils::read.delim(
        file.path(mutation, "assay-inputs.tsv"), colClasses = "character",
        check.names = FALSE
    )
    mapping_row <- which(inputs$mapping_object_id != "none")[[1L]]
    inputs$mapping_target_field[[mapping_row]] <-
        inputs$mapping_source_field[[mapping_row]]
    validation_test_write_tsv(inputs, file.path(mutation, "assay-inputs.tsv"))
    mapping_selector_rejected <- validation_test_rejected(
        validation_selection_validate(mutation, root)
    )

    copy_bundle()
    inputs <- utils::read.delim(
        file.path(mutation, "assay-inputs.tsv"), colClasses = "character",
        check.names = FALSE
    )
    mapping_row <- which(inputs$mapping_object_id != "none")[[1L]]
    inputs$mapping_container[[mapping_row]] <- "zip"
    validation_test_write_tsv(inputs, file.path(mutation, "assay-inputs.tsv"))
    mapping_container_rejected <- validation_test_rejected(
        validation_selection_validate(mutation, root)
    )

    copy_bundle()
    inputs <- utils::read.delim(
        file.path(mutation, "assay-inputs.tsv"), colClasses = "character",
        check.names = FALSE
    )
    objects <- utils::read.delim(
        file.path(mutation, "objects-planned.tsv"), colClasses = "character",
        check.names = FALSE
    )
    pseudo_row <- which(inputs$cell_metadata_object_id != "none")[[1L]]
    sample_object <- inputs$sample_metadata_object_id[[pseudo_row]]
    inputs$sample_metadata_object_id[[pseudo_row]] <-
        inputs$cell_metadata_object_id[[pseudo_row]]
    inputs$sample_metadata_member[[pseudo_row]] <-
        inputs$cell_metadata_member[[pseudo_row]]
    inputs$sample_metadata_container[[pseudo_row]] <-
        inputs$cell_metadata_container[[pseudo_row]]
    objects <- objects[objects$object_id != sample_object, , drop = FALSE]
    objects$object_order <- as.character(seq_len(nrow(objects)))
    validation_test_write_tsv(inputs, file.path(mutation, "assay-inputs.tsv"))
    validation_test_write_tsv(
        objects, file.path(mutation, "objects-planned.tsv")
    )
    metadata_reference_alias_rejected <- validation_test_rejected(
        validation_selection_validate(mutation, root)
    )

    copy_bundle()
    inputs <- utils::read.delim(
        file.path(mutation, "assay-inputs.tsv"), colClasses = "character",
        check.names = FALSE
    )
    inputs$sample_id_field[[1L]] <- "source_arm"
    validation_test_write_tsv(inputs, file.path(mutation, "assay-inputs.tsv"))
    sample_join_filter_alias_rejected <- validation_test_rejected(
        validation_selection_validate(mutation, root)
    )

    copy_bundle()
    tasks <- utils::read.delim(
        file.path(mutation, "tasks.tsv"), colClasses = "character",
        check.names = FALSE
    )
    design <- utils::read.delim(
        file.path(mutation, "task-design.tsv"), colClasses = "character",
        check.names = FALSE
    )
    inputs <- utils::read.delim(
        file.path(mutation, "assay-inputs.tsv"), colClasses = "character",
        check.names = FALSE
    )
    source_id <- inputs$source_id[[1L]]
    task_id <- tasks$task_id[tasks$source_id == source_id][[1L]]
    block <- design[design$task_id == task_id, , drop = FALSE]
    block$field_order <- "2"
    block$field_role <- "block"
    block$metadata_field <- inputs$sample_id_field[[1L]]
    block$field_type <- "categorical"
    design <- rbind(design, block)
    validation_test_write_tsv(design, file.path(mutation, "task-design.tsv"))
    sample_join_model_alias_rejected <- validation_test_rejected(
        validation_selection_validate(mutation, root)
    )

    copy_bundle()
    inputs <- utils::read.delim(
        file.path(mutation, "assay-inputs.tsv"), colClasses = "character",
        check.names = FALSE
    )
    pseudo_row <- which(inputs$cell_metadata_object_id != "none")[[1L]]
    inputs$cell_id_field[[pseudo_row]] <- "source_cell_type"
    validation_test_write_tsv(inputs, file.path(mutation, "assay-inputs.tsv"))
    cell_join_filter_alias_rejected <- validation_test_rejected(
        validation_selection_validate(mutation, root)
    )

    copy_bundle()
    tasks <- utils::read.delim(
        file.path(mutation, "tasks.tsv"), colClasses = "character",
        check.names = FALSE
    )
    citation_rows <- which(tasks$citation_id == "PMID:100")
    tasks$citation_locator[citation_rows[[2L]]] <-
        "https://pubmed.ncbi.nlm.nih.gov/999/"
    validation_test_write_tsv(tasks, file.path(mutation, "tasks.tsv"))
    citation_binding_rejected <- validation_test_rejected(
        validation_selection_validate(mutation, root)
    )

    copy_bundle()
    tasks <- utils::read.delim(
        file.path(mutation, "tasks.tsv"), colClasses = "character",
        check.names = FALSE
    )
    cited <- which(tasks$task_kind == "positive")[[1L]]
    tasks$citation_id[[cited]] <- "none"
    validation_test_write_tsv(tasks, file.path(mutation, "tasks.tsv"))
    citation_sentinel_rejected <- validation_test_rejected(
        validation_selection_validate(mutation, root)
    )

    copy_bundle()
    tasks <- utils::read.delim(
        file.path(mutation, "tasks.tsv"), colClasses = "character",
        check.names = FALSE
    )
    filters <- utils::read.delim(
        file.path(mutation, "task-filters.tsv"), colClasses = "character",
        check.names = FALSE
    )
    design <- utils::read.delim(
        file.path(mutation, "task-design.tsv"), colClasses = "character",
        check.names = FALSE
    )
    resampling <- utils::read.delim(
        file.path(mutation, "task-resampling.tsv"), colClasses = "character",
        check.names = FALSE
    )
    source_id <- tasks$task_id[
        tasks$source_id == "bulk_h1" & tasks$contrast_id == "contrast_a"
    ][[1L]]
    duplicate_id <- tasks$task_id[
        tasks$source_id == "bulk_h1" & tasks$contrast_id == "contrast_b"
    ][[1L]]
    source_row <- match(source_id, tasks$task_id)
    duplicate_row <- match(duplicate_id, tasks$task_id)
    for (field in c(
        "population_id", "reference_condition_id", "intervention_condition_id",
        "exposure_id", "collection_time_id", "dependence_group_id"
    )) tasks[[field]][[duplicate_row]] <- tasks[[field]][[source_row]]
    tasks$unit_frame_id[[duplicate_row]] <- "frame_duplicate"
    block <- design[design$task_id == duplicate_id, , drop = FALSE]
    block$field_order <- "2"
    block$field_role <- "block"
    block$metadata_field <- "registered_duplicate_block"
    block$field_type <- "categorical"
    design <- rbind(design, block)
    copied_filters <- filters[filters$task_id == source_id, , drop = FALSE]
    copied_filters <- copied_filters[rev(seq_len(nrow(copied_filters))), , drop = FALSE]
    copied_filters$task_id <- duplicate_id
    copied_filters$filter_order <- as.character(seq_len(nrow(copied_filters)))
    filters <- rbind(
        filters[filters$task_id != duplicate_id, , drop = FALSE], copied_filters
    )
    copied_pools <- resampling[
        resampling$task_id == source_id, , drop = FALSE
    ]
    copied_pools <- copied_pools[rev(seq_len(nrow(copied_pools))), , drop = FALSE]
    copied_pools$task_id <- duplicate_id
    copied_pools$pool_order <- as.character(seq_len(nrow(copied_pools)))
    copied_pools$resample_pool_id <- paste0(
        copied_pools$resample_pool_id, "_duplicate"
    )
    resampling <- rbind(
        resampling[resampling$task_id != duplicate_id, , drop = FALSE],
        copied_pools
    )
    validation_test_write_tsv(tasks, file.path(mutation, "tasks.tsv"))
    validation_test_write_tsv(design, file.path(mutation, "task-design.tsv"))
    validation_test_write_tsv(
        filters, file.path(mutation, "task-filters.tsv")
    )
    validation_test_write_tsv(
        resampling, file.path(mutation, "task-resampling.tsv")
    )
    reordered_duplicate_rejected <- validation_test_rejected(
        validation_selection_validate(mutation, root)
    )

    copy_bundle()
    tasks <- utils::read.delim(
        file.path(mutation, "tasks.tsv"), colClasses = "character",
        check.names = FALSE
    )
    filters <- utils::read.delim(
        file.path(mutation, "task-filters.tsv"), colClasses = "character",
        check.names = FALSE
    )
    resampling <- utils::read.delim(
        file.path(mutation, "task-resampling.tsv"), colClasses = "character",
        check.names = FALSE
    )
    forward_id <- tasks$task_id[
        tasks$source_id == "bulk_h1" & tasks$contrast_id == "contrast_a"
    ][[1L]]
    reverse_id <- tasks$task_id[
        tasks$source_id == "bulk_h1" & tasks$contrast_id == "contrast_b"
    ][[1L]]
    forward_row <- match(forward_id, tasks$task_id)
    reverse_row <- match(reverse_id, tasks$task_id)
    for (field in c(
        "population_id", "exposure_id", "collection_time_id",
        "unit_frame_id", "dependence_group_id"
    )) tasks[[field]][[reverse_row]] <- tasks[[field]][[forward_row]]
    tasks$reference_condition_id[[reverse_row]] <-
        tasks$intervention_condition_id[[forward_row]]
    tasks$intervention_condition_id[[reverse_row]] <-
        tasks$reference_condition_id[[forward_row]]
    reverse_filters <- filters[
        filters$task_id == forward_id, , drop = FALSE
    ]
    reverse_filters$task_id <- reverse_id
    condition <- reverse_filters$filter_role == "condition"
    reverse_filters$metadata_value[
        condition & reverse_filters$arm == "control"
    ] <- tasks$reference_condition_id[[reverse_row]]
    reverse_filters$metadata_value[
        condition & reverse_filters$arm == "intervention"
    ] <- tasks$intervention_condition_id[[reverse_row]]
    filters <- rbind(
        filters[filters$task_id != reverse_id, , drop = FALSE],
        reverse_filters
    )
    forward_pools <- resampling[
        resampling$task_id == forward_id, , drop = FALSE
    ]
    reverse_pools <- forward_pools
    reverse_pools$task_id <- reverse_id
    reverse_pools$resample_pool_id <- rev(forward_pools$resample_pool_id)
    resampling <- rbind(
        resampling[resampling$task_id != reverse_id, , drop = FALSE],
        reverse_pools
    )
    validation_test_write_tsv(tasks, file.path(mutation, "tasks.tsv"))
    validation_test_write_tsv(
        filters, file.path(mutation, "task-filters.tsv")
    )
    validation_test_write_tsv(
        resampling, file.path(mutation, "task-resampling.tsv")
    )
    reversed_contrast_rejected <- validation_test_rejected(
        validation_selection_validate(mutation, root)
    )

    copy_bundle()
    tasks <- utils::read.delim(
        file.path(mutation, "tasks.tsv"), colClasses = "character",
        check.names = FALSE
    )
    design <- utils::read.delim(
        file.path(mutation, "task-design.tsv"), colClasses = "character",
        check.names = FALSE
    )
    filters <- utils::read.delim(
        file.path(mutation, "task-filters.tsv"), colClasses = "character",
        check.names = FALSE
    )
    resampling <- utils::read.delim(
        file.path(mutation, "task-resampling.tsv"), colClasses = "character",
        check.names = FALSE
    )
    positive_ids <- tasks$task_id[
        tasks$source_id == "bulk_h1" & tasks$contrast_id == "contrast_a"
    ]
    positive_id <- positive_ids[[1L]]
    paired_id <- tasks$task_id[
        tasks$source_id == "bulk_h1" & tasks$task_kind == "null"
    ][[1L]]
    positive_row <- match(positive_id, tasks$task_id)
    paired_row <- match(paired_id, tasks$task_id)
    for (field in c(
        "population_id", "reference_condition_id", "intervention_condition_id",
        "exposure_id", "collection_time_id", "unit_frame_id",
        "dependence_group_id", "citation_id", "citation_locator"
    )) tasks[[field]][[paired_row]] <- tasks[[field]][[positive_row]]
    tasks$null_assignment_mode[[paired_row]] <- "paired_swap"
    tasks$label_seed[[paired_row]] <- "0"
    paired_tasks <- c(positive_ids, paired_id)
    design <- design[!design$task_id %in% paired_tasks, , drop = FALSE]
    paired_design <- do.call(rbind, lapply(paired_tasks, function(id) {
        data.frame(
            validation_test_identity(2L), task_id = id,
            field_order = c("1", "2"), field_role = c("unit", "pair"),
            metadata_field = c("biological_unit_id", "pair_id"),
            field_type = c("identifier", "categorical"),
            stringsAsFactors = FALSE, check.names = FALSE
        )
    }))
    paired_block <- paired_design[
        paired_design$task_id == paired_id &
            paired_design$field_role == "unit", , drop = FALSE
    ]
    paired_block$field_order <- "3"
    paired_block$field_role <- "block"
    paired_block$metadata_field <- "registered_null_block"
    paired_block$field_type <- "categorical"
    paired_design <- rbind(paired_design, paired_block)
    design <- rbind(design, paired_design)
    copied_filters <- filters[
        filters$task_id == positive_id, , drop = FALSE
    ]
    copied_filters$task_id <- paired_id
    filters <- rbind(
        filters[filters$task_id != paired_id, , drop = FALSE], copied_filters
    )
    resampling <- resampling[
        !resampling$task_id %in% paired_tasks, , drop = FALSE
    ]
    paired_pools <- data.frame(
        validation_test_identity(length(paired_tasks)),
        task_id = paired_tasks, pool_order = "1",
        resample_pool_id = ifelse(
            paired_tasks == paired_id,
            "bulk_h1_null_paired_pool", "bulk_h1_exact_paired_pool"
        ),
        pool_role = "joint", arm = "joint", stringsAsFactors = FALSE,
        check.names = FALSE
    )
    resampling <- rbind(resampling, paired_pools)
    validation_test_write_tsv(tasks, file.path(mutation, "tasks.tsv"))
    validation_test_write_tsv(design, file.path(mutation, "task-design.tsv"))
    validation_test_write_tsv(
        filters, file.path(mutation, "task-filters.tsv")
    )
    validation_test_write_tsv(
        resampling, file.path(mutation, "task-resampling.tsv")
    )
    positive_null_overlap_rejected <- validation_test_rejected(
        validation_selection_validate(mutation, root)
    )

    copy_bundle()
    tasks <- utils::read.delim(
        file.path(mutation, "tasks.tsv"), colClasses = "character",
        check.names = FALSE
    )
    null_row <- which(tasks$task_kind == "null")[[1L]]
    tasks$permutation_seed[[null_row]] <- "2147483648"
    validation_test_write_tsv(tasks, file.path(mutation, "tasks.tsv"))
    seed_bound_rejected <- validation_test_rejected(
        validation_selection_validate(mutation, root)
    )

    manifest_path <- file.path(scratch, "bytes.tsv")
    manifest <- data.frame(
        schema_version = VALIDATION_BYTE_SCHEMA,
        protocol_version = VALIDATION_AMENDMENT_VERSION,
        protocol_sha256 = VALIDATION_AMENDMENT_SHA256,
        selection_commit = paste(rep("0", 40L), collapse = ""),
        selection_sha256 = valid$selection_sha256,
        object_id = valid$objects$object_id,
        status = "unavailable", resolved_locator = "none",
        terminal_http_status = "none",
        retrieved_utc = "2026-07-31T00:00:00Z", bytes = "0",
        sha256 = "none", failure_code = "transport_error",
        stringsAsFactors = FALSE
    )
    validation_test_write_tsv(manifest, manifest_path)
    byte_valid <- validation_byte_manifest_validate(
        bundle, manifest_path, root, verify_git = FALSE
    )

    bad_manifest <- manifest
    bad_manifest$resolved_locator[[1L]] <- "https://objects.example/missing.bin"
    bad_manifest$terminal_http_status[[1L]] <- "200"
    bad_manifest$failure_code[[1L]] <- "not_found"
    validation_test_write_tsv(bad_manifest, manifest_path)
    failure_mapping_rejected <- validation_test_rejected(
        validation_byte_manifest_validate(
            bundle, manifest_path, root, verify_git = FALSE
        )
    )
    validation_test_write_tsv(manifest, manifest_path)

    short_manifest <- manifest[-nrow(manifest), , drop = FALSE]
    validation_test_write_tsv(short_manifest, manifest_path)
    closure_rejected <- validation_test_rejected(
        validation_byte_manifest_validate(bundle, manifest_path, root)
    )

    manifest$bytes[[2L]] <- "0e0"
    validation_test_write_tsv(manifest, manifest_path)
    numeric_rejected <- validation_test_rejected(
        validation_byte_manifest_validate(bundle, manifest_path, root)
    )

    object_root <- file.path(scratch, "objects")
    verified_rows <- c(
        1L, 2L,
        which(valid$objects$access_class == "metadata_allowlisted")[[1L]]
    )
    source_role <- setNames(
        valid$sources$access_role, valid$sources$source_id
    )
    object_source <- setNames(
        valid$objects$source_id, valid$objects$object_id
    )
    local_paths <- vapply(verified_rows, function(i) {
        object_id <- manifest$object_id[[i]]
        file.path(
            object_root, source_role[[object_source[[object_id]]]],
            paste0(object_id, ".bin")
        )
    }, character(1L))
    for (i in seq_along(local_paths)) {
        dir.create(
            dirname(local_paths[[i]]), recursive = TRUE, showWarnings = FALSE
        )
        connection <- file(local_paths[[i]], open = "wb")
        writeBin(charToRaw(paste0("opaque fixture bytes ", i, "\n")), connection)
        close(connection)
    }
    manifest$bytes <- "0"
    manifest$status[verified_rows] <- "verified"
    manifest$resolved_locator[verified_rows] <- paste0(
        "https://objects.example/resolved-", verified_rows, ".bin"
    )
    manifest$terminal_http_status[verified_rows] <- "200"
    manifest$bytes[verified_rows] <- as.character(file.info(local_paths)$size)
    manifest$sha256[verified_rows] <- unname(tools::sha256sum(local_paths))
    manifest$failure_code[verified_rows] <- "none"
    validation_test_write_tsv(manifest, manifest_path)
    local_valid <- validation_byte_manifest_validate(
        bundle, manifest_path, root, object_root = object_root,
        verify_git = FALSE
    )

    cross_plan_manifest <- manifest
    heldout_row <- verified_rows[[2L]]
    metadata_row <- verified_rows[[3L]]
    cross_plan_manifest$status[[heldout_row]] <- "unavailable"
    cross_plan_manifest$resolved_locator[[heldout_row]] <- "none"
    cross_plan_manifest$terminal_http_status[[heldout_row]] <- "none"
    cross_plan_manifest$bytes[[heldout_row]] <- "0"
    cross_plan_manifest$sha256[[heldout_row]] <- "none"
    cross_plan_manifest$failure_code[[heldout_row]] <- "transport_error"
    cross_plan_manifest$resolved_locator[[metadata_row]] <-
        valid$objects$request_locator[[heldout_row]]
    heldout_backup <- file.path(scratch, "heldout-bytes.backup")
    moved <- file.rename(local_paths[[2L]], heldout_backup)
    validation_test_write_tsv(cross_plan_manifest, manifest_path)
    terminal_planned_alias_rejected <- moved && validation_test_rejected(
        validation_byte_manifest_validate(
            bundle, manifest_path, root, object_root = object_root,
            verify_git = FALSE
        )
    )
    if (!file.rename(heldout_backup, local_paths[[2L]])) {
        stop("Cannot restore held-out fixture bytes.", call. = FALSE)
    }
    validation_test_write_tsv(manifest, manifest_path)

    verified_without_root_rejected <- validation_test_rejected(
        validation_byte_manifest_validate(
            bundle, manifest_path, root, verify_git = FALSE
        )
    )
    alias_manifest <- manifest
    metadata_row <- verified_rows[[3L]]
    metadata_bytes <- readBin(
        local_paths[[3L]], what = "raw", n = file.info(local_paths[[3L]])$size
    )
    file.copy(local_paths[[1L]], local_paths[[3L]], overwrite = TRUE)
    alias_manifest$bytes[[metadata_row]] <- alias_manifest$bytes[[1L]]
    alias_manifest$sha256[[metadata_row]] <- alias_manifest$sha256[[1L]]
    validation_test_write_tsv(alias_manifest, manifest_path)
    cross_access_alias_rejected <- validation_test_rejected(
        validation_byte_manifest_validate(
            bundle, manifest_path, root, object_root = object_root,
            verify_git = FALSE
        )
    )
    connection <- file(local_paths[[3L]], open = "wb")
    writeBin(metadata_bytes, connection)
    close(connection)
    epoch_alias <- manifest
    heldout_bytes <- readBin(
        local_paths[[2L]], what = "raw", n = file.info(local_paths[[2L]])$size
    )
    file.copy(local_paths[[1L]], local_paths[[2L]], overwrite = TRUE)
    epoch_alias$bytes[[2L]] <- epoch_alias$bytes[[1L]]
    epoch_alias$sha256[[2L]] <- epoch_alias$sha256[[1L]]
    validation_test_write_tsv(epoch_alias, manifest_path)
    cross_epoch_alias_rejected <- validation_test_rejected(
        validation_byte_manifest_validate(
            bundle, manifest_path, root, object_root = object_root,
            verify_git = FALSE
        )
    )
    connection <- file(local_paths[[2L]], open = "wb")
    writeBin(heldout_bytes, connection)
    close(connection)
    validation_test_write_tsv(manifest, manifest_path)

    manifest$terminal_http_status[[1L]] <- "206"
    validation_test_write_tsv(manifest, manifest_path)
    partial_http_rejected <- validation_test_rejected(
        validation_byte_manifest_validate(
            bundle, manifest_path, root, object_root = object_root,
            verify_git = FALSE
        )
    )
    manifest$terminal_http_status[[1L]] <- "200"
    validation_test_write_tsv(manifest, manifest_path)

    extra_path <- file.path(object_root, "development", "unexpected.bin")
    connection <- file(extra_path, open = "wb")
    writeBin(charToRaw("unexpected\n"), connection)
    close(connection)
    extra_object_rejected <- validation_test_rejected(
        validation_byte_manifest_validate(
            bundle, manifest_path, root, object_root = object_root,
            verify_git = FALSE
        )
    )
    unlink(extra_path)

    fifo_object_backup <- file.path(scratch, "fifo-object.backup")
    moved <- file.rename(local_paths[[1L]], fifo_object_backup)
    fifo_status <- if (moved) {
        suppressWarnings(system2("mkfifo", shQuote(local_paths[[1L]])))
    } else 1L
    fifo_object_rejected <- moved && identical(fifo_status, 0L) &&
        validation_test_rejected(validation_byte_manifest_validate(
            bundle, manifest_path, root, object_root = object_root,
            verify_git = FALSE
        ))
    unlink(local_paths[[1L]])
    if (!file.rename(fifo_object_backup, local_paths[[1L]])) {
        stop("Cannot restore regular fixture bytes.", call. = FALSE)
    }

    expected <- validation_access_expected(valid, manifest)
    event_count <- length(expected$event_type)
    selection_commit <- unique(manifest$selection_commit)
    byte_commit <- paste(rep("1", 40L), collapse = "")
    implementation_commit <- paste(rep("2", 40L), collapse = "")
    repository_head <- rep(selection_commit, event_count)
    repository_head[
        expected$byte_freeze:(expected$implementation_freeze - 1L)
    ] <- byte_commit
    repository_head[expected$implementation_freeze:event_count] <-
        implementation_commit
    byte_sha256 <- validation_contract_sha256(manifest_path)
    byte_manifest_sha256 <- rep("none", event_count)
    byte_manifest_sha256[expected$byte_freeze:event_count] <- byte_sha256
    implementation <- rep("none", event_count)
    implementation[expected$implementation_freeze:event_count] <-
        implementation_commit
    execution_closure <- rep("none", event_count)
    execution_closure[expected$implementation_freeze:event_count] <-
        paste(rep("e", 64L), collapse = "")
    event_utc <- rep("2026-07-31T00:00:00Z", event_count)
    event_utc[expected$stream_rows] <- manifest$retrieved_utc
    access <- data.frame(
        schema_version = VALIDATION_ACCESS_SCHEMA,
        protocol_version = VALIDATION_AMENDMENT_VERSION,
        protocol_sha256 = VALIDATION_AMENDMENT_SHA256,
        event_order = as.character(seq_len(event_count)),
        event_type = expected$event_type,
        repository_head = repository_head, worktree_state = "clean",
        selection_sha256 = valid$selection_sha256,
        byte_manifest_sha256 = byte_manifest_sha256,
        implementation_commit = implementation,
        execution_closure_sha256 = execution_closure,
        execution_attempt_id = paste0("attempt_", paste(rep(
            "a", 32L
        ), collapse = "")),
        object_id = expected$object_id, event_utc = event_utc,
        bytes = expected$bytes, sha256 = expected$sha256,
        stringsAsFactors = FALSE
    )
    access_path <- file.path(scratch, "access.tsv")
    validation_test_write_tsv(access, access_path)
    access_valid <- validation_access_log_validate(
        bundle, manifest_path, access_path, root, object_root,
        verify_git = FALSE
    )

    bad_access <- access
    swap <- c(expected$execution_start, expected$heldout_rows[[1L]])
    bad_access[swap, ] <- bad_access[rev(swap), ]
    bad_access$event_order <- as.character(seq_len(nrow(bad_access)))
    validation_test_write_tsv(bad_access, access_path)
    access_order_rejected <- validation_test_rejected(
        validation_access_log_validate(
            bundle, manifest_path, access_path, root, object_root,
            verify_git = FALSE
        )
    )
    bad_access <- access
    bad_access$execution_attempt_id[[2L]] <- paste0(
        "attempt_", paste(rep("b", 32L), collapse = "")
    )
    validation_test_write_tsv(bad_access, access_path)
    access_attempt_rejected <- validation_test_rejected(
        validation_access_log_validate(
            bundle, manifest_path, access_path, root, object_root,
            verify_git = FALSE
        )
    )
    bad_access <- access
    bad_access$event_utc[[nrow(bad_access)]] <- "2026-07-30T23:59:59Z"
    validation_test_write_tsv(bad_access, access_path)
    access_time_rejected <- validation_test_rejected(
        validation_access_log_validate(
            bundle, manifest_path, access_path, root, object_root,
            verify_git = FALSE
        )
    )
    bad_access <- access[-nrow(access), , drop = FALSE]
    validation_test_write_tsv(bad_access, access_path)
    access_prefix_rejected <- validation_test_rejected(
        validation_access_log_validate(
            bundle, manifest_path, access_path, root, object_root,
            verify_git = FALSE
        )
    )

    connection <- file(local_paths[[1L]], open = "ab")
    writeBin(as.raw(1L), connection)
    close(connection)
    local_mutation_rejected <- validation_test_rejected(
        validation_byte_manifest_validate(
            bundle, manifest_path, root, object_root = object_root,
            verify_git = FALSE
        )
    )

    production <- validation_test_production_access(root, scratch)
    runtime_drift <- validation_test_production_access(
        root, scratch, runtime_drift = TRUE
    )
    selection_tree_drift <- validation_test_production_access(
        root, scratch, selection_tree_drift = TRUE
    )
    early_closure <- validation_test_production_access(
        root, scratch, early_closure = TRUE
    )
    runtime_symlink_mode <- validation_test_production_access(
        root, scratch, runtime_symlink_mode = TRUE
    )

    checks <- c(
        role_rejected = role_rejected,
        bundle_extra_rejected = bundle_extra_rejected,
        git_spaced_path_valid = git_spaced_path_valid,
        git_env_ignored = git_env_ignored,
        git_env_name_ignored = git_env_name_ignored,
        include_config_rejected = include_config_rejected,
        worktree_config_rejected = worktree_config_rejected,
        fifo_contract_rejected = fifo_contract_rejected,
        license_binding_rejected = license_binding_rejected,
        source_identity_rejected = source_identity_rejected,
        source_sentinel_rejected = source_sentinel_rejected,
        independence_sentinel_rejected = independence_sentinel_rejected,
        source_id_sentinel_rejected = source_id_sentinel_rejected,
        population_rejected = population_rejected,
        metadata_scope_rejected = metadata_scope_rejected,
        contrast_rejected = contrast_rejected,
        dependence_study_rejected = dependence_study_rejected,
        dependence_split_rejected = dependence_split_rejected,
        seed_rejected = seed_rejected,
        forbidden_rejected = forbidden_rejected,
        implicit_rownames_rejected = implicit_rownames_rejected,
        roster_rejected = roster_rejected,
        extra_condition_rejected = extra_condition_rejected,
        filtered_model_field_rejected = filtered_model_field_rejected,
        design_field_rejected = design_field_rejected,
        design_type_rejected = design_type_rejected,
        design_role_alias_rejected = design_role_alias_rejected,
        multiple_block_rejected = multiple_block_rejected,
        nuisance_pool_fragment_rejected = nuisance_pool_fragment_rejected,
        frame_unit_binding_rejected = frame_unit_binding_rejected,
        membership_pool_fragment_rejected = membership_pool_fragment_rejected,
        filter_field_rejected = filter_field_rejected,
        predicate_role_alias_rejected = predicate_role_alias_rejected,
        contradictory_filter_rejected = contradictory_filter_rejected,
        archive_member_rejected = archive_member_rejected,
        unsafe_member_rejected = unsafe_member_rejected,
        unused_object_rejected = unused_object_rejected,
        reserved_object_rejected = reserved_object_rejected,
        fragment_alias_rejected = fragment_alias_rejected,
        canonical_locator_rejected = canonical_locator_rejected,
        illegal_uri_character_rejected = illegal_uri_character_rejected,
        long_authority_rejected = long_authority_rejected,
        evidence_alias_rejected = evidence_alias_rejected,
        canonical_evidence_alias_rejected = canonical_evidence_alias_rejected,
        input_semantics_rejected = input_semantics_rejected,
        zero_missing_rejected = zero_missing_rejected,
        underflow_missing_rejected = underflow_missing_rejected,
        terminal_missing_token_rejected = terminal_missing_token_rejected,
        adapter_rejected = adapter_rejected,
        metadata_adapter_rejected = metadata_adapter_rejected,
        selector_rejected = selector_rejected,
        axis_alias_rejected = axis_alias_rejected,
        quantity_alias_rejected = quantity_alias_rejected,
        cell_alias_rejected = cell_alias_rejected,
        container_rejected = container_rejected,
        quantity_selector_rejected = quantity_selector_rejected,
        orientation_selector_rejected = orientation_selector_rejected,
        maxlfq_state_rejected = maxlfq_state_rejected,
        duplicate_operator_rejected = duplicate_operator_rejected,
        technical_sample_alias_rejected = technical_sample_alias_rejected,
        technical_model_alias_rejected = technical_model_alias_rejected,
        technical_filter_alias_rejected = technical_filter_alias_rejected,
        mapping_selector_rejected = mapping_selector_rejected,
        mapping_container_rejected = mapping_container_rejected,
        metadata_reference_alias_rejected = metadata_reference_alias_rejected,
        sample_join_filter_alias_rejected = sample_join_filter_alias_rejected,
        sample_join_model_alias_rejected = sample_join_model_alias_rejected,
        cell_join_filter_alias_rejected = cell_join_filter_alias_rejected,
        citation_binding_rejected = citation_binding_rejected,
        citation_sentinel_rejected = citation_sentinel_rejected,
        reordered_duplicate_rejected = reordered_duplicate_rejected,
        reversed_contrast_rejected = reversed_contrast_rejected,
        positive_null_overlap_rejected = positive_null_overlap_rejected,
        seed_bound_rejected = seed_bound_rejected,
        closure_rejected = closure_rejected,
        numeric_rejected = numeric_rejected,
        failure_mapping_rejected = failure_mapping_rejected,
        partial_http_rejected = partial_http_rejected,
        verified_without_root_rejected = verified_without_root_rejected,
        terminal_planned_alias_rejected = terminal_planned_alias_rejected,
        cross_access_alias_rejected = cross_access_alias_rejected,
        cross_epoch_alias_rejected = cross_epoch_alias_rejected,
        extra_object_rejected = extra_object_rejected,
        fifo_object_rejected = fifo_object_rejected,
        access_order_rejected = access_order_rejected,
        access_attempt_rejected = access_attempt_rejected,
        access_time_rejected = access_time_rejected,
        access_prefix_rejected = access_prefix_rejected,
        local_mutation_rejected = local_mutation_rejected,
        production_git_valid = production$production_git_valid,
        dirty_worktree_rejected = production$dirty_worktree_rejected,
        descendant_head_rejected = production$descendant_head_rejected,
        ambient_closure_ignored = production$ambient_closure_ignored,
        closure_hash_rejected = production$closure_hash_rejected,
        replace_ref_rejected = production$replace_ref_rejected,
        runtime_drift_rejected = runtime_drift$fault_rejected,
        selection_tree_drift_rejected = selection_tree_drift$fault_rejected,
        early_closure_rejected = early_closure$fault_rejected,
        runtime_symlink_mode_rejected = runtime_symlink_mode$fault_rejected
    )
    if (!all(checks) || byte_valid$unavailable != nrow(manifest) ||
        local_valid$verified != length(verified_rows) ||
        access_valid$events != event_count) {
        stop("Validation selection adversaries failed: ",
            paste(names(checks)[!checks], collapse = ", "), call. = FALSE)
    }
    list(
        selection_sha256 = valid$selection_sha256,
        sources = nrow(valid$sources), tasks = nrow(valid$tasks),
        null_tasks = sum(valid$tasks$task_kind == "null"),
        eligible_targets = valid$eligible_targets,
        planned_objects = nrow(valid$objects), access_events = event_count,
        adversaries = length(checks)
    )
}
