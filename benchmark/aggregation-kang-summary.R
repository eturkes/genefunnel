# Assisted-by: OpenAI Codex.

aggregation_kang_quantile <- function(values, probability) {
    values <- values[is.finite(values)]
    if (!length(values)) return(NA_real_)
    unname(stats::quantile(
        values, probs = probability, names = FALSE, type = 8L
    ))
}

aggregation_kang_spearman <- function(first, second) {
    retained <- is.finite(first) & is.finite(second)
    first <- first[retained]
    second <- second[retained]
    if (length(first) < 3L || length(unique(first)) < 2L ||
        length(unique(second)) < 2L) {
        return(NA_real_)
    }
    unname(stats::cor(first, second, method = "spearman"))
}

aggregation_kang_audit_vector <- function(
    audits,
    donor,
    condition,
    view,
    pathways
) {
    selected <- audits$donor == donor & audits$condition == condition &
        audits$view == view
    value <- audits[selected, , drop = FALSE]
    if (anyDuplicated(value$pathway_id)) {
        stop("Kang audit keys are duplicated.", call. = FALSE)
    }
    index <- match(pathways, value$pathway_id)
    present <- !is.na(index)
    normalized <- rep.int(NA_real_, length(pathways))
    defined <- rep.int(FALSE, length(pathways))
    normalized[present] <- value$normalized_gap[index[present]]
    defined[present] <- value$eligible[index[present]] &
        value$normalized_status[index[present]] == "defined" &
        is.finite(normalized[present])
    list(normalized = normalized, defined = defined)
}

aggregation_kang_stability_one <- function(
    audits,
    donor,
    condition,
    pathways
) {
    odd <- aggregation_kang_audit_vector(
        audits, donor, condition, "odd", pathways
    )
    even <- aggregation_kang_audit_vector(
        audits, donor, condition, "even", pathways
    )
    common <- odd$defined & even$defined
    correlation <- aggregation_kang_spearman(
        odd$normalized[common], even$normalized[common]
    )
    data.frame(
        donor = donor, condition = condition,
        common_defined_pathways = sum(common), split_spearman = correlation,
        defined = is.finite(correlation), stringsAsFactors = FALSE,
        check.names = FALSE
    )
}

aggregation_kang_stability <- function(audits, catalogue, registry) {
    grid <- expand.grid(
        condition = aggregation_registry_vector(registry, "kang", "conditions"),
        donor = aggregation_kang_donors(registry), KEEP.OUT.ATTRS = FALSE,
        stringsAsFactors = FALSE
    )
    rows <- lapply(seq_len(nrow(grid)), function(index) {
        aggregation_kang_stability_one(
            audits, grid$donor[[index]], grid$condition[[index]],
            names(catalogue$sets)
        )
    })
    result <- do.call(rbind, rows)
    rownames(result) <- NULL
    result
}

aggregation_kang_contrast_one <- function(
    audits,
    pathway,
    donor,
    registry
) {
    conditions <- aggregation_registry_vector(registry, "kang", "conditions")
    values <- lapply(conditions, function(condition) {
        aggregation_kang_audit_vector(
            audits, donor, condition, "full", pathway
        )
    })
    defined <- vapply(values, function(value) value$defined[[1L]], logical(1L))
    normalized <- vapply(values, function(value) {
        value$normalized[[1L]]
    }, numeric(1L))
    contrast <- if (all(defined)) normalized[[2L]] - normalized[[1L]] else NA
    data.frame(
        pathway_id = pathway, donor = donor,
        cohort = if (donor %in% aggregation_registry_vector(
            registry, "kang", "training_donors"
        )) "training" else "heldout",
        ctrl_normalized_gap = normalized[[1L]],
        stim_normalized_gap = normalized[[2L]], defined = all(defined),
        contrast = contrast, stringsAsFactors = FALSE, check.names = FALSE
    )
}

aggregation_kang_contrasts <- function(audits, registry) {
    pathways <- aggregation_registry_vector(registry, "kang", "primary_pathways")
    donors <- aggregation_kang_donors(registry)
    grid <- expand.grid(
        donor = donors, pathway_id = pathways, KEEP.OUT.ATTRS = FALSE,
        stringsAsFactors = FALSE
    )
    rows <- lapply(seq_len(nrow(grid)), function(index) {
        aggregation_kang_contrast_one(
            audits, grid$pathway_id[[index]], grid$donor[[index]], registry
        )
    })
    result <- do.call(rbind, rows)
    rownames(result) <- NULL
    result
}

aggregation_kang_exact_sign_p <- function(values) {
    values <- values[is.finite(values) & values != 0]
    if (!length(values)) return(1)
    stats::binom.test(
        sum(values > 0), length(values), p = 0.5, alternative = "two.sided"
    )$p.value
}

aggregation_kang_pathway_decision <- function(contrasts, pathway) {
    value <- contrasts[contrasts$pathway_id == pathway, , drop = FALSE]
    training <- value[value$cohort == "training", , drop = FALSE]
    heldout <- value[value$cohort == "heldout", , drop = FALSE]
    training_complete <- nrow(training) == 4L && all(training$defined)
    training_median <- if (training_complete) {
        stats::median(training$contrast)
    } else NA_real_
    direction <- if (is.finite(training_median) && training_median != 0) {
        sign(training_median)
    } else NA_real_
    matching <- if (is.finite(direction)) {
        sum(heldout$defined & heldout$contrast != 0 &
            sign(heldout$contrast) == direction)
    } else 0L
    heldout_complete <- nrow(heldout) == 4L && all(heldout$defined)
    data.frame(
        pathway_id = pathway, training_defined = sum(training$defined),
        training_median_contrast = training_median,
        training_direction = if (is.na(direction)) "undefined" else if (
            direction > 0
        ) "positive" else "negative",
        heldout_defined = sum(heldout$defined), heldout_matching = matching,
        heldout_pass = training_complete && heldout_complete && matching >= 3L,
        all_donors_defined = nrow(value) == 8L && all(value$defined),
        nonzero_donors = sum(value$defined & value$contrast != 0),
        positive_donors = sum(value$defined & value$contrast > 0),
        exact_sign_p = aggregation_kang_exact_sign_p(value$contrast),
        stringsAsFactors = FALSE, check.names = FALSE
    )
}

aggregation_kang_pathway_decisions <- function(contrasts, registry) {
    pathways <- aggregation_registry_vector(registry, "kang", "primary_pathways")
    result <- do.call(rbind, lapply(pathways, function(pathway) {
        aggregation_kang_pathway_decision(contrasts, pathway)
    }))
    result$holm_adjusted_p <- stats::p.adjust(result$exact_sign_p, method = "holm")
    result$biological_pass <- result$heldout_pass &
        result$all_donors_defined & result$holm_adjusted_p <= 0.05
    rownames(result) <- NULL
    result
}

aggregation_kang_endpoint <- function(
    gate,
    endpoint,
    stratum,
    estimate,
    comparison,
    threshold,
    count,
    complete
) {
    complete <- isTRUE(complete)
    passed <- complete && is.finite(estimate) && switch(
        comparison,
        `<=` = estimate <= threshold,
        `>=` = estimate >= threshold,
        stop("Kang endpoint comparison is invalid.", call. = FALSE)
    )
    data.frame(
        gate = gate, endpoint = endpoint, stratum = stratum,
        estimate = estimate, comparison = comparison, threshold = threshold,
        count = as.integer(count), complete = complete, passed = passed,
        stringsAsFactors = FALSE, check.names = FALSE
    )
}

aggregation_kang_stability_endpoints <- function(stability) {
    complete <- nrow(stability) == 16L && all(stability$defined)
    values <- stability$split_spearman
    rbind(
        aggregation_kang_endpoint(
            "technical_stability", "median_split_spearman",
            "all_donor_conditions", aggregation_kang_quantile(values, 0.5),
            ">=", 0.70, nrow(stability), complete
        ),
        aggregation_kang_endpoint(
            "technical_stability", "quantile10_split_spearman",
            "all_donor_conditions", aggregation_kang_quantile(values, 0.1),
            ">=", 0.50, nrow(stability), complete
        )
    )
}

aggregation_kang_decision_endpoints <- function(decisions) {
    heldout <- do.call(rbind, lapply(seq_len(nrow(decisions)), function(index) {
        value <- decisions[index, , drop = FALSE]
        aggregation_kang_endpoint(
            "heldout_replication", "matching_heldout_donors",
            value$pathway_id, value$heldout_matching, ">=", 3, 4L,
            value$training_direction != "undefined" && value$heldout_defined == 4L
        )
    }))
    biological <- do.call(rbind, lapply(seq_len(nrow(decisions)), function(index) {
        value <- decisions[index, , drop = FALSE]
        aggregation_kang_endpoint(
            "biological_effect", "holm_adjusted_sign_p", value$pathway_id,
            value$holm_adjusted_p, "<=", 0.05, value$nonzero_donors,
            value$all_donors_defined && value$heldout_pass
        )
    }))
    rbind(heldout, biological)
}

aggregation_kang_endpoints <- function(stability, decisions) {
    result <- rbind(
        aggregation_kang_stability_endpoints(stability),
        aggregation_kang_decision_endpoints(decisions)
    )
    rownames(result) <- NULL
    result
}

aggregation_kang_summary <- function(
    audits,
    units,
    catalogue,
    stability,
    decisions,
    endpoints
) {
    technical <- endpoints$gate == "technical_stability"
    heldout <- endpoints$gate == "heldout_replication"
    biological <- endpoints$gate == "biological_effect"
    data.frame(
        protocol_version = AGGREGATION_PROTOCOL_VERSION,
        fixed_units = nrow(units$manifest), eligible_units = sum(units$manifest$eligible),
        retained_cells = nrow(units$assignments), retained_pathways = length(catalogue$sets),
        audit_rows = nrow(audits), defined_audits = sum(
            audits$eligible & audits$normalized_status == "defined"
        ), stable_donor_conditions = sum(stability$defined),
        technical_stability_pass = all(endpoints$passed[technical]),
        heldout_replication_pass = all(endpoints$passed[heldout]),
        biological_effect_pass = all(endpoints$passed[biological]),
        endpoint_count = nrow(endpoints), failed_endpoints = sum(!endpoints$passed),
        all_kang_gates_pass = all(endpoints$passed), stringsAsFactors = FALSE,
        check.names = FALSE
    )
}

aggregation_kang_summarize <- function(audits, units, catalogue, registry) {
    stability <- aggregation_kang_stability(audits, catalogue, registry)
    contrasts <- aggregation_kang_contrasts(audits, registry)
    decisions <- aggregation_kang_pathway_decisions(contrasts, registry)
    endpoints <- aggregation_kang_endpoints(stability, decisions)
    list(
        stability = stability, contrasts = contrasts, decisions = decisions,
        endpoints = endpoints, summary = aggregation_kang_summary(
            audits, units, catalogue, stability, decisions, endpoints
        )
    )
}

aggregation_kang_smoke_genes <- function() {
    symbols <- sprintf("G%03d", seq_len(160L))
    symbols[[159L]] <- ""
    symbols[[160L]] <- "G096"
    data.frame(
        ensembl_id = sprintf("ENSG%011d", seq_len(160L)),
        gene_symbol = symbols,
        stringsAsFactors = FALSE, check.names = FALSE
    )
}

aggregation_kang_smoke_catalogue <- function(genes, registry) {
    primary <- aggregation_registry_vector(registry, "kang", "primary_pathways")
    ids <- c(primary, sprintf("R-HSA-990%03d", seq_len(10L)))
    lines <- vapply(seq_along(ids), function(index) {
        members <- genes$gene_symbol[
            seq.int(8L * index - 7L, 8L * index)
        ]
        paste(c(paste("Smoke pathway", index), ids[[index]], members),
            collapse = "\t")
    }, character(1L))
    aggregation_kang_parse_pathways(lines, genes$gene_symbol, registry)
}

aggregation_kang_smoke_metadata <- function(registry) {
    grid <- expand.grid(
        cell_number = seq_len(40L),
        cell = aggregation_registry_vector(registry, "kang", "cell_types")[1:2],
        stim = aggregation_registry_vector(registry, "kang", "conditions"),
        ind = aggregation_kang_donors(registry), KEEP.OUT.ATTRS = FALSE,
        stringsAsFactors = FALSE
    )
    value <- data.frame(
        ind = grid$ind, stim = grid$stim, cell = grid$cell,
        multiplets = aggregation_registry_value(registry, "kang", "multiplets"),
        stringsAsFactors = FALSE, check.names = FALSE
    )
    rownames(value) <- sprintf(
        "BC_%s_%s_%s_%03d", grid$ind, grid$stim,
        match(grid$cell, unique(grid$cell)), grid$cell_number
    )
    value
}

aggregation_kang_smoke_counts <- function(genes, metadata, catalogue, registry) {
    value <- matrix(
        3, nrow = nrow(genes), ncol = nrow(metadata),
        dimnames = list(genes$ensembl_id, rownames(metadata))
    )
    first_type <- metadata$cell == aggregation_registry_vector(
        registry, "kang", "cell_types"
    )[[1L]]
    stimulated <- metadata$stim == "stim"
    donor_rank <- match(metadata$ind, aggregation_kang_donors(registry))
    for (index in seq_along(catalogue$sets)) {
        rows <- match(catalogue$sets[[index]], genes$gene_symbol)
        high <- 160 + index
        low <- 10 + index
        delta <- rep.int(12 * index, nrow(metadata))
        if (index <= 2L) {
            delta[stimulated] <- 120 + donor_rank[stimulated]
        }
        value[rows[1:4], first_type] <- high
        value[rows[5:8], first_type] <- low
        second <- !first_type
        value[rows[1:4], second] <- rep(high - delta[second], each = 4L)
        value[rows[5:8], second] <- rep(low + delta[second], each = 4L)
    }
    Matrix::Matrix(value, sparse = TRUE)
}

aggregation_kang_smoke_fixture <- function(registry) {
    genes <- aggregation_kang_smoke_genes()
    metadata <- aggregation_kang_smoke_metadata(registry)
    catalogue <- aggregation_kang_smoke_catalogue(genes, registry)
    list(
        counts = aggregation_kang_smoke_counts(
            genes, metadata, catalogue, registry
        ),
        genes = genes, metadata = metadata, catalogue = catalogue
    )
}

aggregation_kang_expect_error <- function(expression, pattern) {
    message <- tryCatch(
        {
            force(expression)
            NA_character_
        },
        error = conditionMessage
    )
    if (is.na(message) || !grepl(pattern, message, fixed = TRUE)) {
        stop("Kang malformed-input smoke did not fail closed.", call. = FALSE)
    }
    invisible(message)
}

aggregation_kang_write_matrix_smoke <- function(
    path,
    negative = FALSE,
    header = "real",
    trailing = FALSE
) {
    connection <- gzfile(path, open = "wt")
    on.exit(close(connection), add = TRUE)
    entry <- if (negative) "1 1 -2" else "1 1 2"
    lines <- c(
        paste("%%MatrixMarket matrix coordinate", header, "general"), "% smoke",
        "2 2 3", entry, "1 1 3", "2 2 0"
    )
    if (trailing) lines <- c(lines, "1 2 4")
    writeLines(lines, connection, useBytes = TRUE)
    invisible(path)
}

aggregation_kang_validate_matrix_smoke <- function(directory) {
    path <- file.path(directory, "duplicate.mtx.gz")
    aggregation_kang_write_matrix_smoke(path)
    parsed <- aggregation_kang_read_matrix(path, c(2L, 2L), "smoke")
    stopifnot(
        parsed$matrix[1L, 1L] == 5,
        parsed$facts$input_entries == 3L,
        parsed$facts$canonical_nonzero == 1L
    )
    negative <- file.path(directory, "negative.mtx.gz")
    aggregation_kang_write_matrix_smoke(negative, negative = TRUE)
    aggregation_kang_expect_error(
        aggregation_kang_read_matrix(negative, c(2L, 2L), "negative"),
        "Matrix Market counts are invalid"
    )
    integer <- file.path(directory, "integer.mtx.gz")
    aggregation_kang_write_matrix_smoke(integer, header = "integer")
    aggregation_kang_expect_error(
        aggregation_kang_read_matrix(integer, c(2L, 2L), "integer"),
        "Matrix Market declaration is invalid"
    )
    trailing <- file.path(directory, "trailing.mtx.gz")
    aggregation_kang_write_matrix_smoke(trailing, trailing = TRUE)
    aggregation_kang_expect_error(
        aggregation_kang_read_matrix(trailing, c(2L, 2L), "trailing"),
        "Matrix Market counts are invalid"
    )
    invisible(TRUE)
}

aggregation_kang_validate_parser_smoke <- function(directory) {
    joined <- aggregation_kang_join_barcodes(c("A", "B"), c("A", "C"))
    stopifnot(
        identical(joined$joined, c("A", "B", "A1", "C")),
        identical(joined$raw_duplicates, 1L)
    )
    aggregation_kang_expect_error(
        aggregation_kang_join_barcodes(c("A"), c("B"), c("B", "A")),
        "barcode-to-metadata join failed"
    )
    aggregation_kang_validate_matrix_smoke(directory)
    invisible(TRUE)
}

aggregation_kang_validate_result_smoke <- function(result, registry) {
    stopifnot(
        nrow(result$units$manifest) == 96L,
        sum(result$units$manifest$eligible) == 32L,
        nrow(result$units$assignments) == 1280L,
        nrow(result$design$columns) == 96L,
        result$collapsed$facts$nonempty_gene_rows == 159L,
        result$collapsed$facts$unique_symbols == 158L,
        result$collapsed$facts$duplicate_symbol_rows == 1L,
        match("G096", rownames(result$collapsed$matrix)) == 96L,
        nrow(result$audits) == 576L,
        all(result$audits$normalized_status == "defined"),
        nrow(result$stability) == 16L,
        all(result$stability$defined),
        nrow(result$contrasts) == 16L,
        nrow(result$decisions) == 2L,
        all(result$decisions$exact_sign_p == 0.0078125),
        all(result$decisions$holm_adjusted_p == 0.015625),
        nrow(result$endpoints) == 6L,
        all(result$endpoints$passed),
        isTRUE(result$summary$all_kang_gates_pass)
    )
    broken <- result$audits[!(
        result$audits$view == "odd" & result$audits$donor ==
            aggregation_kang_donors(registry)[[1L]] &
            result$audits$condition == "ctrl"
    ), , drop = FALSE]
    failed <- aggregation_kang_summarize(
        broken, result$units, list(sets = stats::setNames(
            vector("list", 12L), unique(result$audits$pathway_id)
        )), registry
    )
    stopifnot(!failed$stability$defined[[1L]], !failed$endpoints$passed[[1L]])
    invisible(result)
}

aggregation_kang_validate_decision_smoke <- function(result, registry) {
    contrasts <- result$contrasts
    pathway <- aggregation_registry_vector(
        registry, "kang", "primary_pathways"
    )[[1L]]
    selected <- contrasts$pathway_id == pathway &
        contrasts$cohort == "training"
    contrasts$contrast[selected] <- c(-0.2, -0.1, 0.1, 0.2)
    decisions <- aggregation_kang_pathway_decisions(contrasts, registry)
    first <- decisions$pathway_id == pathway
    stopifnot(
        decisions$training_median_contrast[first] == 0,
        decisions$training_direction[first] == "undefined",
        !decisions$heldout_pass[first],
        !decisions$biological_pass[first]
    )
    invisible(decisions)
}

aggregation_kang_validate_empty_smoke <- function(fixture, registry, audit) {
    retained <- seq_len(16L)
    fixture$metadata <- fixture$metadata[retained, , drop = FALSE]
    fixture$counts <- fixture$counts[, retained, drop = FALSE]
    result <- aggregation_kang_analyze(fixture, registry, audit)
    stopifnot(
        !any(result$units$manifest$eligible),
        nrow(result$audits) == 0L,
        nrow(result$stability) == 16L,
        !any(result$stability$defined),
        nrow(result$endpoints) == 6L,
        !any(result$endpoints$passed),
        !result$summary$all_kang_gates_pass
    )
    invisible(result)
}

aggregation_validate_kang_smoke <- function(registry, audit) {
    old_locale <- Sys.getlocale("LC_COLLATE")
    on.exit(suppressWarnings(Sys.setlocale("LC_COLLATE", old_locale)), add = TRUE)
    if (is.na(suppressWarnings(Sys.setlocale("LC_COLLATE", "C")))) {
        stop("Kang smoke cannot set C collation.", call. = FALSE)
    }
    directory <- tempfile("aggregation-kang-smoke-")
    dir.create(directory)
    on.exit(unlink(directory, recursive = TRUE, force = TRUE), add = TRUE)
    aggregation_kang_validate_parser_smoke(directory)
    fixture <- aggregation_kang_smoke_fixture(registry)
    result <- aggregation_kang_analyze(fixture, registry, audit)
    misaligned <- fixture
    misaligned$metadata <- misaligned$metadata[rev(seq_len(nrow(
        misaligned$metadata
    ))), , drop = FALSE]
    aggregation_kang_expect_error(
        aggregation_kang_analyze(misaligned, registry, audit),
        "analysis inputs are malformed or misaligned"
    )
    aggregation_kang_validate_empty_smoke(fixture, registry, audit)
    aggregation_kang_validate_result_smoke(result, registry)
    aggregation_kang_validate_decision_smoke(result, registry)
}
