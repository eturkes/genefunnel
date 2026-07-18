# Assisted-by: OpenAI Codex.

aggregation_kang_registry_integer <- function(registry, key) {
    value <- suppressWarnings(as.integer(aggregation_registry_value(
        registry,
        "kang",
        key
    )))
    if (length(value) != 1L || is.na(value) || value < 1L) {
        stop("Kang integer registry field is invalid: ", key, call. = FALSE)
    }
    value
}

aggregation_kang_donors <- function(registry) {
    value <- c(
        aggregation_registry_vector(registry, "kang", "training_donors"),
        aggregation_registry_vector(registry, "kang", "heldout_donors")
    )
    if (length(value) != 8L || anyDuplicated(value)) {
        stop("Kang donor registry is invalid.", call. = FALSE)
    }
    value
}

aggregation_kang_validate_contract <- function(registry) {
    keys <- c(
        "barcode_join", "unit_eligibility", "stability_gate",
        "training_direction", "biological_test"
    )
    expected <- c(
        paste0(
            "concatenate_barcode_files_in_matrix_order;",
            "base_R_make.unique_sep_empty;must_identically_equal_metadata_order"
        ),
        paste0(
            "at_least_40_cells_total_and_20_each_half;",
            "remove_failed_cell_type_from_full_and_both_halves"
        ),
        "median_split_spearman>=0.70;quantile10_split_spearman>=0.50",
        "sign_of_training_donor_median_contrast;exact_zero_or_undefined_fails",
        paste0(
            "zero_contrasts_omitted;two_sided_exact_binomial_p0.5;",
            "p=1_if_no_nonzero;Holm_across_two_pathways"
        )
    )
    observed <- vapply(keys, function(key) {
        aggregation_registry_value(registry, "kang", key)
    }, character(1L))
    if (!identical(unname(observed), expected)) {
        stop("Kang runner disagrees with B-1.0.3.", call. = FALSE)
    }
    if (!Sys.getlocale("LC_COLLATE") %in% c("C", "POSIX")) {
        stop("Kang execution requires C collation.", call. = FALSE)
    }
    invisible(registry)
}

aggregation_kang_extract_tar <- function(path, directory, registry) {
    members <- aggregation_registry_vector(registry, "kang", "tar_members")
    listed <- utils::untar(path, list = TRUE)
    if (any(vapply(members, function(member) sum(listed == member) != 1L,
        logical(1L)))) {
        stop("Kang raw tar members are absent or duplicated.", call. = FALSE)
    }
    status <- utils::untar(path, files = members, exdir = directory)
    paths <- file.path(directory, members)
    if (!identical(status, 0L) || any(!file.exists(paths))) {
        stop("Kang raw tar extraction failed.", call. = FALSE)
    }
    stats::setNames(paths, c("matrix_1", "barcodes_1", "matrix_2", "barcodes_2"))
}

aggregation_kang_read_barcodes <- function(path) {
    value <- readLines(gzfile(path), warn = FALSE)
    if (!length(value) || anyNA(value) || any(!nzchar(value)) ||
        anyDuplicated(value)) {
        stop("Kang barcode file is malformed.", call. = FALSE)
    }
    value
}

aggregation_kang_join_barcodes <- function(first, second, expected = NULL) {
    raw <- c(first, second)
    joined <- make.unique(raw, sep = "")
    if (anyDuplicated(joined) || (!is.null(expected) &&
        !identical(joined, expected))) {
        stop("Kang barcode-to-metadata join failed.", call. = FALSE)
    }
    list(
        joined = joined,
        raw_duplicates = length(raw) - length(unique(raw)),
        batch_sizes = c(length(first), length(second))
    )
}

aggregation_kang_read_metadata <- function(path, joined) {
    value <- utils::read.delim(
        gzfile(path), row.names = 1L, check.names = FALSE,
        stringsAsFactors = FALSE
    )
    required <- c("ind", "stim", "cell", "multiplets")
    valid <- all(required %in% names(value)) &&
        identical(rownames(value), joined)
    if (!valid) {
        stop("Kang metadata schema or barcode alignment failed.", call. = FALSE)
    }
    value$ind <- as.character(value$ind)
    value$stim <- as.character(value$stim)
    value$cell <- as.character(value$cell)
    value$multiplets <- as.character(value$multiplets)
    value
}

aggregation_kang_read_genes <- function(path) {
    value <- utils::read.delim(
        gzfile(path), header = FALSE, check.names = FALSE,
        stringsAsFactors = FALSE, quote = "", comment.char = ""
    )
    if (ncol(value) != 2L || !nrow(value)) {
        stop("Kang gene table schema is invalid.", call. = FALSE)
    }
    names(value) <- c("ensembl_id", "gene_symbol")
    value$ensembl_id <- as.character(value$ensembl_id)
    value$gene_symbol <- as.character(value$gene_symbol)
    value$gene_symbol[is.na(value$gene_symbol)] <- ""
    if (anyNA(value$ensembl_id) || any(!nzchar(value$ensembl_id)) ||
        anyDuplicated(value$ensembl_id)) {
        stop("Kang Ensembl identifiers are invalid.", call. = FALSE)
    }
    value
}

aggregation_kang_matrix_declaration <- function(path) {
    connection <- gzfile(path, open = "rt")
    on.exit(close(connection), add = TRUE)
    header <- readLines(connection, n = 1L, warn = FALSE)
    if (!identical(
        header,
        "%%MatrixMarket matrix coordinate real general"
    )) {
        stop("Kang Matrix Market declaration is invalid.", call. = FALSE)
    }
    line <- readLines(connection, n = 1L, warn = FALSE)
    while (length(line) && startsWith(line, "%")) {
        line <- readLines(connection, n = 1L, warn = FALSE)
    }
    if (length(line) != 1L || !nzchar(trimws(line))) {
        stop("Kang Matrix Market dimensions are invalid.", call. = FALSE)
    }
    fields <- suppressWarnings(as.numeric(strsplit(
        trimws(line), "[[:space:]]+"
    )[[1L]]))
    if (length(fields) != 3L || anyNA(fields) || any(!is.finite(fields)) ||
        any(fields < 0) || any(fields != floor(fields))) {
        stop("Kang Matrix Market dimensions are invalid.", call. = FALSE)
    }
    fields
}

aggregation_kang_read_matrix <- function(path, dimensions, batch) {
    declaration <- aggregation_kang_matrix_declaration(path)
    connection <- gzfile(path, open = "rt")
    on.exit(close(connection), add = TRUE)
    value <- Matrix::readMM(connection)
    trailing <- scan(connection, what = character(), nmax = 1L, quiet = TRUE)
    if (!inherits(value, "TsparseMatrix") ||
        length(trailing) > 0L || length(value@x) != declaration[[3L]] ||
        !identical(as.numeric(dim(value)), declaration[1:2]) ||
        !identical(as.integer(dim(value)), as.integer(dimensions)) ||
        anyNA(value@x) || any(!is.finite(value@x)) || any(value@x < 0) ||
        any(value@x != floor(value@x))) {
        stop("Kang Matrix Market counts are invalid.", call. = FALSE)
    }
    entries <- length(value@x)
    canonical <- Matrix::sparseMatrix(
        i = value@i, j = value@j, x = value@x,
        dims = dim(value), index1 = FALSE
    )
    canonical <- Matrix::drop0(canonical)
    if (any(canonical@x < 0) || any(canonical@x != floor(canonical@x))) {
        stop("Kang coalesced counts are invalid.", call. = FALSE)
    }
    list(
        matrix = canonical,
        facts = data.frame(
            batch = batch, rows = nrow(canonical), columns = ncol(canonical),
            input_entries = entries, canonical_nonzero = length(canonical@x),
            coalesced_or_zero_entries = entries - length(canonical@x),
            stringsAsFactors = FALSE, check.names = FALSE
        )
    )
}

aggregation_kang_read_counts <- function(paths, genes, barcodes) {
    dimensions <- list(
        c(nrow(genes), barcodes$batch_sizes[[1L]]),
        c(nrow(genes), barcodes$batch_sizes[[2L]])
    )
    first <- aggregation_kang_read_matrix(
        paths[["matrix_1"]], dimensions[[1L]], "GSM2560248"
    )
    second <- aggregation_kang_read_matrix(
        paths[["matrix_2"]], dimensions[[2L]], "GSM2560249"
    )
    value <- cbind(first$matrix, second$matrix)
    rownames(value) <- genes$ensembl_id
    colnames(value) <- barcodes$joined
    list(matrix = value, facts = rbind(first$facts, second$facts))
}

aggregation_kang_unit_grid <- function(registry) {
    value <- expand.grid(
        cell_type = aggregation_registry_vector(registry, "kang", "cell_types"),
        condition = aggregation_registry_vector(registry, "kang", "conditions"),
        donor = aggregation_kang_donors(registry), KEEP.OUT.ATTRS = FALSE,
        stringsAsFactors = FALSE
    )
    value <- value[c("donor", "condition", "cell_type")]
    value$unit_id <- paste(value$donor, value$condition, value$cell_type, sep = "::")
    value
}

aggregation_kang_retained_cells <- function(metadata, registry) {
    complete <- stats::complete.cases(metadata[c(
        "ind", "stim", "cell", "multiplets"
    )])
    complete & metadata$ind %in% aggregation_kang_donors(registry) &
        metadata$stim %in% aggregation_registry_vector(
            registry, "kang", "conditions"
        ) & metadata$cell %in% aggregation_registry_vector(
            registry, "kang", "cell_types"
        ) & metadata$multiplets == aggregation_registry_value(
            registry, "kang", "multiplets"
        )
}

aggregation_kang_unit_one <- function(indices, barcodes, minimum) {
    ordered <- indices[order(barcodes[indices], method = "radix")]
    half <- ifelse(seq_along(ordered) %% 2L == 1L, "odd", "even")
    counts <- c(total = length(ordered), odd = sum(half == "odd"),
        even = sum(half == "even"))
    list(
        indices = ordered, half = half, counts = counts,
        eligible = counts[["total"]] >= 2L * minimum &&
            counts[["odd"]] >= minimum && counts[["even"]] >= minimum
    )
}

aggregation_kang_unit_result <- function(
    row,
    groups,
    barcodes,
    minimum
) {
    indices <- groups[[row$unit_id]]
    if (is.null(indices)) indices <- integer()
    unit <- aggregation_kang_unit_one(indices, barcodes, minimum)
    assignment <- if (unit$eligible) data.frame(
        cell_index = unit$indices, unit_id = row$unit_id,
        half = unit$half, stringsAsFactors = FALSE, check.names = FALSE
    ) else NULL
    fact <- data.frame(
        total_cells = unit$counts[["total"]],
        odd_cells = unit$counts[["odd"]],
        even_cells = unit$counts[["even"]], eligible = unit$eligible,
        reason = if (unit$eligible) "eligible" else "insufficient_cells",
        stringsAsFactors = FALSE, check.names = FALSE
    )
    list(fact = fact, assignment = assignment)
}

aggregation_kang_units <- function(metadata, registry) {
    grid <- aggregation_kang_unit_grid(registry)
    retained <- aggregation_kang_retained_cells(metadata, registry)
    key <- paste(metadata$ind, metadata$stim, metadata$cell, sep = "::")
    groups <- split(which(retained), key[retained], drop = TRUE)
    minimum <- aggregation_kang_registry_integer(
        registry, "minimum_cells_per_half"
    )
    results <- lapply(seq_len(nrow(grid)), function(index) {
        aggregation_kang_unit_result(
            grid[index, , drop = FALSE], groups, rownames(metadata), minimum
        )
    })
    manifest <- cbind(grid, do.call(rbind, lapply(results, `[[`, "fact")))
    rownames(manifest) <- NULL
    assignments <- lapply(results, `[[`, "assignment")
    assignments <- assignments[!vapply(assignments, is.null, logical(1L))]
    assignment <- if (length(assignments)) do.call(rbind, assignments) else {
        data.frame(cell_index = integer(), unit_id = character(), half = character())
    }
    rownames(assignment) <- NULL
    list(manifest = manifest, assignments = assignment)
}

aggregation_kang_design_columns <- function(manifest) {
    eligible <- manifest[manifest$eligible, , drop = FALSE]
    grid <- expand.grid(
        view = c("full", "odd", "even"), unit_id = eligible$unit_id,
        KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
    )
    grid$column_id <- paste(grid$unit_id, grid$view, sep = "::")
    grid
}

aggregation_kang_cell_design <- function(units, cell_count) {
    columns <- aggregation_kang_design_columns(units$manifest)
    assignment <- units$assignments
    full <- paste(assignment$unit_id, "full", sep = "::")
    half <- paste(assignment$unit_id, assignment$half, sep = "::")
    value <- if (nrow(columns)) Matrix::sparseMatrix(
        i = rep(assignment$cell_index, 2L),
        j = match(c(full, half), columns$column_id), x = 1,
        dims = c(cell_count, nrow(columns))
    ) else Matrix::sparseMatrix(
        i = integer(), j = integer(), x = numeric(),
        dims = c(cell_count, 0L)
    )
    colnames(value) <- columns$column_id
    observed <- as.numeric(Matrix::colSums(value))
    manifest <- units$manifest[match(columns$unit_id, units$manifest$unit_id), ]
    expected <- ifelse(
        columns$view == "full", manifest$total_cells,
        ifelse(columns$view == "odd", manifest$odd_cells, manifest$even_cells)
    )
    if (anyNA(value@x) || !identical(observed, as.numeric(expected))) {
        stop("Kang cell aggregation design failed.", call. = FALSE)
    }
    list(matrix = value, columns = columns)
}

aggregation_kang_collapse_symbols <- function(counts, design, genes) {
    gene_sums <- counts %*% design$matrix
    retained <- nzchar(genes$gene_symbol)
    symbols <- unique(genes$gene_symbol[retained])
    mapping <- Matrix::sparseMatrix(
        i = match(genes$gene_symbol[retained], symbols),
        j = which(retained), x = 1,
        dims = c(length(symbols), nrow(counts))
    )
    value <- Matrix::drop0(mapping %*% gene_sums)
    rownames(value) <- symbols
    colnames(value) <- design$columns$column_id
    if (anyNA(value@x) || any(!is.finite(value@x)) || any(value@x < 0) ||
        any(value@x != floor(value@x))) {
        stop("Kang gene-symbol collapse failed.", call. = FALSE)
    }
    list(
        matrix = value,
        facts = data.frame(
            gene_rows = nrow(genes), nonempty_gene_rows = sum(retained),
            unique_symbols = length(symbols),
            duplicate_symbol_rows = sum(retained) - length(symbols),
            aggregate_columns = ncol(value),
            collapsed_nonzero = length(value@x),
            stringsAsFactors = FALSE, check.names = FALSE
        )
    )
}

aggregation_kang_pathway_records <- function(lines) {
    fields <- strsplit(lines, "\t", fixed = TRUE)
    if (!length(fields) || any(lengths(fields) < 3L)) {
        stop("Reactome GMT records are malformed.", call. = FALSE)
    }
    ids <- vapply(fields, `[[`, character(1L), 2L)
    selected <- startsWith(ids, "R-HSA-")
    fields <- fields[selected]
    ids <- ids[selected]
    if (!length(fields) || anyDuplicated(ids)) {
        stop("Reactome human pathway identifiers are invalid.", call. = FALSE)
    }
    list(fields = fields, ids = ids)
}

aggregation_kang_parse_pathways <- function(lines, measured_symbols, registry) {
    records <- aggregation_kang_pathway_records(lines)
    limits <- suppressWarnings(as.integer(strsplit(
        aggregation_registry_value(registry, "kang", "pathway_size"),
        ":", fixed = TRUE
    )[[1L]]))
    if (length(limits) != 2L || anyNA(limits)) {
        stop("Kang pathway-size registry is invalid.", call. = FALSE)
    }
    retained <- lapply(records$fields, function(field) {
        members <- unique(field[-c(1L, 2L)])
        members[nzchar(members) & members %in% measured_symbols]
    })
    raw_count <- vapply(records$fields, function(field) {
        length(unique(field[-c(1L, 2L)][nzchar(field[-c(1L, 2L)])]))
    }, integer(1L))
    retained_count <- lengths(retained)
    included <- retained_count >= limits[[1L]] & retained_count <= limits[[2L]]
    manifest <- data.frame(
        pathway_id = records$ids,
        pathway_name = vapply(records$fields, `[[`, character(1L), 1L),
        raw_member_count = raw_count, retained_member_count = retained_count,
        included = included, stringsAsFactors = FALSE, check.names = FALSE
    )
    sets <- stats::setNames(retained[included], records$ids[included])
    rows <- lapply(names(sets), function(id) data.frame(
        pathway_id = id, member_rank = seq_along(sets[[id]]),
        gene_symbol = sets[[id]], stringsAsFactors = FALSE,
        check.names = FALSE
    ))
    membership <- if (length(rows)) do.call(rbind, rows) else data.frame(
        pathway_id = character(), member_rank = integer(),
        gene_symbol = character(), stringsAsFactors = FALSE,
        check.names = FALSE
    )
    rownames(membership) <- NULL
    list(manifest = manifest, membership = membership, sets = sets)
}

aggregation_kang_read_pathways <- function(path, measured_symbols, registry) {
    member <- "ReactomePathways.gmt"
    listed <- utils::unzip(path, list = TRUE)
    if (sum(listed$Name == member) != 1L) {
        stop("Reactome archive member is absent or duplicated.", call. = FALSE)
    }
    connection <- unz(path, member, open = "rt")
    on.exit(close(connection), add = TRUE)
    aggregation_kang_parse_pathways(
        readLines(connection, warn = FALSE), measured_symbols, registry
    )
}

aggregation_kang_view_profiles <- function(collapsed, units, view) {
    manifest <- units$manifest[units$manifest$eligible, , drop = FALSE]
    columns <- paste(manifest$unit_id, view, sep = "::")
    counts <- switch(
        view,
        full = manifest$total_cells,
        odd = manifest$odd_cells,
        even = manifest$even_cells,
        stop("Unknown Kang audit view.", call. = FALSE)
    )
    value <- as.matrix(collapsed[, columns, drop = FALSE])
    value <- sweep(value, 2L, counts, "/")
    colnames(value) <- manifest$unit_id
    group <- paste(manifest$donor, manifest$condition, sep = "::")
    list(
        matrix = value, groups = stats::setNames(group, manifest$unit_id),
        weights = stats::setNames(as.numeric(counts), manifest$unit_id),
        manifest = manifest
    )
}

aggregation_kang_audit_view <- function(
    collapsed,
    units,
    catalogue,
    view,
    audit
) {
    profiles <- aggregation_kang_view_profiles(collapsed, units, view)
    value <- audit(
        profiles$matrix, catalogue$sets, profiles$groups, profiles$weights,
        missing = "reject"
    )$summary
    map <- unique(data.frame(
        group = unname(profiles$groups), donor = profiles$manifest$donor,
        condition = profiles$manifest$condition, stringsAsFactors = FALSE
    ))
    map <- map[match(value$group, map$group), , drop = FALSE]
    data.frame(
        view = view, donor = map$donor, condition = map$condition,
        pathway_id = value$gene_set, eligible = value$eligible,
        reason = value$reason, aggregate_score = value$aggregate_score,
        weighted_unit_score = value$weighted_unit_score,
        aggregation_gap = value$aggregation_gap,
        normalized_gap = value$normalized_gap,
        normalized_status = value$normalized_status,
        retained_size = value$retained_size,
        active_cell_types = value$active_unit_count,
        stringsAsFactors = FALSE, check.names = FALSE
    )
}

aggregation_kang_audits <- function(collapsed, units, catalogue, audit) {
    if (!any(units$manifest$eligible)) {
        return(data.frame(
            view = character(), donor = character(), condition = character(),
            pathway_id = character(), eligible = logical(), reason = character(),
            aggregate_score = numeric(), weighted_unit_score = numeric(),
            aggregation_gap = numeric(), normalized_gap = numeric(),
            normalized_status = character(), retained_size = integer(),
            active_cell_types = integer(), stringsAsFactors = FALSE,
            check.names = FALSE
        ))
    }
    result <- do.call(rbind, lapply(c("full", "odd", "even"), function(view) {
        aggregation_kang_audit_view(
            collapsed, units, catalogue, view, audit
        )
    }))
    rownames(result) <- NULL
    result
}

aggregation_kang_external_inputs <- function(files, registry, directory) {
    paths <- stats::setNames(files$path, files$data_id)
    extracted <- aggregation_kang_extract_tar(
        paths[["kang_batch2_counts"]], directory, registry
    )
    first <- aggregation_kang_read_barcodes(extracted[["barcodes_1"]])
    second <- aggregation_kang_read_barcodes(extracted[["barcodes_2"]])
    raw_metadata <- utils::read.delim(
        gzfile(paths[["kang_batch2_metadata"]]), row.names = 1L,
        check.names = FALSE, stringsAsFactors = FALSE
    )
    barcodes <- aggregation_kang_join_barcodes(
        first, second, rownames(raw_metadata)
    )
    metadata <- aggregation_kang_read_metadata(
        paths[["kang_batch2_metadata"]], barcodes$joined
    )
    genes <- aggregation_kang_read_genes(paths[["kang_batch2_genes"]])
    counts <- aggregation_kang_read_counts(extracted, genes, barcodes)
    symbols <- unique(genes$gene_symbol[nzchar(genes$gene_symbol)])
    catalogue <- aggregation_kang_read_pathways(
        paths[["reactome_v97"]], symbols, registry
    )
    list(
        counts = counts$matrix, matrix_facts = counts$facts,
        genes = genes, metadata = metadata, barcodes = barcodes,
        catalogue = catalogue
    )
}

aggregation_kang_validate_inputs <- function(inputs) {
    counts <- inputs$counts
    valid <- (is.matrix(counts) || inherits(counts, "dMatrix")) &&
        nrow(counts) == nrow(inputs$genes) &&
        ncol(counts) == nrow(inputs$metadata) &&
        identical(rownames(counts), inputs$genes$ensembl_id) &&
        identical(colnames(counts), rownames(inputs$metadata)) &&
        is.list(inputs$catalogue$sets) && length(inputs$catalogue$sets) > 0L &&
        !is.null(names(inputs$catalogue$sets)) &&
        !anyDuplicated(names(inputs$catalogue$sets))
    if (!valid) {
        stop("Kang analysis inputs are malformed or misaligned.", call. = FALSE)
    }
    invisible(inputs)
}

aggregation_kang_analyze <- function(inputs, registry, audit) {
    aggregation_kang_validate_contract(registry)
    aggregation_kang_validate_inputs(inputs)
    units <- aggregation_kang_units(inputs$metadata, registry)
    design <- aggregation_kang_cell_design(units, nrow(inputs$metadata))
    collapsed <- aggregation_kang_collapse_symbols(
        inputs$counts, design, inputs$genes
    )
    audits <- aggregation_kang_audits(
        collapsed$matrix, units, inputs$catalogue, audit
    )
    decision <- aggregation_kang_summarize(
        audits, units, inputs$catalogue, registry
    )
    c(list(
        units = units, design = design, collapsed = collapsed,
        audits = audits
    ), decision)
}
