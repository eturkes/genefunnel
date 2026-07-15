# Assisted-by: OpenAI Codex.

benchmark_assert_count <- function(value, label, minimum = 1L) {
    valid <- is.numeric(value) && length(value) == 1L && !is.na(value) &&
        is.finite(value) && value == floor(value) && value >= minimum
    if (!valid) {
        stop(label, " must be an integer >= ", minimum, ".", call. = FALSE)
    }
    as.integer(value)
}

benchmark_assert_fraction <- function(value, label, allow_na = FALSE) {
    if (allow_na && length(value) == 1L && is.na(value)) {
        return(NA_real_)
    }
    valid <- is.numeric(value) && length(value) == 1L && !is.na(value) &&
        is.finite(value) && value >= 0 && value <= 1
    if (!valid) {
        stop(label, " must be between zero and one.", call. = FALSE)
    }
    as.double(value)
}

benchmark_scenario_row <- function(
    preset,
    storage,
    overlap,
    backend,
    backend_workers,
    n_features,
    n_samples,
    n_sets,
    set_size,
    matrix_seed,
    set_seed,
    density = NA_real_,
    zero_fraction = NA_real_,
    cell_missing_fraction = NA_real_,
    stored_missing_fraction = NA_real_
) {
    density_label <- if (is.na(density)) "dense" else {
        paste0("d", formatC(1000 * density, width = 3L, flag = "0"))
    }
    fixture_id <- paste(storage, overlap, density_label, sep = "-")
    data.frame(
        scenario_id = paste(fixture_id, backend, sep = "-"),
        fixture_id = fixture_id,
        preset = preset,
        storage = storage,
        overlap = overlap,
        backend = backend,
        backend_workers = as.integer(if (backend == "serial") 1L else backend_workers),
        n_features = as.integer(n_features),
        n_samples = as.integer(n_samples),
        n_sets = as.integer(n_sets),
        set_size = as.integer(set_size),
        matrix_seed = as.integer(matrix_seed),
        set_seed = as.integer(set_seed),
        density = as.double(density),
        zero_fraction = as.double(zero_fraction),
        cell_missing_fraction = as.double(cell_missing_fraction),
        stored_missing_fraction = as.double(stored_missing_fraction),
        stringsAsFactors = FALSE,
        check.names = FALSE
    )
}

benchmark_scenarios <- function(preset = "full", workers = 2L) {
    workers <- benchmark_assert_count(workers, "`workers`")
    if (!preset %in% c("smoke", "full", "hotspot")) {
        stop("`preset` must be smoke, full, or hotspot.", call. = FALSE)
    }

    if (preset == "hotspot") {
        rows <- list()
        for (density in c(0.001, 0.01, 0.03, 0.1)) {
            for (overlap in c("low", "high")) {
                rows[[length(rows) + 1L]] <- benchmark_scenario_row(
                    preset = preset,
                    storage = "sparse",
                    overlap = overlap,
                    backend = "serial",
                    backend_workers = 1L,
                    n_features = 20000L,
                    n_samples = 200L,
                    n_sets = 1000L,
                    set_size = 20L,
                    matrix_seed = 1731L,
                    set_seed = if (overlap == "low") 1732L else 1733L,
                    density = density,
                    stored_missing_fraction = 0
                )
            }
        }
        return(do.call(rbind, rows))
    }

    dimensions <- if (preset == "full") {
        list(
            dense = c(features = 20000L, samples = 60L, sets = 1000L, size = 20L),
            sparse = c(features = 20000L, samples = 600L, sets = 1000L, size = 20L)
        )
    } else {
        list(
            dense = c(features = 400L, samples = 8L, sets = 20L, size = 8L),
            sparse = c(features = 400L, samples = 16L, sets = 20L, size = 8L)
        )
    }

    overlaps <- if (preset == "full") c("low", "high") else c("low", "high")
    backends <- c("serial", "snow")
    rows <- list()
    for (storage in c("dense", "sparse")) {
        dims <- dimensions[[storage]]
        for (overlap in overlaps) {
            for (backend in backends) {
                rows[[length(rows) + 1L]] <- benchmark_scenario_row(
                    preset = preset,
                    storage = storage,
                    overlap = overlap,
                    backend = backend,
                    backend_workers = workers,
                    n_features = dims[["features"]],
                    n_samples = dims[["samples"]],
                    n_sets = dims[["sets"]],
                    set_size = dims[["size"]],
                    matrix_seed = if (storage == "dense") 1741L else 1742L,
                    set_seed = if (overlap == "low") 1743L else 1744L,
                    density = if (storage == "sparse") 0.03 else NA_real_,
                    zero_fraction = if (storage == "dense") 0.35 else NA_real_,
                    cell_missing_fraction = if (storage == "dense") 0.01 else NA_real_,
                    stored_missing_fraction = if (storage == "sparse") 0.05 else NA_real_
                )
            }
        }
    }
    do.call(rbind, rows)
}

benchmark_feature_ids <- function(n_features) {
    sprintf("feature_%06d", seq_len(n_features))
}

benchmark_sample_ids <- function(n_samples) {
    sprintf("sample_%04d", seq_len(n_samples))
}

benchmark_gene_sets <- function(
    features,
    n_sets,
    set_size,
    overlap,
    seed
) {
    n_sets <- benchmark_assert_count(n_sets, "`n_sets`")
    set_size <- benchmark_assert_count(set_size, "`set_size`", minimum = 2L)
    if (set_size > length(features)) {
        stop("`set_size` exceeds the available feature count.", call. = FALSE)
    }
    if (!overlap %in% c("low", "high")) {
        stop("`overlap` must be low or high.", call. = FALSE)
    }

    set.seed(seed)
    if (overlap == "low") {
        member_count <- as.double(n_sets) * set_size
        if (member_count > length(features)) {
            stop(
                "Low-overlap fixtures require n_sets * set_size <= n_features.",
                call. = FALSE
            )
        }
        selected <- sample(features, member_count, replace = FALSE)
        sets <- split(selected, rep.int(seq_len(n_sets), rep.int(set_size, n_sets)))
    } else {
        common_size <- max(1L, min(set_size - 1L, floor(0.75 * set_size)))
        common <- sample(features, common_size, replace = FALSE)
        remaining <- setdiff(features, common)
        sets <- lapply(seq_len(n_sets), function(index) {
            sample(c(
                common,
                sample(remaining, set_size - common_size, replace = FALSE)
            ))
        })
    }
    names(sets) <- sprintf("set_%04d", seq_len(n_sets))
    sets
}

benchmark_dense_matrix <- function(
    n_features,
    n_samples,
    zero_fraction,
    missing_fraction,
    seed
) {
    zero_fraction <- benchmark_assert_fraction(zero_fraction, "`zero_fraction`")
    missing_fraction <- benchmark_assert_fraction(
        missing_fraction,
        "`missing_fraction`"
    )
    total <- as.double(n_features) * n_samples
    set.seed(seed)
    mat <- matrix(stats::rexp(total), nrow = n_features, ncol = n_samples)

    zero_count <- as.integer(round(total * zero_fraction))
    if (zero_count > 0L) {
        mat[sample.int(total, zero_count, replace = FALSE)] <- 0
    }
    missing_count <- as.integer(round(total * missing_fraction))
    if (missing_count > 0L) {
        positions <- sample.int(total, missing_count, replace = FALSE)
        mat[positions] <- rep(c(NA_real_, NaN), length.out = missing_count)
    }
    dimnames(mat) <- list(
        benchmark_feature_ids(n_features),
        benchmark_sample_ids(n_samples)
    )
    mat
}

benchmark_sparse_matrix <- function(
    n_features,
    n_samples,
    density,
    stored_missing_fraction,
    seed
) {
    density <- benchmark_assert_fraction(density, "`density`")
    stored_missing_fraction <- benchmark_assert_fraction(
        stored_missing_fraction,
        "`stored_missing_fraction`"
    )
    total <- as.double(n_features) * n_samples
    stored_count <- as.integer(round(total * density))
    set.seed(seed)
    positions <- if (stored_count == 0L) {
        numeric()
    } else {
        sample.int(total, stored_count, replace = FALSE)
    }
    rows <- as.integer((positions - 1) %% n_features + 1)
    columns <- as.integer((positions - 1) %/% n_features + 1)
    values <- stats::rexp(stored_count)

    missing_count <- as.integer(round(stored_count * stored_missing_fraction))
    if (missing_count > 0L) {
        missing <- sample.int(stored_count, missing_count, replace = FALSE)
        values[missing] <- rep(c(NA_real_, NaN), length.out = missing_count)
    }

    Matrix::sparseMatrix(
        i = rows,
        j = columns,
        x = values,
        dims = c(n_features, n_samples),
        dimnames = list(
            benchmark_feature_ids(n_features),
            benchmark_sample_ids(n_samples)
        ),
        repr = "C"
    )
}

benchmark_fixture <- function(scenario) {
    required <- c(
        "storage", "overlap", "n_features", "n_samples", "n_sets",
        "set_size", "matrix_seed", "set_seed"
    )
    missing <- setdiff(required, names(scenario))
    if (length(missing) > 0L) {
        stop("Scenario lacks fields: ", paste(missing, collapse = ", "), call. = FALSE)
    }

    mat <- switch(
        scenario$storage,
        dense = benchmark_dense_matrix(
            scenario$n_features,
            scenario$n_samples,
            scenario$zero_fraction,
            scenario$cell_missing_fraction,
            scenario$matrix_seed
        ),
        sparse = benchmark_sparse_matrix(
            scenario$n_features,
            scenario$n_samples,
            scenario$density,
            scenario$stored_missing_fraction,
            scenario$matrix_seed
        ),
        stop("Unknown scenario storage.", call. = FALSE)
    )
    gene_sets <- benchmark_gene_sets(
        rownames(mat),
        scenario$n_sets,
        scenario$set_size,
        scenario$overlap,
        scenario$set_seed
    )
    list(mat = mat, gene_sets = gene_sets)
}

benchmark_fixture_stats <- function(fixture) {
    mat <- fixture$mat
    sparse <- inherits(mat, "sparseMatrix")
    stored_values <- if (sparse) methods::slot(mat, "x") else mat
    total_cells <- as.double(nrow(mat)) * ncol(mat)
    stored_entries <- if (sparse) length(stored_values) else total_cells
    missing_cells <- sum(is.na(stored_values))
    zero_cells <- if (sparse) {
        total_cells - stored_entries + sum(stored_values == 0, na.rm = TRUE)
    } else {
        sum(stored_values == 0, na.rm = TRUE)
    }

    members <- unlist(fixture$gene_sets, use.names = FALSE)
    unique_members <- length(unique(members))
    data.frame(
        input_bytes = as.numeric(object.size(mat)),
        logical_dense_bytes = total_cells * 8,
        total_cells = total_cells,
        stored_entries = as.double(stored_entries),
        missing_cells = as.double(missing_cells),
        zero_cells = as.double(zero_cells),
        catalog_memberships = as.double(length(members)),
        catalog_unique_features = as.double(unique_members),
        catalog_reuse = length(members) / unique_members,
        check.names = FALSE
    )
}
