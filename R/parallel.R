.column_chunk_ranges <- function(n_columns, n_workers) {
    valid_column_count <- is.numeric(n_columns) &&
        length(n_columns) == 1L &&
        !is.na(n_columns) &&
        is.finite(n_columns) &&
        n_columns >= 1 &&
        n_columns == floor(n_columns)
    if (!valid_column_count) {
        stop("Internal matrix column count is invalid.", call. = FALSE)
    }

    valid_worker_count <- is.numeric(n_workers) &&
        length(n_workers) == 1L &&
        !is.na(n_workers) &&
        is.finite(n_workers) &&
        n_workers >= 1 &&
        n_workers == floor(n_workers)
    if (!valid_worker_count) {
        stop(
            "`BPPARAM` must provide at least one worker.",
            call. = FALSE
        )
    }

    # Two tasks per parallel worker allow limited load balancing while keeping
    # the task count bounded by the worker count.
    target_tasks <- if (n_workers == 1) 1 else 2 * as.double(n_workers)
    task_count <- as.integer(min(n_columns, target_tasks))
    boundaries <- floor(
        seq.int(0, task_count) * as.double(n_columns) / task_count
    )

    data.frame(
        id = seq_len(task_count),
        first = as.integer(boundaries[-length(boundaries)] + 1),
        last = as.integer(boundaries[-1L]),
        row.names = seq_len(task_count)
    )
}

.matrix_chunk_iterator <- function(mat, ranges) {
    next_id <- 0L

    function() {
        if (next_id >= nrow(ranges)) {
            return(NULL)
        }
        next_id <<- next_id + 1L
        first <- ranges$first[[next_id]]
        last <- ranges$last[[next_id]]
        chunk <- if (first == 1L && last == ncol(mat)) {
            mat
        } else {
            mat[, seq.int(first, last), drop = FALSE]
        }

        list(
            id = ranges$id[[next_id]],
            first = first,
            last = last,
            mat = chunk
        )
    }
}

.score_matrix_task <- function(task, gene_indices, storage, ...) {
    tryCatch(
        list(
            id = task$id,
            first = task$first,
            last = task$last,
            scores = .score_matrix_chunk(
                task$mat,
                gene_indices,
                storage
            )
        ),
        error = function(condition) {
            stop(
                sprintf(
                    "GeneFunnel scoring failed in chunk %d (%s): %s",
                    task$id,
                    .format_chunk_context(task),
                    conditionMessage(condition)
                ),
                call. = FALSE
            )
        }
    )
}

.format_chunk_context <- function(task) {
    columns <- if (task$first == task$last) {
        sprintf("matrix column %d", task$first)
    } else {
        sprintf("matrix columns %d-%d", task$first, task$last)
    }
    sample_names <- colnames(task$mat)
    if (is.null(sample_names)) {
        return(columns)
    }

    format_name <- function(name) {
        if (is.na(name)) "<NA>" else encodeString(name, quote = '"')
    }
    samples <- if (length(sample_names) == 1L) {
        sprintf("sample %s", format_name(sample_names[[1L]]))
    } else {
        sprintf(
            "samples %s through %s",
            format_name(sample_names[[1L]]),
            format_name(sample_names[[length(sample_names)]])
        )
    }
    paste(columns, samples, sep = ", ")
}

.assemble_score_chunks <- function(chunks, n_gene_sets, n_columns) {
    if (length(chunks) == 0L) {
        stop(
            "Internal parallel scoring returned no chunks.",
            call. = FALSE
        )
    }
    chunk_ids <- vapply(chunks, `[[`, integer(1), "id")
    expected_ids <- seq_along(chunks)
    if (!identical(unname(sort(chunk_ids)), expected_ids)) {
        stop(
            "Internal parallel scoring returned invalid chunk identifiers.",
            call. = FALSE
        )
    }
    chunks <- chunks[order(chunk_ids)]

    first <- vapply(chunks, `[[`, integer(1), "first")
    last <- vapply(chunks, `[[`, integer(1), "last")
    valid_ranges <- first[[1L]] == 1L &&
        last[[length(last)]] == n_columns &&
        all(last >= first) &&
        (length(first) == 1L ||
            identical(first[-1L], last[-length(last)] + 1L))
    if (!valid_ranges) {
        stop(
            "Internal parallel scoring returned invalid column ranges.",
            call. = FALSE
        )
    }

    scores <- matrix(NA_real_, nrow = n_gene_sets, ncol = n_columns)
    for (chunk in chunks) {
        expected_dimensions <- c(
            n_gene_sets,
            chunk$last - chunk$first + 1L
        )
        if (!identical(dim(chunk$scores), expected_dimensions)) {
            stop(
                "Internal parallel scoring returned invalid score dimensions.",
                call. = FALSE
            )
        }
        scores[, seq.int(chunk$first, chunk$last)] <- chunk$scores
    }
    scores
}
