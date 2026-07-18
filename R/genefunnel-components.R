# Assisted-by: OpenAI Codex.

#' Calculate GeneFunnel score components
#'
#' Calculate the authoritative GeneFunnel score together with aligned
#' magnitude, balance, penalty, and observed-support diagnostics. This is an
#' additive audit API: it does not change [genefunnel()] or reconstruct its
#' score from rounded components.
#'
#' @inheritParams genefunnel
#'
#' @details
#' For each retained gene-set/sample cell, let `effective_size` be the number
#' of matched members remaining after sample-specific `NA`/`NaN` omission and
#' let `observed_sum` be their sum. `observed_fraction` divides effective size
#' by the set's globally matched size. It is distinct from declared-set
#' coverage returned by [gene_set_coverage()].
#'
#' When at least two values remain and their sum is positive, `balance` is
#' Bulla's evenness index for the within-set shares and the score factorizes as
#' `score = observed_sum * balance`. `penalty = observed_sum - score` in exact
#' arithmetic. All-zero cells have score and penalty zero but undefined
#' balance. Cells with fewer than two observed values have factual support and
#' sum diagnostics but undefined score, penalty, and balance.
#'
#' Finite inputs can make a mathematical component overflow or underflow an R
#' double even when the authoritative score is representable. Such a component
#' has ordinary value `NA`, status `"scaled"`, and a canonical sidecar pair
#' `(mantissa, exponent)` representing `mantissa * 2^exponent`, with
#' `0.5 <= mantissa < 1`. Ordinary and undefined components have `NA` sidecar
#' pairs. Status `"unavailable"` marks an undefined or uncertifiable
#' diagnostic; it never changes the score.
#'
#' `conditioning` is `"safe"` only in the protocol's conservative binary64
#' region where both component reconstructions are numerically checkable.
#' `"ill_conditioned"` directs callers to the scaled representation rather
#' than rounded subtraction. `"not_applicable"` covers zero totals and cells
#' with fewer than two observed members.
#'
#' @return A fixed-order named list of aligned base matrices:
#'
#' - `score`, `observed_sum`, `penalty`, and `balance`: double matrices;
#' - `effective_size`: an integer matrix;
#' - `observed_fraction`: a double matrix;
#' - `status`: character matrices `semantic`, `observed_sum`, `penalty`,
#'   `balance`, and `conditioning`; and
#' - `scaled`: entries `observed_sum`, `penalty`, and `balance`, each containing
#'   double `mantissa` and integer `exponent` matrices.
#'
#' Every matrix uses the retained set/sample order and names of [genefunnel()].
#' Semantic status is `"scoreable"`, `"too_few_observed"`, or `"zero_total"`.
#' Component availability is `"ordinary"`, `"scaled"`, or `"unavailable"`.
#'
#' @section Interpretation:
#' Balance describes arithmetic evenness within one observed cell. It is not
#' co-expression, correlation, coordinated regulation, pathway truth, causal
#' attribution, uncertainty, or preprocessing invariance. It depends on the
#' effective member support, feature-specific abundance and dynamic range,
#' measurement noise, and meaningful-zero assumption.
#'
#' @references
#' Bulla L (1994). "An index of evenness and its associated diversity
#' measure." *Oikos*, 70(1), 167-171. \doi{10.2307/3545713}.
#'
#' @seealso [genefunnel()], [gene_set_coverage()],
#'   `system.file("COMPONENTS_SPEC.md", package = "genefunnel")`
#' @export
#'
#' @examples
#' mat <- rbind(
#'     A = c(equal = 2, unequal = 6, missing = 1),
#'     B = c(equal = 2, unequal = 2, missing = 2),
#'     C = c(equal = NA, unequal = NA, missing = NA)
#' )
#' components <- genefunnel_components(
#'     mat,
#'     list(pair = c("A", "B"), triple = c("A", "B", "C")),
#'     BPPARAM = BiocParallel::SerialParam()
#' )
#' components$score
#' components$observed_sum
#' components$balance
#' components$effective_size
genefunnel_components <- function(
    mat,
    gene_sets,
    BPPARAM = BiocParallel::bpparam()
) {
    .validate_score_matrix(mat)
    prepared <- .prepare_gene_sets(gene_sets, rownames(mat))
    .validate_bpparam(BPPARAM)

    retained <- prepared$coverage$scoreable
    if (any(!retained)) {
        .warn_unscoreable_sets(prepared$coverage$gene_set[!retained])
    }

    retained_names <- prepared$coverage$gene_set[retained]
    result_dimnames <- list(retained_names, colnames(mat))
    if (!any(retained)) {
        return(.new_component_result(
            n_gene_sets = 0L,
            n_columns = ncol(mat),
            dimnames = result_dimnames
        ))
    }

    indices <- prepared$indices[retained]
    storage <- .matrix_storage(mat)
    ranges <- .column_chunk_ranges(
        ncol(mat),
        BiocParallel::bpnworkers(BPPARAM)
    )
    chunks <- BiocParallel::bpiterate(
        .matrix_chunk_iterator(mat, ranges),
        .component_matrix_task,
        gene_indices = indices,
        storage = storage,
        BPPARAM = BPPARAM
    )

    .assemble_component_chunks(
        chunks,
        n_gene_sets = length(indices),
        n_columns = ncol(mat),
        dimnames = result_dimnames
    )
}

.component_matrix_chunk <- function(mat, gene_indices, storage) {
    switch(
        storage,
        dense = calculateComponentsDense(
            .as_dense_numeric_chunk(mat),
            gene_indices
        ),
        sparse = calculateComponentsSparse(
            .as_sparse_numeric_chunk(mat),
            gene_indices
        ),
        stop("Internal matrix storage kind is invalid.", call. = FALSE)
    )
}

.component_matrix_task <- function(task, gene_indices, storage, ...) {
    tryCatch(
        list(
            id = task$id,
            first = task$first,
            last = task$last,
            components = .component_matrix_chunk(
                task$mat,
                gene_indices,
                storage
            )
        ),
        error = function(condition) {
            stop(
                sprintf(
                    "GeneFunnel component calculation failed in chunk %d (%s): %s",
                    task$id,
                    .format_chunk_context(task),
                    conditionMessage(condition)
                ),
                call. = FALSE
            )
        }
    )
}

.new_component_result <- function(n_gene_sets, n_columns, dimnames = NULL) {
    double_matrix <- function() {
        matrix(
            NA_real_,
            nrow = n_gene_sets,
            ncol = n_columns,
            dimnames = dimnames
        )
    }
    integer_matrix <- function() {
        matrix(
            NA_integer_,
            nrow = n_gene_sets,
            ncol = n_columns,
            dimnames = dimnames
        )
    }
    character_matrix <- function() {
        matrix(
            NA_character_,
            nrow = n_gene_sets,
            ncol = n_columns,
            dimnames = dimnames
        )
    }
    scaled_component <- function() {
        list(mantissa = double_matrix(), exponent = integer_matrix())
    }

    list(
        score = double_matrix(),
        observed_sum = double_matrix(),
        penalty = double_matrix(),
        balance = double_matrix(),
        effective_size = integer_matrix(),
        observed_fraction = double_matrix(),
        status = list(
            semantic = character_matrix(),
            observed_sum = character_matrix(),
            penalty = character_matrix(),
            balance = character_matrix(),
            conditioning = character_matrix()
        ),
        scaled = list(
            observed_sum = scaled_component(),
            penalty = scaled_component(),
            balance = scaled_component()
        )
    )
}

.validate_component_matrix <- function(value, type, dimensions, label) {
    if (!is.matrix(value) || is.object(value) ||
        !identical(typeof(value), type) ||
        !identical(dim(value), dimensions)) {
        stop(
            "Internal component calculation returned invalid ",
            label,
            ".",
            call. = FALSE
        )
    }
}

.component_result_fields <- function() {
    c(
        "score",
        "observed_sum",
        "penalty",
        "balance",
        "effective_size",
        "observed_fraction",
        "status",
        "scaled"
    )
}

.component_status_fields <- function() {
    c(
        "semantic",
        "observed_sum",
        "penalty",
        "balance",
        "conditioning"
    )
}

.component_scaled_fields <- function() {
    c("observed_sum", "penalty", "balance")
}

.validate_component_schema <- function(value) {
    if (!is.list(value) ||
        !identical(names(value), .component_result_fields()) ||
        !is.list(value$status) ||
        !identical(names(value$status), .component_status_fields()) ||
        !is.list(value$scaled) ||
        !identical(names(value$scaled), .component_scaled_fields())) {
        stop(
            "Internal component calculation returned an invalid schema.",
            call. = FALSE
        )
    }
}

.validate_component_ordinary_matrices <- function(value, dimensions) {
    for (field in c(
        "score",
        "observed_sum",
        "penalty",
        "balance",
        "observed_fraction"
    )) {
        .validate_component_matrix(
            value[[field]],
            "double",
            dimensions,
            paste0("`", field, "` matrix")
        )
    }
    .validate_component_matrix(
        value$effective_size,
        "integer",
        dimensions,
        "`effective_size` matrix"
    )
    for (field in .component_status_fields()) {
        .validate_component_matrix(
            value$status[[field]],
            "character",
            dimensions,
            paste0("status `", field, "` matrix")
        )
    }
}

.validate_component_scaled_matrices <- function(value, dimensions) {
    for (field in .component_scaled_fields()) {
        component <- value$scaled[[field]]
        if (!is.list(component) ||
            !identical(names(component), c("mantissa", "exponent"))) {
            stop(
                "Internal component calculation returned an invalid scaled `",
                field,
                "` schema.",
                call. = FALSE
            )
        }
        .validate_component_matrix(
            component$mantissa,
            "double",
            dimensions,
            paste0("scaled `", field, "` mantissa matrix")
        )
        .validate_component_matrix(
            component$exponent,
            "integer",
            dimensions,
            paste0("scaled `", field, "` exponent matrix")
        )
    }
}

.validate_component_chunk <- function(value, dimensions) {
    .validate_component_schema(value)
    .validate_component_ordinary_matrices(value, dimensions)
    .validate_component_scaled_matrices(value, dimensions)
    invisible(value)
}

.order_component_chunks <- function(chunks, n_columns) {
    if (length(chunks) == 0L) {
        stop(
            "Internal parallel component calculation returned no chunks.",
            call. = FALSE
        )
    }
    chunk_ids <- vapply(chunks, `[[`, integer(1L), "id")
    expected_ids <- seq_along(chunks)
    if (!identical(unname(sort(chunk_ids)), expected_ids)) {
        stop(
            "Internal parallel component calculation returned invalid ",
            "chunk identifiers.",
            call. = FALSE
        )
    }
    chunks <- chunks[order(chunk_ids)]

    first <- vapply(chunks, `[[`, integer(1L), "first")
    last <- vapply(chunks, `[[`, integer(1L), "last")
    valid_ranges <- first[[1L]] == 1L &&
        last[[length(last)]] == n_columns &&
        all(last >= first) &&
        (length(first) == 1L ||
            identical(first[-1L], last[-length(last)] + 1L))
    if (!valid_ranges) {
        stop(
            "Internal parallel component calculation returned invalid ",
            "column ranges.",
            call. = FALSE
        )
    }
    chunks
}

.copy_component_chunk <- function(result, chunk, n_gene_sets) {
    expected_dimensions <- c(
        as.integer(n_gene_sets),
        chunk$last - chunk$first + 1L
    )
    .validate_component_chunk(chunk$components, expected_dimensions)
    columns <- seq.int(chunk$first, chunk$last)
    matrix_fields <- c(
        "score",
        "observed_sum",
        "penalty",
        "balance",
        "effective_size",
        "observed_fraction"
    )
    for (field in matrix_fields) {
        result[[field]][, columns] <- chunk$components[[field]]
    }
    for (field in .component_status_fields()) {
        result$status[[field]][, columns] <-
            chunk$components$status[[field]]
    }
    for (field in .component_scaled_fields()) {
        result$scaled[[field]]$mantissa[, columns] <-
            chunk$components$scaled[[field]]$mantissa
        result$scaled[[field]]$exponent[, columns] <-
            chunk$components$scaled[[field]]$exponent
    }
    result
}

.assemble_component_chunks <- function(
    chunks,
    n_gene_sets,
    n_columns,
    dimnames
) {
    chunks <- .order_component_chunks(chunks, n_columns)
    result <- .new_component_result(n_gene_sets, n_columns, dimnames)
    for (chunk in chunks) {
        result <- .copy_component_chunk(result, chunk, n_gene_sets)
    }
    result
}
