reference_score <- function(values) {
    stopifnot(
        is.numeric(values),
        !any(is.infinite(values)),
        !any(values < 0, na.rm = TRUE)
    )

    observed <- values[!is.na(values)]
    n_observed <- length(observed)
    if (n_observed < 2L) {
        return(NA_real_)
    }

    total <- sum(observed)
    center <- total / n_observed
    below_center <- if (center == 0 && total > 0) {
        # `total / n_observed` can underflow for subnormal values. In that
        # case only exact zeros are below the positive mathematical mean.
        observed == 0
    } else {
        observed < center
    }

    # Algebraically identical to the normative subtraction, but all terms are
    # non-negative. This avoids catastrophic cancellation for valid inputs.
    (
        (n_observed - 1L - sum(below_center)) * total +
            n_observed * sum(observed[below_center])
    ) / (n_observed - 1L)
}

reference_memberships <- function(gene_sets, features) {
    unique_members <- lapply(gene_sets, function(members) {
        members[!duplicated(members)]
    })
    indices <- lapply(unique_members, function(members) {
        matched <- match(members, features, nomatch = 0L)
        unname(matched[matched > 0L])
    })
    indices[lengths(indices) >= 2L]
}

reference_cells <- function(mat, gene_sets, FUN) {
    indices <- reference_memberships(gene_sets, rownames(mat))
    result <- matrix(
        NA_real_,
        nrow = length(indices),
        ncol = ncol(mat),
        dimnames = list(names(indices), colnames(mat))
    )

    for (set_index in seq_along(indices)) {
        for (sample_index in seq_len(ncol(mat))) {
            result[set_index, sample_index] <- FUN(
                mat[indices[[set_index]], sample_index]
            )
        }
    }
    result
}

reference_scores <- function(mat, gene_sets) {
    reference_cells(mat, gene_sets, reference_score)
}
