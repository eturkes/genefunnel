# Assisted-by: OpenAI Codex.

sensitivity_score_reference <- function(values) {
    observed <- values[!is.na(values)]
    n_observed <- length(observed)
    if (n_observed < 2L) {
        return(NA_real_)
    }

    total <- sum(observed)
    center <- total / n_observed
    total - n_observed / (2 * (n_observed - 1L)) *
        sum(abs(observed - center))
}

sensitivity_members_reference <- function(members, features, values) {
    canonical <- members[!duplicated(members)]
    matched <- canonical[match(canonical, features, nomatch = 0L) > 0L]
    matched_values <- unname(values[match(matched, features)])
    observed <- !is.na(matched_values)

    list(
        members = matched[observed],
        values = matched_values[observed]
    )
}

sensitivity_deltas_reference <- function(values, members = NULL) {
    if (is.null(members)) {
        members <- names(values)
    }
    observed <- !is.na(values)
    values <- unname(values[observed])
    members <- members[observed]
    if (length(values) < 3L) {
        return(setNames(numeric(), character()))
    }

    full <- sensitivity_score_reference(values)
    deleted <- vapply(seq_along(values), function(index) {
        sensitivity_score_reference(values[-index])
    }, numeric(1))
    setNames(full - deleted, members)
}

sensitivity_summary_reference <- function(values, members = NULL) {
    if (is.null(members)) {
        members <- names(values)
    }
    effective_size <- sum(!is.na(values))
    deltas <- sensitivity_deltas_reference(values, members)
    if (effective_size < 3L) {
        return(list(
            effective_size = effective_size,
            largest_member = NA_character_,
            largest_absolute_delta = NA_real_,
            largest_delta = NA_real_,
            largest_delta_over_sum = NA_real_,
            median_absolute_delta = NA_real_
        ))
    }

    largest <- which.max(abs(deltas))
    total <- sum(values, na.rm = TRUE)
    list(
        effective_size = effective_size,
        largest_member = names(deltas)[largest],
        largest_absolute_delta = abs(unname(deltas[largest])),
        largest_delta = unname(deltas[largest]),
        largest_delta_over_sum = if (total > 0) {
            unname(deltas[largest]) / total
        } else {
            NA_real_
        },
        median_absolute_delta = stats::median(abs(deltas))
    )
}
