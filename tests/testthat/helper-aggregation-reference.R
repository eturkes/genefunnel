# Assisted-by: OpenAI Codex.

reference_aggregation_gap <- function(units, weights) {
    stopifnot(
        is.matrix(units),
        is.numeric(units),
        ncol(units) >= 2L,
        is.numeric(weights),
        length(weights) == nrow(units),
        !anyNA(weights),
        all(is.finite(weights)),
        all(weights >= 0),
        any(weights > 0)
    )
    weight_tolerance <- 32 * .Machine$double.eps * max(1, length(weights))
    stopifnot(abs(sum(weights) - 1) <= weight_tolerance)

    active <- weights > 0
    active_units <- units[active, , drop = FALSE]
    active_weights <- weights[active]
    stopifnot(
        !anyNA(active_units),
        all(is.finite(active_units)),
        all(active_units >= 0)
    )

    centers <- rowMeans(active_units)
    centered <- sweep(active_units, 1L, centers, "-")
    weighted <- function(values) {
        colSums(sweep(values, 1L, active_weights, "*"))
    }
    pooled <- weighted(active_units)
    member_count <- ncol(active_units)
    coefficient <- member_count / (2 * (member_count - 1L))
    coordinate_gap <- weighted(abs(centered)) - abs(weighted(centered))
    positive_mass <- weighted(pmax(centered, 0))
    negative_mass <- weighted(pmax(-centered, 0))
    formula_gap <- coefficient * sum(coordinate_gap)
    cancellation_gap <- 2 * coefficient * sum(pmin(
        positive_mass,
        negative_mass
    ))
    unit_scores <- apply(active_units, 1L, reference_score)
    weighted_unit_score <- sum(active_weights * unit_scores)
    aggregate_score <- reference_score(pooled)
    score_gap <- aggregate_score - weighted_unit_score

    list(
        aggregate = pooled,
        aggregate_score = aggregate_score,
        unit_scores = unit_scores,
        weighted_unit_score = weighted_unit_score,
        score_gap = score_gap,
        formula_gap = formula_gap,
        cancellation_gap = cancellation_gap,
        normalized_gap = if (aggregate_score > 0) {
            formula_gap / aggregate_score
        } else {
            NA_real_
        },
        coordinate_gap = coordinate_gap,
        positive_mass = positive_mass,
        negative_mass = negative_mass,
        no_opposing_deviations = all(
            positive_mass == 0 | negative_mass == 0
        ),
        active = active,
        weights = active_weights
    )
}
