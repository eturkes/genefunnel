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
    below_center <- observed < center

    # Algebraically identical to the normative subtraction, but all terms are
    # non-negative. This avoids catastrophic cancellation for valid inputs.
    (
        (n_observed - 1L - sum(below_center)) * total +
            n_observed * sum(observed[below_center])
    ) / (n_observed - 1L)
}
