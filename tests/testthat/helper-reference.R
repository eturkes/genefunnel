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
    total -
        n_observed / (2 * (n_observed - 1L)) *
            sum(abs(observed - center))
}
