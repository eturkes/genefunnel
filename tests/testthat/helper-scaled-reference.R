# Assisted-by: OpenAI Codex.

# Test-only double-double significand + explicit binary exponent. Arithmetic
# stays independent of the native scorer and never relies on platform long
# double width.

scaled_ref_two_sum <- function(left, right) {
    total <- left + right
    right_virtual <- total - left
    error <- (left - (total - right_virtual)) + (right - right_virtual)
    c(hi = total, lo = error)
}

scaled_ref_two_product <- function(left, right) {
    product <- left * right
    splitter <- 134217729
    left_split <- splitter * left
    right_split <- splitter * right
    left_hi <- left_split - (left_split - left)
    right_hi <- right_split - (right_split - right)
    left_lo <- left - left_hi
    right_lo <- right - right_hi
    error <- ((left_hi * right_hi - product) +
        left_hi * right_lo + left_lo * right_hi) +
        left_lo * right_lo
    c(hi = product, lo = error)
}

scaled_ref_pow2 <- function(value, shift) {
    stopifnot(length(value) == 1L, length(shift) == 1L, shift == floor(shift))
    if (value == 0 || shift == 0) {
        return(value)
    }
    while (shift > 0) {
        step <- min(shift, 1023)
        value <- value * 2^step
        shift <- shift - step
    }
    while (shift < 0) {
        step <- max(shift, -1022)
        value <- value * 2^step
        shift <- shift - step
    }
    value
}

scaled_ref_normalize <- function(hi, lo = 0, exponent = 0) {
    combined <- scaled_ref_two_sum(hi, lo)
    hi <- unname(combined[["hi"]])
    lo <- unname(combined[["lo"]])
    if (hi == 0) {
        if (lo == 0) {
            return(list(hi = 0, lo = 0, exponent = 0))
        }
        hi <- lo
        lo <- 0
    }

    repeat {
        magnitude <- abs(hi)
        if (magnitude >= 0.5 && magnitude < 1) {
            break
        }
        adjustment <- floor(log2(magnitude)) + 1
        hi <- scaled_ref_pow2(hi, -adjustment)
        lo <- scaled_ref_pow2(lo, -adjustment)
        exponent <- exponent + adjustment
        combined <- scaled_ref_two_sum(hi, lo)
        hi <- unname(combined[["hi"]])
        lo <- unname(combined[["lo"]])
    }
    list(hi = hi, lo = lo, exponent = exponent)
}

scaled_ref_from_double <- function(value) {
    stopifnot(length(value) == 1L, is.finite(value))
    scaled_ref_normalize(value)
}

scaled_ref_negate <- function(value) {
    list(
        hi = -value$hi,
        lo = -value$lo,
        exponent = value$exponent
    )
}

scaled_ref_add <- function(left, right) {
    if (left$hi == 0) {
        return(right)
    }
    if (right$hi == 0) {
        return(left)
    }
    exponent <- max(left$exponent, right$exponent)
    left_shift <- left$exponent - exponent
    right_shift <- right$exponent - exponent
    left_hi <- scaled_ref_pow2(left$hi, left_shift)
    left_lo <- scaled_ref_pow2(left$lo, left_shift)
    right_hi <- scaled_ref_pow2(right$hi, right_shift)
    right_lo <- scaled_ref_pow2(right$lo, right_shift)

    leading <- scaled_ref_two_sum(left_hi, right_hi)
    trailing <- scaled_ref_two_sum(left_lo, right_lo)
    middle <- scaled_ref_two_sum(
        unname(leading[["lo"]]),
        unname(trailing[["hi"]])
    )
    result <- scaled_ref_two_sum(
        unname(leading[["hi"]]),
        unname(middle[["hi"]])
    )
    lo <- unname(result[["lo"]]) + unname(middle[["lo"]]) +
        unname(trailing[["lo"]])
    scaled_ref_normalize(unname(result[["hi"]]), lo, exponent)
}

scaled_ref_compare <- function(left, right) {
    if (left$hi < 0 || right$hi < 0) {
        stop("scaled reference comparison requires non-negative values")
    }
    if (left$hi == 0 && right$hi == 0) {
        return(0L)
    }
    if (left$hi == 0) {
        return(-1L)
    }
    if (right$hi == 0) {
        return(1L)
    }
    if (left$exponent != right$exponent) {
        return(if (left$exponent < right$exponent) -1L else 1L)
    }
    if (left$hi != right$hi) {
        return(if (left$hi < right$hi) -1L else 1L)
    }
    if (left$lo != right$lo) {
        return(if (left$lo < right$lo) -1L else 1L)
    }
    0L
}

scaled_ref_subtract <- function(left, right) {
    if (scaled_ref_compare(left, right) < 0L) {
        stop("scaled reference subtraction would be negative")
    }
    scaled_ref_add(left, scaled_ref_negate(right))
}

scaled_ref_multiply <- function(left, right) {
    if (left$hi == 0 || right$hi == 0) {
        return(scaled_ref_from_double(0))
    }
    leading <- scaled_ref_two_product(left$hi, right$hi)
    cross <- left$hi * right$lo + left$lo * right$hi
    trailing <- unname(leading[["lo"]]) + cross + left$lo * right$lo
    combined <- scaled_ref_two_sum(unname(leading[["hi"]]), trailing)
    scaled_ref_normalize(
        unname(combined[["hi"]]),
        unname(combined[["lo"]]),
        left$exponent + right$exponent
    )
}

scaled_ref_multiply_double <- function(value, multiplier) {
    stopifnot(length(multiplier) == 1L, is.finite(multiplier), multiplier >= 0)
    scaled_ref_multiply(value, scaled_ref_from_double(multiplier))
}

scaled_ref_divide <- function(numerator, denominator) {
    stopifnot(numerator$hi >= 0, denominator$hi > 0)
    if (numerator$hi == 0) {
        return(scaled_ref_from_double(0))
    }
    quotient <- scaled_ref_normalize(
        numerator$hi / denominator$hi,
        exponent = numerator$exponent - denominator$exponent
    )
    for (iteration in seq_len(3L)) {
        product <- scaled_ref_multiply(denominator, quotient)
        comparison <- scaled_ref_compare(numerator, product)
        if (comparison == 0L) {
            break
        }
        residual <- if (comparison > 0L) {
            scaled_ref_subtract(numerator, product)
        } else {
            scaled_ref_subtract(product, numerator)
        }
        correction <- scaled_ref_normalize(
            residual$hi / denominator$hi,
            exponent = residual$exponent - denominator$exponent
        )
        quotient <- if (comparison > 0L) {
            scaled_ref_add(quotient, correction)
        } else {
            scaled_ref_subtract(quotient, correction)
        }
    }
    quotient
}

scaled_ref_divide_double <- function(value, divisor) {
    stopifnot(length(divisor) == 1L, is.finite(divisor), divisor > 0)
    scaled_ref_divide(value, scaled_ref_from_double(divisor))
}

scaled_ref_abs_difference <- function(left, right) {
    if (scaled_ref_compare(left, right) >= 0L) {
        scaled_ref_subtract(left, right)
    } else {
        scaled_ref_subtract(right, left)
    }
}

scaled_ref_sum <- function(values) {
    Reduce(scaled_ref_add, values, init = scaled_ref_from_double(0))
}

scaled_ref_pair <- function(value) {
    normalized <- scaled_ref_normalize(value$hi, value$lo, value$exponent)
    if (normalized$hi == 0) {
        return(c(mantissa = 0, exponent = 0))
    }
    mantissa <- normalized$hi + normalized$lo
    if (abs(mantissa) >= 1) {
        mantissa <- mantissa / 2
        normalized$exponent <- normalized$exponent + 1
    }
    c(mantissa = mantissa, exponent = normalized$exponent)
}

scaled_ref_as_double <- function(value) {
    pair <- scaled_ref_pair(value)
    ordinary <- scaled_ref_pow2(pair[["mantissa"]], pair[["exponent"]])
    if (pair[["mantissa"]] != 0 && (!is.finite(ordinary) || ordinary == 0)) {
        return(NA_real_)
    }
    ordinary
}

scaled_ref_relative_error <- function(observed, expected) {
    difference <- scaled_ref_abs_difference(observed, expected)
    scaled_ref_as_double(scaled_ref_divide(difference, expected))
}

reference_scaled_components <- function(values) {
    stopifnot(
        is.numeric(values),
        length(values) >= 2L,
        !any(is.infinite(values)),
        !any(values < 0, na.rm = TRUE)
    )
    observed <- values[!is.na(values)]
    effective_size <- length(observed)
    observed_fraction <- effective_size / length(values)
    scaled_values <- lapply(observed, scaled_ref_from_double)
    observed_sum <- scaled_ref_sum(scaled_values)

    if (effective_size < 2L) {
        return(list(
            semantic = "too_few_observed",
            effective_size = effective_size,
            observed_fraction = observed_fraction,
            observed_sum = observed_sum,
            score = NULL,
            penalty = NULL,
            balance = NULL,
            conditioning = "not_applicable"
        ))
    }
    if (observed_sum$hi == 0) {
        zero <- scaled_ref_from_double(0)
        return(list(
            semantic = "zero_total",
            effective_size = effective_size,
            observed_fraction = observed_fraction,
            observed_sum = observed_sum,
            score = zero,
            penalty = zero,
            balance = NULL,
            conditioning = "not_applicable"
        ))
    }

    center <- scaled_ref_divide_double(observed_sum, effective_size)
    deviations <- lapply(
        scaled_values,
        scaled_ref_abs_difference,
        right = center
    )
    penalty <- scaled_ref_divide_double(
        scaled_ref_multiply_double(
            scaled_ref_sum(deviations),
            effective_size
        ),
        2 * (effective_size - 1L)
    )

    below <- vapply(
        scaled_values,
        scaled_ref_compare,
        integer(1L),
        right = center
    ) < 0L
    below_sum <- scaled_ref_sum(scaled_values[below])
    score <- scaled_ref_divide_double(
        scaled_ref_add(
            scaled_ref_multiply_double(
                observed_sum,
                effective_size - 1L - sum(below)
            ),
            scaled_ref_multiply_double(below_sum, effective_size)
        ),
        effective_size - 1L
    )
    balance <- scaled_ref_divide(score, observed_sum)
    ordinary <- vapply(
        list(observed_sum, penalty, balance, score),
        scaled_ref_as_double,
        numeric(1L)
    )
    kappa <- if (score$hi == 0) {
        NULL
    } else {
        scaled_ref_divide(
            scaled_ref_add(observed_sum, penalty),
            score
        )
    }
    ordinary_kappa <- if (is.null(kappa)) Inf else scaled_ref_as_double(kappa)
    unit_roundoff <- .Machine$double.eps / 2
    allowance <- 64 * (effective_size + 1) * unit_roundoff
    safe <- all(is.finite(ordinary)) &&
        ordinary[[4L]] >= .Machine$double.xmin &&
        is.finite(ordinary_kappa) &&
        allowance < 1 && allowance * ordinary_kappa <= 2^-20

    list(
        semantic = "scoreable",
        effective_size = effective_size,
        observed_fraction = observed_fraction,
        observed_sum = observed_sum,
        score = score,
        penalty = penalty,
        balance = balance,
        kappa = kappa,
        allowance = allowance,
        conditioning = if (safe) "safe" else "ill_conditioned"
    )
}
