# Assisted-by: OpenAI Codex.

.GF_BIG_BASE_BITS <- 26L
.GF_BIG_BASE <- 2^.GF_BIG_BASE_BITS

.gf_big_trim <- function(value) {
    if (length(value) == 0L) {
        return(0)
    }
    last <- length(value)
    while (last > 1L && value[[last]] == 0) {
        last <- last - 1L
    }
    if (last == length(value)) {
        return(value)
    }
    unname(value[seq_len(last)])
}

.gf_big_is_zero <- function(value) {
    value <- .gf_big_trim(value)
    length(value) == 1L && value[[1L]] == 0
}

.gf_big_from_integer <- function(value) {
    if (length(value) != 1L || !is.finite(value) || value < 0 ||
        value != floor(value)) {
        stop("Internal exact integer is invalid.", call. = FALSE)
    }
    if (value == 0) {
        return(0)
    }
    limbs <- numeric()
    while (value > 0) {
        quotient <- floor(value / .GF_BIG_BASE)
        limbs <- c(limbs, value - quotient * .GF_BIG_BASE)
        value <- quotient
    }
    .gf_big_trim(limbs)
}

.gf_big_compare <- function(left, right) {
    left <- .gf_big_trim(left)
    right <- .gf_big_trim(right)
    if (length(left) != length(right)) {
        return(if (length(left) < length(right)) -1L else 1L)
    }
    for (index in seq.int(length(left), 1L)) {
        if (left[[index]] != right[[index]]) {
            return(if (left[[index]] < right[[index]]) -1L else 1L)
        }
    }
    0L
}

.gf_big_add <- function(left, right) {
    size <- max(length(left), length(right))
    result <- numeric(size + 1L)
    carry <- 0
    for (index in seq_len(size)) {
        left_value <- if (index <= length(left)) left[[index]] else 0
        right_value <- if (index <= length(right)) right[[index]] else 0
        total <- left_value + right_value + carry
        carry <- floor(total / .GF_BIG_BASE)
        result[[index]] <- total - carry * .GF_BIG_BASE
    }
    result[[size + 1L]] <- carry
    .gf_big_trim(result)
}

.gf_big_subtract <- function(left, right) {
    if (.gf_big_compare(left, right) < 0L) {
        stop("Internal exact subtraction would be negative.", call. = FALSE)
    }
    result <- numeric(length(left))
    borrow <- 0
    for (index in seq_along(left)) {
        right_value <- if (index <= length(right)) right[[index]] else 0
        difference <- left[[index]] - right_value - borrow
        if (difference < 0) {
            difference <- difference + .GF_BIG_BASE
            borrow <- 1
        } else {
            borrow <- 0
        }
        result[[index]] <- difference
    }
    if (borrow != 0) {
        stop("Internal exact subtraction borrow escaped.", call. = FALSE)
    }
    .gf_big_trim(result)
}

.gf_big_multiply <- function(left, right) {
    if (.gf_big_is_zero(left) || .gf_big_is_zero(right)) {
        return(0)
    }
    result <- numeric(length(left) + length(right) + 1L)
    for (left_index in seq_along(left)) {
        carry <- 0
        for (right_index in seq_along(right)) {
            index <- left_index + right_index - 1L
            total <- result[[index]] +
                left[[left_index]] * right[[right_index]] + carry
            carry <- floor(total / .GF_BIG_BASE)
            result[[index]] <- total - carry * .GF_BIG_BASE
        }
        index <- left_index + length(right)
        while (carry > 0) {
            total <- result[[index]] + carry
            carry <- floor(total / .GF_BIG_BASE)
            result[[index]] <- total - carry * .GF_BIG_BASE
            index <- index + 1L
        }
    }
    .gf_big_trim(result)
}

.gf_big_multiply_small <- function(value, multiplier) {
    .gf_big_multiply(value, .gf_big_from_integer(multiplier))
}

.gf_big_shift_left <- function(value, bits) {
    if (bits < 0 || bits != floor(bits)) {
        stop("Internal exact shift is invalid.", call. = FALSE)
    }
    if (.gf_big_is_zero(value) || bits == 0) {
        return(value)
    }
    words <- floor(bits / .GF_BIG_BASE_BITS)
    remainder <- bits - words * .GF_BIG_BASE_BITS
    result <- c(rep.int(0, words), value, 0)
    carry <- 0
    factor <- 2^remainder
    for (index in seq_along(value) + words) {
        total <- result[[index]] * factor + carry
        carry <- floor(total / .GF_BIG_BASE)
        result[[index]] <- total - carry * .GF_BIG_BASE
    }
    result[[length(value) + words + 1L]] <- carry
    .gf_big_trim(result)
}

.gf_big_sum <- function(values) {
    Reduce(.gf_big_add, values, init = 0)
}

.gf_scale_pow2 <- function(value, shift) {
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

.gf_dyadic_value <- function(value) {
    if (!is.finite(value) || value < 0) {
        stop("Internal dyadic input is invalid.", call. = FALSE)
    }
    if (value == 0) {
        return(list(magnitude = 0, exponent = 0L))
    }
    exponent <- floor(log2(value))
    normalized <- .gf_scale_pow2(value, -exponent)
    while (normalized < 1) {
        normalized <- normalized * 2
        exponent <- exponent - 1
    }
    while (normalized >= 2) {
        normalized <- normalized / 2
        exponent <- exponent + 1
    }
    significand <- normalized * 2^52
    rounded <- round(significand)
    if (significand != rounded || rounded < 1 || rounded > 2^53) {
        stop("Internal dyadic decomposition failed.", call. = FALSE)
    }
    exponent <- exponent - 52
    while (rounded %% 2 == 0) {
        rounded <- rounded / 2
        exponent <- exponent + 1
    }
    list(
        magnitude = .gf_big_from_integer(rounded),
        exponent = as.integer(exponent)
    )
}

.gf_dyadic_vector <- function(values) {
    pieces <- lapply(values, .gf_dyadic_value)
    positive <- !vapply(
        pieces,
        function(value) .gf_big_is_zero(value$magnitude),
        logical(1L)
    )
    common_exponent <- if (any(positive)) {
        min(vapply(pieces[positive], `[[`, integer(1L), "exponent"))
    } else {
        0L
    }
    magnitudes <- lapply(pieces, function(value) {
        if (.gf_big_is_zero(value$magnitude)) {
            return(0)
        }
        .gf_big_shift_left(
            value$magnitude,
            value$exponent - common_exponent
        )
    })
    list(magnitudes = magnitudes, exponent = common_exponent)
}

.gf_exact_score_numerator <- function(magnitudes) {
    size <- length(magnitudes)
    if (size < 2L) {
        stop("Internal exact score support is invalid.", call. = FALSE)
    }
    total <- .gf_big_sum(magnitudes)
    below <- vapply(magnitudes, function(value) {
        .gf_big_compare(
            .gf_big_multiply_small(value, size),
            total
        ) < 0L
    }, logical(1L))
    below_sum <- .gf_big_sum(magnitudes[below])
    numerator <- .gf_big_add(
        .gf_big_multiply_small(total, size - 1L - sum(below)),
        .gf_big_multiply_small(below_sum, size)
    )
    list(numerator = numerator, total = total)
}

.gf_signed_difference <- function(left, right) {
    comparison <- .gf_big_compare(left, right)
    if (comparison == 0L) {
        return(list(sign = 0L, magnitude = 0))
    }
    if (comparison > 0L) {
        return(list(
            sign = 1L,
            magnitude = .gf_big_subtract(left, right)
        ))
    }
    list(sign = -1L, magnitude = .gf_big_subtract(right, left))
}

.gf_exact_deltas <- function(values) {
    size <- length(values)
    if (size < 3L) {
        stop("Internal exact delta support is invalid.", call. = FALSE)
    }
    dyadic <- .gf_dyadic_vector(values)
    full <- .gf_exact_score_numerator(dyadic$magnitudes)
    denominator <- .gf_big_multiply(
        .gf_big_from_integer(size - 1L),
        .gf_big_from_integer(size - 2L)
    )
    deltas <- lapply(seq_len(size), function(index) {
        deleted <- .gf_exact_score_numerator(dyadic$magnitudes[-index])
        .gf_signed_difference(
            .gf_big_multiply_small(full$numerator, size - 2L),
            .gf_big_multiply_small(deleted$numerator, size - 1L)
        )
    })
    list(
        deltas = deltas,
        denominator = denominator,
        total = full$total,
        exponent = dyadic$exponent
    )
}

.gf_big_prefix_sums <- function(values) {
    prefix <- vector("list", length(values) + 1L)
    prefix[[1L]] <- 0
    for (index in seq_along(values)) {
        prefix[[index + 1L]] <- .gf_big_add(prefix[[index]], values[[index]])
    }
    prefix
}

.gf_sorted_below_count <- function(values, total, multiplier) {
    lower <- 0L
    upper <- length(values)
    while (lower < upper) {
        middle <- (lower + upper + 1L) %/% 2L
        below <- .gf_big_compare(
            .gf_big_multiply_small(values[[middle]], multiplier),
            total
        ) < 0L
        if (below) {
            lower <- middle
        } else {
            upper <- middle - 1L
        }
    }
    lower
}

.gf_deleted_score_numerator <- function(
    index,
    magnitudes,
    positions,
    sorted,
    prefix,
    total
) {
    remaining_size <- length(magnitudes) - 1L
    remaining_total <- .gf_big_subtract(total, magnitudes[[index]])
    boundary <- .gf_sorted_below_count(
        sorted,
        remaining_total,
        remaining_size
    )
    included <- positions[[index]] <= boundary
    below_sum <- prefix[[boundary + 1L]]
    if (included) {
        below_sum <- .gf_big_subtract(below_sum, magnitudes[[index]])
    }
    below_count <- boundary - as.integer(included)
    .gf_big_add(
        .gf_big_multiply_small(
            remaining_total,
            remaining_size - 1L - below_count
        ),
        .gf_big_multiply_small(below_sum, remaining_size)
    )
}

.gf_exact_deltas_sorted <- function(values) {
    size <- length(values)
    if (size < 3L) {
        stop("Internal exact delta support is invalid.", call. = FALSE)
    }
    dyadic <- .gf_dyadic_vector(values)
    full <- .gf_exact_score_numerator(dyadic$magnitudes)
    ordering <- order(values, method = "radix")
    sorted <- dyadic$magnitudes[ordering]
    prefix <- .gf_big_prefix_sums(sorted)
    if (.gf_big_compare(prefix[[size + 1L]], full$total) != 0L) {
        stop("Internal exact prefix sum is invalid.", call. = FALSE)
    }
    positions <- integer(size)
    positions[ordering] <- seq_len(size)
    deltas <- lapply(seq_len(size), function(index) {
        deleted <- .gf_deleted_score_numerator(
            index, dyadic$magnitudes, positions, sorted, prefix, full$total
        )
        .gf_signed_difference(
            .gf_big_multiply_small(full$numerator, size - 2L),
            .gf_big_multiply_small(deleted, size - 1L)
        )
    })
    list(
        deltas = deltas,
        denominator = .gf_big_multiply(
            .gf_big_from_integer(size - 1L),
            .gf_big_from_integer(size - 2L)
        ),
        total = full$total,
        exponent = dyadic$exponent
    )
}

.gf_big_merge_block <- function(indices, values, left, middle, right) {
    left_indices <- indices[seq.int(left, middle)]
    right_indices <- indices[seq.int(middle + 1L, right)]
    merged <- integer(length(left_indices) + length(right_indices))
    left_position <- 1L
    right_position <- 1L
    for (position in seq_along(merged)) {
        use_left <- right_position > length(right_indices) ||
            (left_position <= length(left_indices) && .gf_big_compare(
                values[[left_indices[[left_position]]]],
                values[[right_indices[[right_position]]]]
            ) <= 0L)
        if (use_left) {
            merged[[position]] <- left_indices[[left_position]]
            left_position <- left_position + 1L
        } else {
            merged[[position]] <- right_indices[[right_position]]
            right_position <- right_position + 1L
        }
    }
    merged
}

.gf_big_order <- function(values) {
    size <- length(values)
    indices <- seq_len(size)
    width <- 1L
    while (width < size) {
        next_indices <- integer(size)
        for (left in seq.int(1L, size, by = 2L * width)) {
            middle <- min(left + width - 1L, size)
            right <- min(left + 2L * width - 1L, size)
            next_indices[seq.int(left, right)] <- if (middle == right) {
                indices[seq.int(left, right)]
            } else {
                .gf_big_merge_block(indices, values, left, middle, right)
            }
        }
        indices <- next_indices
        width <- 2L * width
    }
    indices
}

.gf_big_scaled_head <- function(value) {
    value <- .gf_big_trim(value)
    if (.gf_big_is_zero(value)) {
        return(list(mantissa = 0, exponent = 0L))
    }
    first <- max(1L, length(value) - 2L)
    head <- 0
    for (index in seq.int(length(value), first)) {
        head <- head * .GF_BIG_BASE + value[[index]]
    }
    exponent <- floor(log2(head))
    mantissa <- .gf_scale_pow2(head, -exponent)
    while (mantissa >= 2) {
        mantissa <- mantissa / 2
        exponent <- exponent + 1
    }
    list(
        mantissa = mantissa,
        exponent = as.integer(
            exponent + .GF_BIG_BASE_BITS * (first - 1L)
        )
    )
}

.gf_ratio_candidate <- function(numerator, denominator, power_exponent) {
    top <- .gf_big_scaled_head(numerator)
    bottom <- .gf_big_scaled_head(denominator)
    mantissa <- top$mantissa / bottom$mantissa
    exponent <- top$exponent - bottom$exponent + power_exponent
    while (mantissa < 1) {
        mantissa <- mantissa * 2
        exponent <- exponent - 1L
    }
    while (mantissa >= 2) {
        mantissa <- mantissa / 2
        exponent <- exponent + 1L
    }
    .gf_scale_pow2(mantissa, exponent)
}

.gf_ratio_error_ok <- function(
    numerator,
    denominator,
    power_exponent,
    rounded,
    effective_size
) {
    rounded_dyadic <- .gf_dyadic_value(abs(rounded))
    common_exponent <- min(power_exponent, rounded_dyadic$exponent)
    exact_numerator <- .gf_big_shift_left(
        numerator,
        power_exponent - common_exponent
    )
    rounded_numerator <- .gf_big_shift_left(
        .gf_big_multiply(rounded_dyadic$magnitude, denominator),
        rounded_dyadic$exponent - common_exponent
    )
    error <- .gf_signed_difference(
        exact_numerator,
        rounded_numerator
    )$magnitude
    .gf_big_compare(
        .gf_big_shift_left(error, 49L),
        .gf_big_multiply_small(exact_numerator, effective_size)
    ) <= 0L
}

.gf_validate_ratio_inputs <- function(
    denominator,
    power_exponent,
    sign,
    effective_size
) {
    valid_integer <- function(value) {
        length(value) == 1L && !is.na(value) && is.finite(value) &&
            value == floor(value)
    }
    if (.gf_big_is_zero(denominator) ||
        !valid_integer(power_exponent) ||
        !valid_integer(sign) || !sign %in% -1:1 ||
        !valid_integer(effective_size) || effective_size < 1) {
        stop("Internal exact ratio is invalid.", call. = FALSE)
    }
}

.gf_big_ratio_double <- function(
    numerator,
    denominator,
    power_exponent,
    sign,
    effective_size
) {
    .gf_validate_ratio_inputs(
        denominator,
        power_exponent,
        sign,
        effective_size
    )
    if (.gf_big_is_zero(numerator) || sign == 0L) {
        return(0)
    }
    candidate <- .gf_ratio_candidate(
        numerator,
        denominator,
        power_exponent
    )
    if (!is.finite(candidate) || candidate == 0 ||
        !.gf_ratio_error_ok(
            numerator,
            denominator,
            power_exponent,
            candidate,
            effective_size
        )) {
        return(NA_real_)
    }
    sign * candidate
}
