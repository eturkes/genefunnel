# Assisted-by: OpenAI Codex.

.GF_CATALOGUE_SCHEMA_VERSION <- 1L
.GF_FORMULA_VERSION <- "GF-1"
.GF_CATALOGUE_ENCODING <- "GFCAT-1"
.GF_CATALOGUE_MAGIC <- charToRaw("GFCAT")
.GF_CATALOGUE_DOMAIN <- charToRaw("GFCAN")
.GF_FEATURE_DOMAIN <- charToRaw("GFFEA")
.GF_SHA256_BYTES <- 32L

.catalogue_wire_stop <- function(detail) {
    stop("Malformed GFCAT-1 payload: ", detail, ".", call. = FALSE)
}

.catalogue_concat_raw <- function(chunks) {
    if (length(chunks) == 0L) {
        return(raw())
    }
    unname(do.call(c, chunks))
}

.catalogue_encode_u32 <- function(values) {
    if (!is.numeric(values) || anyNA(values) ||
        any(!is.finite(values)) || any(values < 0) ||
        any(values != floor(values)) ||
        any(values > .Machine$integer.max)) {
        stop(
            "Internal catalogue integer is outside the wire range.",
            call. = FALSE
        )
    }
    if (length(values) == 0L) {
        return(raw())
    }
    values <- as.double(values)
    bytes <- rbind(
        floor(values / 16777216) %% 256,
        floor(values / 65536) %% 256,
        floor(values / 256) %% 256,
        values %% 256
    )
    as.raw(as.vector(bytes))
}

.catalogue_encode_u32_vector <- function(values) {
    c(.catalogue_encode_u32(length(values)), .catalogue_encode_u32(values))
}

.catalogue_identifier_tag <- function(value) {
    encoding <- Encoding(value)
    match(encoding, c("unknown", "UTF-8", "latin1", "bytes"), nomatch = 0L) - 1L
}

.catalogue_encode_identifier <- function(value) {
    if (!is.character(value) || length(value) != 1L || is.object(value) ||
        !is.null(attributes(value)) || is.na(value) || !nzchar(value)) {
        stop("Internal catalogue identifier is invalid.", call. = FALSE)
    }
    tag <- .catalogue_identifier_tag(value)
    if (tag < 0L) {
        stop(
            "Internal catalogue identifier encoding is invalid.",
            call. = FALSE
        )
    }
    bytes <- charToRaw(value)
    c(as.raw(tag), .catalogue_encode_u32(length(bytes)), bytes)
}

.catalogue_encode_identifiers <- function(values) {
    chunks <- lapply(unname(as.list(values)), .catalogue_encode_identifier)
    c(.catalogue_encode_u32(length(values)), .catalogue_concat_raw(chunks))
}

.catalogue_feature_body <- function(state) {
    c(
        .GF_FEATURE_DOMAIN,
        .catalogue_encode_u32(state$schema_version),
        .catalogue_encode_identifiers(state$features)
    )
}

.catalogue_catalogue_body <- function(state) {
    set_chunks <- vector("list", length(state$members))
    for (index in seq_along(state$members)) {
        set_chunks[[index]] <- c(
            .catalogue_encode_identifier(names(state$members)[[index]]),
            .catalogue_encode_u32(state$duplicate_member_count[[index]]),
            .catalogue_encode_identifiers(state$members[[index]])
        )
    }
    c(
        .GF_CATALOGUE_DOMAIN,
        .catalogue_encode_u32(state$schema_version),
        .catalogue_encode_identifier(state$formula_version),
        .catalogue_encode_u32(length(state$members)),
        .catalogue_concat_raw(set_chunks)
    )
}

.catalogue_content_body <- function(state) {
    set_chunks <- vector("list", length(state$members))
    for (index in seq_along(state$members)) {
        set_chunks[[index]] <- c(
            .catalogue_encode_identifier(names(state$members)[[index]]),
            .catalogue_encode_u32(state$duplicate_member_count[[index]]),
            .catalogue_encode_identifiers(state$members[[index]]),
            .catalogue_encode_u32_vector(state$member_index[[index]])
        )
    }
    c(
        .GF_CATALOGUE_MAGIC,
        .catalogue_encode_u32(state$schema_version),
        .catalogue_encode_identifier(state$formula_version),
        .catalogue_encode_identifiers(state$features),
        .catalogue_encode_u32(length(state$members)),
        .catalogue_concat_raw(set_chunks),
        .catalogue_encode_u32_vector(state$row_adjacency$offsets),
        .catalogue_encode_u32_vector(state$row_adjacency$set_indices)
    )
}

.catalogue_sha256 <- function(bytes) {
    unname(tools::sha256sum(bytes = bytes))
}

.catalogue_hex_to_raw <- function(value) {
    if (!is.character(value) || length(value) != 1L ||
        !grepl("^[0-9a-f]{64}$", value)) {
        stop("Internal SHA-256 text is invalid.", call. = FALSE)
    }
    starts <- seq.int(1L, 63L, by = 2L)
    as.raw(strtoi(substring(value, starts, starts + 1L), base = 16L))
}

.catalogue_raw_to_hex <- function(value) {
    paste(sprintf("%02x", as.integer(value)), collapse = "")
}

.catalogue_append_digest <- function(body) {
    c(body, .catalogue_hex_to_raw(.catalogue_sha256(body)))
}

.catalogue_fingerprints <- function(state) {
    c(
        algorithm = "SHA-256",
        encoding = .GF_CATALOGUE_ENCODING,
        catalogue = .catalogue_sha256(.catalogue_catalogue_body(state)),
        features = .catalogue_sha256(.catalogue_feature_body(state)),
        content = .catalogue_sha256(.catalogue_content_body(state))
    )
}

.catalogue_reader <- function(body) {
    reader <- new.env(parent = emptyenv())
    reader$body <- body
    reader$position <- 1
    reader$limit <- as.double(length(body))
    reader
}

.catalogue_reader_remaining <- function(reader) {
    reader$limit - reader$position + 1L
}

.catalogue_read_raw <- function(reader, count, label) {
    if (!is.numeric(count) || length(count) != 1L || is.na(count) ||
        !is.finite(count) || count < 0 || count != floor(count) ||
        count > .Machine$integer.max) {
        .catalogue_wire_stop(paste0(
            label,
            " length is outside the integer range"
        ))
    }
    count <- as.integer(count)
    if (count > .catalogue_reader_remaining(reader)) {
        .catalogue_wire_stop(paste0(label, " is truncated"))
    }
    if (count == 0L) {
        return(raw())
    }
    first <- reader$position
    last <- first + count - 1L
    reader$position <- last + 1
    reader$body[seq.int(first, last)]
}

.catalogue_decode_u32 <- function(bytes, label) {
    value <- sum(as.integer(bytes) * c(16777216, 65536, 256, 1))
    if (value > .Machine$integer.max) {
        .catalogue_wire_stop(paste0(
            label,
            " is outside the supported integer range"
        ))
    }
    as.integer(value)
}

.catalogue_read_u32 <- function(reader, label) {
    .catalogue_decode_u32(.catalogue_read_raw(reader, 4L, label), label)
}

.catalogue_read_u32_vector <- function(reader, label) {
    count <- .catalogue_read_u32(reader, paste0(label, " count"))
    remaining <- .catalogue_reader_remaining(reader)
    if (count > floor(remaining / 4)) {
        .catalogue_wire_stop(paste0(label, " length exceeds remaining bytes"))
    }
    if (count == 0L) {
        return(integer())
    }
    bytes <- as.integer(.catalogue_read_raw(reader, 4 * count, label))
    dim(bytes) <- c(4L, count)
    values <- colSums(bytes * c(16777216, 65536, 256, 1))
    if (any(values > .Machine$integer.max)) {
        .catalogue_wire_stop(paste0(
            label,
            " contains an integer outside the supported range"
        ))
    }
    as.integer(values)
}

.catalogue_read_identifier <- function(reader, label) {
    tag <- as.integer(.catalogue_read_raw(
        reader,
        1L,
        paste0(label, " encoding")
    ))
    encodings <- c("unknown", "UTF-8", "latin1", "bytes")
    if (tag > 3L) {
        .catalogue_wire_stop(paste0(label, " has an invalid encoding tag"))
    }
    byte_count <- .catalogue_read_u32(reader, paste0(label, " byte length"))
    if (byte_count == 0L) {
        .catalogue_wire_stop(paste0(label, " is empty"))
    }
    bytes <- .catalogue_read_raw(reader, byte_count, label)
    value <- tryCatch(
        rawToChar(bytes),
        error = function(condition) {
            .catalogue_wire_stop(paste0(
                label,
                " contains invalid string bytes"
            ))
        }
    )
    expected_encoding <- encodings[[tag + 1L]]
    Encoding(value) <- expected_encoding
    if (!identical(Encoding(value), expected_encoding)) {
        .catalogue_wire_stop(paste0(label, " has a non-canonical encoding tag"))
    }
    value
}

.catalogue_read_identifiers <- function(reader, label) {
    count <- .catalogue_read_u32(reader, paste0(label, " count"))
    if (count > floor(.catalogue_reader_remaining(reader) / 6)) {
        .catalogue_wire_stop(paste0(label, " count exceeds remaining bytes"))
    }
    values <- character(count)
    for (index in seq_len(count)) {
        values[[index]] <- .catalogue_read_identifier(
            reader,
            sprintf("%s %d", label, index)
        )
    }
    values
}

.catalogue_parse_sets <- function(reader, set_count) {
    set_names <- character(set_count)
    members <- vector("list", set_count)
    duplicate_member_count <- integer(set_count)
    member_index <- vector("list", set_count)
    for (index in seq_len(set_count)) {
        set_names[[index]] <- .catalogue_read_identifier(
            reader,
            sprintf("set name %d", index)
        )
        duplicate_member_count[[index]] <- .catalogue_read_u32(
            reader,
            sprintf("duplicate count %d", index)
        )
        members[[index]] <- .catalogue_read_identifiers(
            reader,
            sprintf("set %d member", index)
        )
        member_index[[index]] <- .catalogue_read_u32_vector(
            reader,
            sprintf("set %d member index", index)
        )
        if (length(member_index[[index]]) != length(members[[index]])) {
            .catalogue_wire_stop(sprintf(
                "set %d member/index lengths disagree",
                index
            ))
        }
    }
    names(members) <- set_names
    names(member_index) <- set_names
    list(
        members = members,
        duplicate_member_count = duplicate_member_count,
        member_index = member_index
    )
}

.catalogue_parse_body <- function(body) {
    reader <- .catalogue_reader(body)
    magic <- .catalogue_read_raw(reader, length(.GF_CATALOGUE_MAGIC), "magic")
    if (!identical(magic, .GF_CATALOGUE_MAGIC)) {
        .catalogue_wire_stop("magic is invalid")
    }
    schema_version <- .catalogue_read_u32(reader, "schema version")
    if (!identical(schema_version, .GF_CATALOGUE_SCHEMA_VERSION)) {
        .catalogue_wire_stop("schema version is unsupported")
    }
    formula_version <- .catalogue_read_identifier(reader, "formula version")
    if (!identical(formula_version, .GF_FORMULA_VERSION)) {
        .catalogue_wire_stop("formula version is unsupported")
    }
    features <- .catalogue_read_identifiers(reader, "feature")
    set_count <- .catalogue_read_u32(reader, "set count")
    maximum_sets <- floor(.catalogue_reader_remaining(reader) / 18)
    if (set_count == 0L || set_count > maximum_sets) {
        .catalogue_wire_stop("set count is invalid for the remaining bytes")
    }
    sets <- .catalogue_parse_sets(reader, set_count)
    offsets <- .catalogue_read_u32_vector(reader, "adjacency offset")
    set_indices <- .catalogue_read_u32_vector(reader, "adjacency set index")
    if (.catalogue_reader_remaining(reader) != 0L) {
        .catalogue_wire_stop("body has trailing bytes")
    }
    list(
        schema_version = schema_version,
        formula_version = formula_version,
        features = features,
        members = sets$members,
        duplicate_member_count = sets$duplicate_member_count,
        member_index = sets$member_index,
        row_adjacency = list(offsets = offsets, set_indices = set_indices)
    )
}

.serialize_compiled_catalogue <- function(catalogue) {
    state <- .catalogue_state(catalogue)
    .catalogue_append_digest(.catalogue_content_body(state))
}

.unserialize_compiled_catalogue <- function(bytes) {
    if (!is.raw(bytes) || !is.null(attributes(bytes))) {
        stop("`bytes` must be an unclassed raw vector.", call. = FALSE)
    }
    if (length(bytes) <= .GF_SHA256_BYTES) {
        .catalogue_wire_stop("payload is truncated")
    }
    if (length(bytes) > .Machine$integer.max) {
        .catalogue_wire_stop(
            "payload length is outside the supported integer range"
        )
    }
    body_length <- length(bytes) - .GF_SHA256_BYTES
    body <- bytes[seq_len(body_length)]
    supplied_digest <- .catalogue_raw_to_hex(
        bytes[seq.int(body_length + 1L, length(bytes))]
    )
    observed_digest <- .catalogue_sha256(body)
    if (!identical(supplied_digest, observed_digest)) {
        .catalogue_wire_stop("content fingerprint does not match")
    }

    parsed <- .catalogue_parse_body(body)
    restored <- tryCatch(
        .build_compiled_catalogue(
            parsed$members,
            parsed$duplicate_member_count,
            parsed$features
        ),
        error = function(condition) {
            .catalogue_wire_stop(paste0(
                "decoded primary state is invalid: ",
                conditionMessage(condition)
            ))
        }
    )
    state <- unclass(restored)
    if (!identical(parsed$member_index, state$member_index) ||
        !identical(parsed$row_adjacency, state$row_adjacency) ||
        !identical(observed_digest, state$fingerprints[["content"]])) {
        .catalogue_wire_stop(
            "decoded mappings or adjacency disagree with identifiers"
        )
    }
    restored
}
