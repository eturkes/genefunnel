# Assisted-by: OpenAI Codex.

SENSITIVITY_CONTROLLED_VERSION <- "E-C-1.0.0"
SENSITIVITY_CONTROLLED_MD5 <- "1c119004ee749b4242f1b73ca6fd8c4c"

sensitivity_controlled_read_protocol <- function(root = ".") {
    path <- file.path(root, "benchmark", "sensitivity-controlled-protocol.tsv")
    if (!identical(unname(tools::md5sum(path)), SENSITIVITY_CONTROLLED_MD5)) {
        stop("Sensitivity controlled protocol identity is invalid.", call. = FALSE)
    }
    protocol <- utils::read.delim(
        path,
        stringsAsFactors = FALSE,
        check.names = FALSE,
        colClasses = "character",
        na.strings = character()
    )
    required <- c("protocol_version", "section", "key", "value")
    identity <- paste(protocol$section, protocol$key, sep = "::")
    valid <- identical(names(protocol), required) && nrow(protocol) == 45L &&
        !anyNA(protocol) && all(nzchar(unlist(protocol))) &&
        all(protocol$protocol_version == SENSITIVITY_CONTROLLED_VERSION) &&
        !anyDuplicated(identity)
    if (!valid) {
        stop("Sensitivity controlled protocol schema is invalid.", call. = FALSE)
    }
    protocol
}

sensitivity_controlled_protocol_value <- function(protocol, section, key) {
    selected <- protocol$section == section & protocol$key == key
    if (sum(selected) != 1L) {
        stop("Sensitivity controlled protocol key is invalid.", call. = FALSE)
    }
    protocol$value[[which(selected)]]
}

sensitivity_controlled_validate_protocol <- function(root = ".") {
    protocol <- sensitivity_controlled_read_protocol(root)
    parent <- sensitivity_controlled_protocol_value(
        protocol, "identity", "parent_protocol"
    )
    parent_md5 <- sensitivity_controlled_protocol_value(
        protocol, "identity", "parent_registry_md5"
    )
    observed <- unname(tools::md5sum(file.path(
        root, "benchmark", "sensitivity-protocol.tsv"
    )))
    if (!identical(parent, SENSITIVITY_PROTOCOL_VERSION) ||
        !identical(parent_md5, observed)) {
        stop("Sensitivity controlled parent identity is invalid.", call. = FALSE)
    }
    list(
        protocol_version = SENSITIVITY_CONTROLLED_VERSION,
        protocol_rows = nrow(protocol),
        protocol_md5 = SENSITIVITY_CONTROLLED_MD5,
        parent_protocol = parent,
        parent_registry_md5 = parent_md5
    )
}
