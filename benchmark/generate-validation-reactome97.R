# Assisted-by: OpenAI Codex.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L ||
    !startsWith(args[[1L]], "--input=") ||
    !startsWith(args[[2L]], "--output=")) {
    stop(
        "Usage: Rscript benchmark/generate-validation-reactome97.R ",
        "--input=REACTOME97_DIR --output=TSV",
        call. = FALSE
    )
}
input_dir <- normalizePath(sub("^--input=", "", args[[1L]]), mustWork = TRUE)
output <- sub("^--output=", "", args[[2L]])
dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
output <- file.path(normalizePath(dirname(output), mustWork = TRUE), basename(output))
if (!identical(Sys.setlocale("LC_COLLATE", "C"), "C")) {
    stop("C collation is required.", call. = FALSE)
}

runner <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(runner) != 1L) stop("Cannot locate generator.", call. = FALSE)
benchmark_dir <- dirname(normalizePath(sub("^--file=", "", runner)))
parent_protocol_version <- "F-2.0.0"
parent_protocol_sha256 <-
    "00a3a9f5e1cd6f709fc3ecdbad111f938b80c849c1c6071d134db41e7c63361c"
protocol_path <- file.path(benchmark_dir, "validation-protocol.tsv")
observed_protocol_sha256 <- unname(tools::sha256sum(protocol_path))
if (!identical(observed_protocol_sha256, parent_protocol_sha256)) {
    stop("Validation protocol identity is invalid.", call. = FALSE)
}
protocol <- utils::read.delim(
    protocol_path, stringsAsFactors = FALSE, check.names = FALSE,
    colClasses = "character", na.strings = character()
)
protocol_fields <- c("protocol_version", "section", "key", "value")
protocol_identity <- paste(protocol$section, protocol$key, sep = "::")
protocol_valid <- identical(names(protocol), protocol_fields) &&
    nrow(protocol) == 184L && !anyNA(protocol) &&
    all(nzchar(unlist(protocol))) &&
    all(protocol$protocol_version == parent_protocol_version) &&
    !anyDuplicated(protocol_identity) &&
    all(grepl("^[a-z][a-z0-9_]*$", protocol$section)) &&
    all(grepl("^[A-Za-z][A-Za-z0-9_]*$", protocol$key))
if (!protocol_valid) {
    stop("Validation protocol rows are invalid.", call. = FALSE)
}
protocol_value <- function(section, key) {
    selected <- protocol$section == section & protocol$key == key
    if (sum(selected) != 1L) {
        stop("Validation protocol key is invalid: ", section, "/", key,
            call. = FALSE)
    }
    protocol$value[[which(selected)]]
}

files <- c(
    rna = "Ensembl2Reactome_All_Levels.txt",
    protein = "UniProt2Reactome_All_Levels.txt",
    pathway = "ReactomePathways.txt"
)
paths <- stats::setNames(file.path(input_dir, files), names(files))
expected <- vapply(names(files), function(kind) {
    protocol_value("catalogue", paste0(kind, "_sha256"))
}, character(1L))
observed <- unname(tools::sha256sum(paths))
if (!identical(observed, unname(expected))) {
    stop("Reactome 97 input identity is invalid.", call. = FALSE)
}

human_rows <- function(path, columns) {
    connection <- file(path, open = "rt")
    on.exit(close(connection))
    result <- list()
    index <- 0L
    repeat {
        lines <- readLines(connection, n = 100000L, warn = FALSE)
        if (!length(lines)) break
        lines <- lines[endsWith(lines, "\tHomo sapiens")]
        if (!length(lines)) next
        fields <- strsplit(lines, "\t", fixed = TRUE)
        if (any(lengths(fields) != columns)) {
            stop("Reactome 97 human row shape is invalid.", call. = FALSE)
        }
        index <- index + 1L
        result[[index]] <- fields
    }
    unlist(result, recursive = FALSE, use.names = FALSE)
}

pathway_fields <- human_rows(paths[["pathway"]], 3L)
pathway_ids <- vapply(pathway_fields, `[[`, character(1L), 1L)
if (any(!grepl("^R-HSA-[0-9]+$", pathway_ids)) || anyDuplicated(pathway_ids)) {
    stop("Reactome 97 human pathway registry is invalid.", call. = FALSE)
}

mapping_pairs <- function(path, member_kind) {
    fields <- human_rows(path, 6L)
    member <- vapply(fields, `[[`, character(1L), 1L)
    target <- vapply(fields, `[[`, character(1L), 2L)
    if (member_kind == "rna") {
        member_valid <- grepl("^ENSG[0-9]+$", member)
    } else {
        member <- sub("-[1-9][0-9]*$", "", member)
        member_valid <- grepl(
            "^([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})$",
            member
        )
    }
    keep <- member_valid & grepl("^R-HSA-[0-9]+$", target) &
        target %in% pathway_ids
    unique(paste(target[keep], member[keep], sep = "\t"))
}

eligible <- function(pairs) {
    target <- sub("\t.*$", "", pairs)
    sizes <- table(target)
    keep <- sizes >= 10L & sizes <= 200L
    data.frame(
        target_id = names(sizes)[keep],
        declared_size = as.integer(sizes[keep]),
        stringsAsFactors = FALSE,
        check.names = FALSE
    )
}

rna <- eligible(mapping_pairs(paths[["rna"]], "rna"))
protein <- eligible(mapping_pairs(paths[["protein"]], "protein"))
roster <- rbind(
    data.frame(assay = "bulk_rnaseq", rna, check.names = FALSE),
    data.frame(assay = "pseudobulk_rnaseq", rna, check.names = FALSE),
    data.frame(assay = "bulk_proteomics", protein, check.names = FALSE)
)
assay_order <- match(roster$assay, c(
    "bulk_rnaseq", "pseudobulk_rnaseq", "bulk_proteomics"
))
roster <- roster[order(assay_order, roster$target_id, method = "radix"), , drop = FALSE]
row.names(roster) <- NULL
if (anyDuplicated(roster[c("assay", "target_id")]) ||
    any(roster$declared_size < 10L | roster$declared_size > 200L)) {
    stop("Generated Reactome 97 roster is invalid.", call. = FALSE)
}

utils::write.table(
    roster, output, sep = "\t", row.names = FALSE, col.names = TRUE,
    quote = FALSE, na = "NA", eol = "\n"
)
counts <- table(factor(roster$assay, levels = c(
    "bulk_rnaseq", "pseudobulk_rnaseq", "bulk_proteomics"
)))
cat(
    paste(names(counts), counts, sep = "="),
    paste0("rows=", nrow(roster)),
    paste0("sha256=", unname(tools::sha256sum(output))),
    sep = "\n"
)
