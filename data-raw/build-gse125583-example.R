#!/usr/bin/env Rscript

series_url <- paste0(
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE125nnn/GSE125583/",
    "matrix/GSE125583_series_matrix.txt.gz"
)
raw_url <- paste0(
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE125nnn/GSE125583/",
    "suppl/GSE125583_RAW.tar"
)
expected_md5 <- c(
    GSE125583_series_matrix.txt.gz = "2a56c2cce8448e1fdfaddfa8f34482cb",
    GSE125583_RAW.tar = "04a3223237c7d0bfc6d94aca76b5830f"
)
expected_fixture_md5 <- "4f3038f701f2015f35f46c938e3b0cdd"
gene_ids <- c(
    SYN1 = "6853",
    SYN2 = "6854",
    SYP = "6855",
    SNAP25 = "6616",
    STX1A = "6804",
    VAMP2 = "6844"
)

cache_dir <- Sys.getenv(
    "GENEFUNNEL_GSE125583_CACHE",
    file.path(tempdir(), "genefunnel-gse125583")
)
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

fetch <- function(url, filename) {
    path <- file.path(cache_dir, filename)
    valid <- file.exists(path) &&
        identical(unname(tools::md5sum(path)), expected_md5[[filename]])
    if (!valid) {
        utils::download.file(url, path, mode = "wb", quiet = FALSE)
    }
    stopifnot(
        identical(unname(tools::md5sum(path)), expected_md5[[filename]])
    )
    path
}

series_path <- fetch(series_url, "GSE125583_series_matrix.txt.gz")
raw_path <- fetch(raw_url, "GSE125583_RAW.tar")

series_connection <- gzfile(series_path, open = "rt")
series_lines <- tryCatch(
    readLines(series_connection, warn = FALSE),
    finally = close(series_connection)
)

parse_values <- function(line) {
    values <- strsplit(line, "\t", fixed = TRUE)[[1L]][-1L]
    sub('^"(.*)"$', "\\1", values)
}

single_row <- function(prefix) {
    rows <- series_lines[startsWith(series_lines, prefix)]
    stopifnot(length(rows) == 1L)
    parse_values(rows)
}

accessions <- single_row("!Sample_geo_accession")
characteristic_rows <- series_lines[
    startsWith(series_lines, "!Sample_characteristics_ch1")
]
characteristics <- lapply(characteristic_rows, parse_values)
characteristic_names <- vapply(
    characteristics,
    function(values) sub(":.*$", "", values[[1L]]),
    character(1L)
)
names(characteristics) <- characteristic_names
characteristics <- lapply(
    characteristics,
    function(values) sub("^[^:]+: ", "", values)
)

braak_levels <- c(I = 1L, II = 2L, III = 3L, IV = 4L, V = 5L, VI = 6L)
braak <- unname(braak_levels[characteristics[["braak.score"]]])
age <- suppressWarnings(as.integer(characteristics[["age"]]))

members <- utils::untar(raw_path, list = TRUE)
member_accessions <- sub("_.*$", "", basename(members))
members <- members[match(accessions, member_accessions)]
stopifnot(length(members) == 289L, !anyNA(members))

extract_dir <- tempfile("gse125583-")
dir.create(extract_dir)
on.exit(unlink(extract_dir, recursive = TRUE, force = TRUE), add = TRUE)
utils::untar(raw_path, files = members, exdir = extract_dir)

read_sample <- function(member) {
    connection <- gzfile(file.path(extract_dir, member), open = "rt")
    sample_table <- tryCatch(
        utils::read.delim(
            connection,
            comment.char = "#",
            check.names = FALSE,
            colClasses = c("character", "numeric", "integer")
        ),
        finally = close(connection)
    )
    rows <- match(unname(gene_ids), sample_table$ID_REF)
    stopifnot(!anyNA(rows), all(sample_table$count[rows] >= 0L))
    stats::setNames(sample_table$VALUE[rows], names(gene_ids))
}

expression <- vapply(members, read_sample, numeric(length(gene_ids)))
stopifnot(
    identical(rownames(expression), names(gene_ids)),
    identical(ncol(expression), length(accessions)),
    all(is.finite(expression)),
    all(expression >= 0)
)

fixture <- data.frame(
    sample = accessions,
    case_id = characteristics[["case id"]],
    diagnosis = characteristics[["diagnosis"]],
    age = age,
    sex = characteristics[["Sex"]],
    braak = braak,
    t(expression),
    check.names = FALSE,
    stringsAsFactors = FALSE
)
stopifnot(
    identical(dim(fixture), c(289L, 12L)),
    identical(
        as.integer(table(fixture$diagnosis)[
            c("Alzheimer's disease", "control")
        ]),
        c(219L, 70L)
    ),
    sum(is.na(fixture$braak)) == 2L
)

output <- file.path("inst", "extdata", "GSE125583_presynaptic_release.csv")
dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
utils::write.csv(fixture, output, row.names = FALSE, na = "")
fixture_md5 <- unname(tools::md5sum(output))
stopifnot(identical(fixture_md5, expected_fixture_md5))
message(output, "  md5=", fixture_md5)
