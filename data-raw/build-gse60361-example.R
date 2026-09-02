#!/usr/bin/env Rscript

required_packages <- c("Matrix", "scRNAseq", "SummarizedExperiment")
missing_packages <- required_packages[!vapply(
    required_packages,
    requireNamespace,
    logical(1L),
    quietly = TRUE
)]
if (length(missing_packages)) {
    stop(
        "Install the data-build dependencies: ",
        paste(missing_packages, collapse = ", "),
        call. = FALSE
    )
}

expected_fixture_md5 <- "93e3022388eefd47dc06dfcad9c194aa"

sce <- scRNAseq::ZeiselBrainData(
    ensembl = FALSE,
    location = FALSE,
    legacy = FALSE
)
counts <- SummarizedExperiment::assay(sce, "counts")
row_data <- as.data.frame(SummarizedExperiment::rowData(sce))
cell_data <- as.data.frame(SummarizedExperiment::colData(sce))
stopifnot(
    identical(dim(sce), c(20006L, 3005L)),
    identical(dim(counts), c(20006L, 3005L)),
    identical(
        as.integer(table(row_data$featureType)[c("endogenous", "mito")]),
        c(19972L, 34L)
    ),
    all(c("level1class", "level2class") %in% names(cell_data)),
    !anyDuplicated(rownames(counts)),
    all(Matrix::colSums(counts) > 0)
)

endogenous <- row_data$featureType == "endogenous"
annotated <- cell_data$level2class != "(none)"
counts <- counts[endogenous, annotated, drop = FALSE]
cell_data <- cell_data[annotated, , drop = FALSE]
stopifnot(
    identical(dim(counts), c(19972L, 2816L)),
    length(unique(cell_data$level1class)) == 7L,
    length(unique(cell_data$level2class)) == 47L
)

aggregate_counts <- function(labels) {
    levels <- sort(unique(labels))
    indices <- match(labels, levels)
    design <- Matrix::sparseMatrix(
        i = seq_along(indices),
        j = indices,
        x = 1L,
        dims = c(length(indices), length(levels)),
        dimnames = list(colnames(counts), levels)
    )
    aggregated <- as.matrix(counts %*% design)
    stopifnot(
        all(is.finite(aggregated)),
        all(aggregated >= 0),
        all(aggregated == floor(aggregated)),
        max(aggregated) <= .Machine$integer.max
    )
    storage.mode(aggregated) <- "integer"
    aggregated
}

broad_class_counts <- aggregate_counts(cell_data$level1class)
subtype_counts <- aggregate_counts(cell_data$level2class)
subtype_levels <- colnames(subtype_counts)
broad_class_levels <- colnames(broad_class_counts)

subtype_metadata <- unique(cell_data[c("level1class", "level2class")])
subtype_metadata <- subtype_metadata[
    match(subtype_levels, subtype_metadata$level2class),
    ,
    drop = FALSE
]
subtype_metadata <- data.frame(
    subtype = subtype_metadata$level2class,
    broad_class = subtype_metadata$level1class,
    cells = as.integer(table(factor(
        cell_data$level2class,
        levels = subtype_levels
    ))),
    check.names = FALSE,
    stringsAsFactors = FALSE
)

broad_class_metadata <- data.frame(
    broad_class = broad_class_levels,
    cells = as.integer(table(factor(
        cell_data$level1class,
        levels = broad_class_levels
    ))),
    subtypes = vapply(
        broad_class_levels,
        function(class_name) {
            length(unique(cell_data$level2class[
                cell_data$level1class == class_name
            ]))
        },
        integer(1L)
    ),
    check.names = FALSE,
    stringsAsFactors = FALSE
)

fixture <- list(
    subtype_counts = subtype_counts,
    broad_class_counts = broad_class_counts,
    subtype_metadata = subtype_metadata,
    broad_class_metadata = broad_class_metadata,
    source = list(
        accession = "GSE60361",
        dataset = "zeisel-brain-2015",
        dataset_version = "2023-12-14",
        original_dimensions = c(genes = 20006L, cells = 3005L),
        endogenous_genes = 19972L,
        annotated_cells = 2816L
    )
)

stopifnot(
    identical(dim(fixture$subtype_counts), c(19972L, 47L)),
    identical(dim(fixture$broad_class_counts), c(19972L, 7L)),
    identical(
        fixture$subtype_metadata$subtype,
        colnames(subtype_counts)
    ),
    identical(
        fixture$broad_class_metadata$broad_class,
        colnames(broad_class_counts)
    ),
    sum(fixture$subtype_metadata$cells) == 2816L,
    sum(fixture$broad_class_metadata$cells) == 2816L
)

output <- file.path(
    "inst",
    "extdata",
    "GSE60361_brain_cell_types.rds"
)
dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
saveRDS(fixture, output, compress = "xz", version = 2)
fixture_md5 <- unname(tools::md5sum(output))
stopifnot(identical(fixture_md5, expected_fixture_md5))
message(output, "  md5=", fixture_md5)
