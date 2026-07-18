# Assisted-by: OpenAI Codex.

aggregation_clean_git_head <- function(root) {
    head <- benchmark_command_output("git", c("rev-parse", "HEAD"), root)
    status <- benchmark_command_output(
        "git",
        c("status", "--porcelain", "--untracked-files=normal"),
        root
    )
    if (length(head) != 1L || !grepl("^[0-9a-f]{40}$", head[[1L]])) {
        stop("Cannot resolve the aggregation runner Git commit.", call. = FALSE)
    }
    if (length(status) > 0L) {
        stop("Full aggregation execution requires a clean tree.", call. = FALSE)
    }
    head[[1L]]
}

aggregation_publish_locked_file <- function(source, target, expected_hash) {
    if (file.exists(target)) {
        if (!identical(unname(tools::sha256sum(target)), expected_hash)) {
            stop("Existing aggregation protocol artifact differs.", call. = FALSE)
        }
        return(invisible(target))
    }
    if (!file.copy(source, target, overwrite = FALSE) ||
        !identical(unname(tools::sha256sum(target)), expected_hash)) {
        stop("Cannot publish aggregation protocol artifact.", call. = FALSE)
    }
    invisible(target)
}

aggregation_publish_locked_tsv <- function(value, target) {
    temporary <- tempfile("aggregation-manifest-", tmpdir = dirname(target))
    on.exit(unlink(temporary, force = TRUE), add = TRUE)
    benchmark_write_tsv(value, temporary)
    if (file.exists(target)) {
        if (!identical(
            unname(tools::sha256sum(target)),
            unname(tools::sha256sum(temporary))
        )) {
            stop("Existing aggregation manifest differs.", call. = FALSE)
        }
    } else if (!file.rename(temporary, target)) {
        stop("Cannot publish aggregation manifest.", call. = FALSE)
    }
    invisible(target)
}

aggregation_install_package <- function(root, output_dir) {
    library <- tempfile("aggregation-library-")
    dir.create(library)
    published <- FALSE
    on.exit({
        if (!published) unlink(library, recursive = TRUE, force = TRUE)
    }, add = TRUE)
    log <- file.path(output_dir, "installation.log")
    status <- system2(
        file.path(R.home("bin"), "R"),
        c(
            "CMD", "INSTALL", "--preclean", "--clean",
            shQuote(paste0("--library=", library)), shQuote(root)
        ),
        stdout = log,
        stderr = log
    )
    if (!identical(status, 0L)) {
        stop("Committed package installation failed; see installation.log.",
            call. = FALSE)
    }
    .libPaths(c(library, .libPaths()))
    namespace <- loadNamespace("genefunnel", lib.loc = library)
    result <- list(
        library = library, namespace = namespace,
        audit = get(".aggregation_audit", envir = namespace, inherits = FALSE),
        scorer = get("genefunnel", envir = namespace, inherits = FALSE)
    )
    published <- TRUE
    result
}

aggregation_artifact_manifest <- function(output_dir, files) {
    paths <- file.path(output_dir, files)
    if (any(!file.exists(paths))) {
        stop("Aggregation evidence artifacts are incomplete.", call. = FALSE)
    }
    data.frame(
        artifact = files, bytes = as.numeric(file.info(paths)$size),
        sha256 = unname(tools::sha256sum(paths)), stringsAsFactors = FALSE,
        check.names = FALSE
    )
}

aggregation_external_file_map <- function(data_dir) {
    basenames <- c(
        cellbench_celseq2_counts = "RNAmix_celseq2.count.csv.gz",
        cellbench_celseq2_metadata = "RNAmix_celseq2.metadata.csv.gz",
        cellbench_sortseq_counts = "RNAmix_sortseq.count.csv.gz",
        cellbench_sortseq_metadata = "RNAmix_sortseq.metadata.csv.gz",
        kang_batch2_counts = "GSE96583_RAW.tar",
        kang_batch2_genes = "GSE96583_batch2.genes.tsv.gz",
        kang_batch2_metadata = "GSE96583_batch2.total.tsne.df.tsv.gz",
        reactome_v97 = "ReactomePathways-v97.gmt.zip"
    )
    stats::setNames(file.path(data_dir, unname(basenames)), names(basenames))
}
