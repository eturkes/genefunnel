# Assisted-by: OpenAI Codex.

if (!nzchar(Sys.getenv("TZ", unset = ""))) {
    Sys.setenv(TZ = "UTC")
}

runner_file <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(runner_file) != 1L) {
    stop("Cannot locate benchmark/prepare-catalogue.R.", call. = FALSE)
}
benchmark_dir <- dirname(normalizePath(sub("^--file=", "", runner_file)))
repo_root <- dirname(benchmark_dir)
source(file.path(benchmark_dir, "catalogue-installation.R"))

CATALOGUE_PROTOCOL_VERSION <- "C-1.0.2"
CATALOGUE_BASELINE_ID <- "a573c124909235e41bdbc3cfae950947465d8755"

catalogue_prepare_options <- function(args) {
    if (identical(args, "--help")) {
        cat(paste0(
            "Usage: Rscript benchmark/prepare-catalogue.R --output=DIR\n",
            "Creates exact baseline/candidate package snapshots and libraries.\n"
        ))
        quit(save = "no", status = 0L)
    }
    if (length(args) != 1L || !startsWith(args, "--output=")) {
        stop("Exactly one `--output=DIR` argument is required.", call. = FALSE)
    }
    output <- substring(args, nchar("--output=") + 1L)
    if (!nzchar(output)) {
        stop("`--output` must not be empty.", call. = FALSE)
    }
    output <- normalizePath(output, mustWork = FALSE)
    root <- normalizePath(repo_root, mustWork = TRUE)
    root_prefix <- paste0(root, .Platform$file.sep)
    if (!startsWith(output, root_prefix)) {
        stop("`--output` must be inside the repository.", call. = FALSE)
    }
    if (file.exists(output)) {
        stop("`--output` already exists.", call. = FALSE)
    }
    dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
    output <- file.path(normalizePath(dirname(output), mustWork = TRUE), basename(output))
    if (!startsWith(output, root_prefix)) {
        stop("`--output` resolves outside the repository.", call. = FALSE)
    }
    output
}

catalogue_prepare_installation <- function(root, method, git_sha) {
    source_dir <- file.path(root, paste0(method, "-source"))
    library <- file.path(root, paste0(method, "-library"))
    archive <- file.path(root, paste0(method, "-source.tar"))
    log <- file.path(root, paste0(method, "-install.log"))
    dir.create(source_dir)
    dir.create(library)
    archive_status <- system2(
        "git",
        c(
            "-C",
            shQuote(repo_root),
            "archive",
            "--format=tar",
            shQuote(paste0("--output=", archive)),
            shQuote(git_sha)
        ),
        stdout = log,
        stderr = log
    )
    if (!identical(archive_status, 0L) || !file.exists(archive)) {
        stop("Failed to create the ", method, " source archive.", call. = FALSE)
    }
    utils::untar(archive, exdir = source_dir)
    if (!file.exists(file.path(source_dir, "DESCRIPTION"))) {
        stop("Prepared catalogue source lacks DESCRIPTION.", call. = FALSE)
    }
    install_status <- system2(
        file.path(R.home("bin"), "R"),
        c(
            "CMD", "INSTALL", "--no-multiarch", "--with-keep.source",
            shQuote(paste0("--library=", library)),
            shQuote(source_dir)
        ),
        stdout = log,
        stderr = log
    )
    package_path <- file.path(library, "genefunnel")
    if (!identical(install_status, 0L) || !dir.exists(package_path)) {
        stop("Failed to install ", method, "; see ", log, ".", call. = FALSE)
    }
    package_version <- utils::packageDescription(
        "genefunnel",
        lib.loc = library
    )[["Version"]]
    provenance <- data.frame(
        protocol_version = CATALOGUE_PROTOCOL_VERSION,
        method = method,
        git_sha = git_sha,
        source_archive_sha256 = unname(tools::sha256sum(archive)),
        package_version = package_version,
        installed_manifest_sha256 = catalogue_installed_digest(package_path),
        stringsAsFactors = FALSE,
        check.names = FALSE
    )
    utils::write.table(
        provenance,
        file.path(library, CATALOGUE_INSTALLATION_MARKER),
        sep = "\t",
        row.names = FALSE,
        col.names = TRUE,
        quote = TRUE,
        na = "NA"
    )
    data.frame(
        method = method,
        library = normalizePath(library, mustWork = TRUE),
        git_sha = git_sha,
        source_archive_sha256 = provenance$source_archive_sha256,
        package_version = package_version,
        installed_manifest_sha256 = provenance$installed_manifest_sha256,
        stringsAsFactors = FALSE,
        check.names = FALSE
    )
}

output <- catalogue_prepare_options(commandArgs(trailingOnly = TRUE))
candidate_id <- system2(
    "git",
    c("-C", shQuote(repo_root), "rev-parse", "HEAD"),
    stdout = TRUE,
    stderr = FALSE
)
if (length(candidate_id) != 1L || !grepl("^[0-9a-f]{40}$", candidate_id)) {
    stop("Cannot resolve the catalogue candidate Git SHA.", call. = FALSE)
}
dir.create(output)
installations <- rbind(
    catalogue_prepare_installation(output, "list", CATALOGUE_BASELINE_ID),
    catalogue_prepare_installation(output, "compiled", candidate_id)
)
utils::write.table(
    installations,
    file.path(output, "installations.tsv"),
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE,
    quote = TRUE,
    na = "NA"
)
cat("Prepared catalogue libraries: ", output, "\n", sep = "")
