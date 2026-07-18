# Assisted-by: OpenAI Codex.

if (!nzchar(Sys.getenv("TZ", unset = ""))) {
    Sys.setenv(TZ = "UTC")
}

runner_file <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(runner_file) != 1L) {
    stop("Cannot locate benchmark/prepare-components.R.", call. = FALSE)
}
benchmark_dir <- dirname(normalizePath(sub("^--file=", "", runner_file)))
repo_root <- dirname(benchmark_dir)
source(file.path(benchmark_dir, "components-installation.R"))

COMPONENT_PROTOCOL_VERSION <- "A2-1.0.0"
COMPONENT_BASELINE_ID <- "9b60a3eb138e5fd267586624ccd8bf51907577e7"

components_prepare_usage <- function() {
    cat(paste0(
        "Usage: Rscript benchmark/prepare-components.R --output=DIR\n",
        "Creates exact baseline/candidate source snapshots and package libraries.\n"
    ))
}

components_prepare_options <- function(args) {
    if (identical(args, "--help")) {
        components_prepare_usage()
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
    parent <- normalizePath(dirname(output), mustWork = TRUE)
    output <- file.path(parent, basename(output))
    if (!startsWith(output, root_prefix)) {
        stop("`--output` resolves outside the repository.", call. = FALSE)
    }
    output
}

components_prepare_installation <- function(root, method, git_sha) {
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
        stop("Prepared source archive lacks DESCRIPTION.", call. = FALSE)
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
        stop("Failed to install the ", method, " package; see ", log, ".", call. = FALSE)
    }

    package_version <- utils::packageDescription(
        "genefunnel",
        lib.loc = library
    )[["Version"]]
    provenance <- data.frame(
        protocol_version = COMPONENT_PROTOCOL_VERSION,
        method = method,
        git_sha = git_sha,
        source_archive_md5 = unname(tools::md5sum(archive)),
        package_version = package_version,
        installed_manifest_md5 = components_manifest_digest(package_path),
        stringsAsFactors = FALSE,
        check.names = FALSE
    )
    utils::write.table(
        provenance,
        file.path(library, COMPONENT_INSTALLATION_MARKER),
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
        source_archive_md5 = provenance$source_archive_md5,
        package_version = package_version,
        installed_manifest_md5 = provenance$installed_manifest_md5,
        stringsAsFactors = FALSE,
        check.names = FALSE
    )
}

output <- components_prepare_options(commandArgs(trailingOnly = TRUE))
candidate_id <- system2(
    "git",
    c("-C", shQuote(repo_root), "rev-parse", "HEAD"),
    stdout = TRUE,
    stderr = FALSE
)
if (length(candidate_id) != 1L || !grepl("^[0-9a-f]{40}$", candidate_id)) {
    stop("Cannot resolve the candidate Git SHA.", call. = FALSE)
}
dir.create(output)
installations <- rbind(
    components_prepare_installation(output, "baseline", COMPONENT_BASELINE_ID),
    components_prepare_installation(output, "candidate", candidate_id)
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
cat("Prepared component libraries: ", output, "\n", sep = "")
