# Assisted-by: OpenAI Codex.

COMPONENT_INSTALLATION_MARKER <- ".genefunnel-components-installation.tsv"

components_git_archive_digest <- function(repo_root, git_sha) {
    archive <- tempfile("genefunnel-components-source-", fileext = ".tar")
    on.exit(unlink(archive), add = TRUE)
    status <- system2(
        "git",
        c(
            "-C",
            shQuote(repo_root),
            "archive",
            "--format=tar",
            shQuote(paste0("--output=", archive)),
            shQuote(git_sha)
        ),
        stdout = FALSE,
        stderr = FALSE
    )
    if (!identical(status, 0L) || !file.exists(archive)) {
        stop("Cannot fingerprint the component source snapshot.", call. = FALSE)
    }
    unname(tools::md5sum(archive))
}

components_installed_manifest <- function(package_path) {
    package_path <- normalizePath(package_path, mustWork = TRUE)
    files <- list.files(
        package_path,
        recursive = TRUE,
        full.names = TRUE,
        all.files = TRUE,
        no.. = TRUE
    )
    files <- sort(files[!dir.exists(files)])
    relative <- substring(files, nchar(package_path) + 2L)
    data.frame(
        path = relative,
        md5 = unname(tools::md5sum(files)),
        stringsAsFactors = FALSE,
        check.names = FALSE
    )
}

components_manifest_digest <- function(package_path) {
    manifest <- components_installed_manifest(package_path)
    path <- tempfile("genefunnel-components-manifest-", fileext = ".tsv")
    on.exit(unlink(path), add = TRUE)
    utils::write.table(
        manifest,
        path,
        sep = "\t",
        row.names = FALSE,
        col.names = TRUE,
        quote = TRUE,
        na = "NA"
    )
    unname(tools::md5sum(path))
}

components_read_installation <- function(library) {
    marker <- file.path(library, COMPONENT_INSTALLATION_MARKER)
    if (!file.exists(marker)) {
        stop("Component package library lacks installation provenance.", call. = FALSE)
    }
    provenance <- utils::read.delim(
        marker,
        stringsAsFactors = FALSE,
        check.names = FALSE
    )
    required <- c(
        "protocol_version", "method", "git_sha", "source_archive_md5",
        "package_version", "installed_manifest_md5"
    )
    if (!identical(names(provenance), required) || nrow(provenance) != 1L) {
        stop("Component installation provenance is malformed.", call. = FALSE)
    }
    provenance
}
