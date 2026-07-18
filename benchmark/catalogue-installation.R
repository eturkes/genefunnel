# Assisted-by: OpenAI Codex.

CATALOGUE_INSTALLATION_MARKER <- ".genefunnel-catalogue-installation.tsv"

catalogue_git_archive_digest <- function(repo_root, git_sha) {
    archive <- tempfile("genefunnel-catalogue-source-", fileext = ".tar")
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
        stop("Cannot fingerprint the catalogue source snapshot.", call. = FALSE)
    }
    unname(tools::sha256sum(archive))
}

catalogue_installed_manifest <- function(package_path) {
    package_path <- normalizePath(package_path, mustWork = TRUE)
    files <- list.files(
        package_path,
        recursive = TRUE,
        full.names = TRUE,
        all.files = TRUE,
        no.. = TRUE
    )
    files <- sort(files[!dir.exists(files)])
    data.frame(
        path = substring(files, nchar(package_path) + 2L),
        sha256 = unname(tools::sha256sum(files)),
        stringsAsFactors = FALSE,
        check.names = FALSE
    )
}

catalogue_installed_digest <- function(package_path) {
    manifest <- catalogue_installed_manifest(package_path)
    bytes <- serialize(manifest, NULL, ascii = FALSE, xdr = TRUE, version = 3L)
    unname(tools::sha256sum(bytes = bytes))
}

catalogue_read_installation <- function(library) {
    marker <- file.path(library, CATALOGUE_INSTALLATION_MARKER)
    if (!file.exists(marker)) {
        stop("Catalogue package library lacks installation provenance.", call. = FALSE)
    }
    provenance <- utils::read.delim(
        marker,
        stringsAsFactors = FALSE,
        check.names = FALSE
    )
    required <- c(
        "protocol_version", "method", "git_sha", "source_archive_sha256",
        "package_version", "installed_manifest_sha256"
    )
    if (!identical(names(provenance), required) || nrow(provenance) != 1L) {
        stop("Catalogue installation provenance is malformed.", call. = FALSE)
    }
    provenance
}

catalogue_package_record <- function(
    repo_root,
    library,
    method,
    expected_id,
    protocol_version,
    require_marker
) {
    library <- normalizePath(library, mustWork = TRUE)
    package_path <- normalizePath(file.path(library, "genefunnel"), mustWork = TRUE)
    description <- utils::packageDescription("genefunnel", lib.loc = library)
    marker <- file.path(library, CATALOGUE_INSTALLATION_MARKER)
    provenance <- if (file.exists(marker)) {
        catalogue_read_installation(library)
    } else if (require_marker) {
        stop("Gate package library lacks catalogue provenance.", call. = FALSE)
    } else {
        NULL
    }
    if (require_marker) {
        archive <- file.path(dirname(library), paste0(method, "-source.tar"))
        archive_digest <- if (file.exists(archive)) {
            unname(tools::sha256sum(archive))
        } else {
            NA_character_
        }
        valid <- provenance$protocol_version[[1L]] == protocol_version &&
            provenance$method[[1L]] == method &&
            provenance$git_sha[[1L]] == expected_id &&
            provenance$source_archive_sha256[[1L]] == archive_digest &&
            provenance$source_archive_sha256[[1L]] ==
                catalogue_git_archive_digest(repo_root, expected_id) &&
            provenance$package_version[[1L]] == description[["Version"]] &&
            provenance$installed_manifest_sha256[[1L]] ==
                catalogue_installed_digest(package_path)
        if (!isTRUE(valid)) {
            stop("Catalogue package installation provenance is invalid.", call. = FALSE)
        }
    }
    data.frame(
        method = method,
        library = library,
        package_path = package_path,
        package_version = description[["Version"]],
        git_sha = if (is.null(provenance)) NA_character_ else provenance$git_sha,
        installed_manifest_sha256 = if (is.null(provenance)) {
            NA_character_
        } else {
            provenance$installed_manifest_sha256
        },
        stringsAsFactors = FALSE,
        check.names = FALSE
    )
}
