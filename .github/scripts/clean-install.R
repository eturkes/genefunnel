# Assisted-by: OpenAI Codex.

tarballs <- list.files(
    pattern = "^genefunnel_.*[.]tar[.]gz$",
    full.names = TRUE
)
if (length(tarballs) != 1L) {
    stop("Expected exactly one GeneFunnel source tarball.")
}

clean_library <- tempfile("genefunnel-clean-library-")
dir.create(clean_library)
on.exit(unlink(clean_library, recursive = TRUE, force = TRUE), add = TRUE)

status <- system2(
    file.path(R.home("bin"), "R"),
    c(
        "CMD",
        "INSTALL",
        "--install-tests",
        paste0("--library=", clean_library),
        tarballs[[1L]]
    )
)
if (!identical(status, 0L)) {
    stop("Clean-library source-tarball installation failed.")
}

.libPaths(c(clean_library, .libPaths()))
package_path <- normalizePath(find.package("genefunnel"))
if (!identical(dirname(package_path), normalizePath(clean_library))) {
    stop("GeneFunnel was not loaded from the clean target library.")
}

# SOCK workers are fresh R processes: propagate the temporary target library
# through their startup environment so they load this tarball installation.
user_libraries <- strsplit(
    Sys.getenv("R_LIBS_USER", unset = ""),
    .Platform$path.sep,
    fixed = TRUE
)[[1L]]
worker_libraries <- unique(c(
    clean_library,
    user_libraries[nzchar(user_libraries)]
))
Sys.setenv(R_LIBS_USER = paste(worker_libraries, collapse = .Platform$path.sep))

worker_paths <- unlist(BiocParallel::bplapply(
    seq_len(2L),
    function(...) normalizePath(find.package("genefunnel")),
    BPPARAM = BiocParallel::SnowParam(workers = 2L, type = "SOCK")
))
if (!all(worker_paths == package_path)) {
    stop("SOCK workers did not load GeneFunnel from the clean target library.")
}

testthat::test_package("genefunnel", reporter = "summary")
cat("Installed-tarball tests passed from ", package_path, ".\n", sep = "")
