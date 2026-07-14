# Assisted-by: OpenAI Codex.

expected <- Sys.getenv("BIOC_VERSION", unset = "")
if (!grepl("^[0-9]+[.][0-9]+$", expected)) {
    stop("BIOC_VERSION must contain a major.minor Bioconductor version.")
}

installed <- unclass(utils::packageVersion("BiocVersion"))[[1L]]
actual <- paste(installed[seq_len(2L)], collapse = ".")
managed <- as.character(BiocManager::version())
if (!identical(actual, expected) || !identical(managed, expected)) {
    stop(
        "Expected Bioconductor ",
        expected,
        "; BiocVersion reports ",
        actual,
        " and BiocManager reports ",
        managed,
        "."
    )
}

cat("Bioconductor version is pinned to ", expected, ".\n", sep = "")
