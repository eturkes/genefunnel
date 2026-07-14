# Assisted-by: OpenAI Codex.

files <- list.files(
    recursive = TRUE,
    all.files = TRUE,
    full.names = TRUE
)
logs <- files[grepl("[.]Rcheck/00check[.]log$", files)]
if (length(logs) != 1L) {
    stop("Expected exactly one R CMD check log; found ", length(logs), ".")
}

status <- grep("^Status:", readLines(logs, warn = FALSE), value = TRUE)
if (length(status) != 1L || !identical(trimws(status), "Status: OK")) {
    stop("R CMD check did not report `Status: OK`.")
}
cat("R CMD check status is clean: ", normalizePath(logs), "\n", sep = "")
