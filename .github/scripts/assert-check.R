# Assisted-by: OpenAI Codex.

logs <- Sys.glob(file.path("*.Rcheck", "00check.log"))
if (length(logs) != 1L) {
    stop(
        "Expected exactly one top-level R CMD check log; found ",
        length(logs),
        "."
    )
}

status <- grep("^Status:", readLines(logs, warn = FALSE), value = TRUE)
if (length(status) != 1L || !identical(trimws(status), "Status: OK")) {
    stop("R CMD check did not report `Status: OK`.")
}
cat("R CMD check status is clean: ", normalizePath(logs), "\n", sep = "")
