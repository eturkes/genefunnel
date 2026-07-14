# Assisted-by: OpenAI Codex.

arguments <- commandArgs(trailingOnly = TRUE)
if (length(arguments) != 1L || !arguments[[1L]] %in% c("git", "tarball")) {
    stop("Usage: Rscript .github/scripts/bioc-check.R git|tarball")
}

assert_clean <- function(result, label) {
    counts <- vapply(
        c("error", "warning"),
        function(level) length(result[[level]]),
        integer(1)
    )
    if (any(counts > 0L)) {
        stop(
            label,
            " produced ",
            counts[["error"]],
            " error(s) and ",
            counts[["warning"]],
            " warning(s)."
        )
    }
    invisible(result)
}

if (identical(arguments[[1L]], "git")) {
    result <- BiocCheck::BiocCheckGitClone(
        ".",
        `quit-with-status` = FALSE
    )
    assert_clean(result, "BiocCheckGitClone")
    quit(save = "no", status = 0L)
}

tarballs <- list.files(
    pattern = "^genefunnel_.*[.]tar[.]gz$",
    full.names = TRUE
)
if (length(tarballs) != 1L) {
    stop("Expected exactly one GeneFunnel source tarball.")
}

# Submission numbering and maintainer-owned Support-site state are explicit
# Session 12 gates. All package-controlled checks remain enabled here.
result <- BiocCheck::BiocCheck(
    tarballs[[1L]],
    `new-package` = TRUE,
    `no-check-version-num` = TRUE,
    `no-check-bioc-help` = TRUE,
    `quit-with-status` = FALSE
)
assert_clean(result, "BiocCheck")
