# Run from any directory; dependencies must already be installed.
args <- commandArgs(trailingOnly = TRUE)
if (any(!args %in% "--bioc") || anyDuplicated(args)) {
    stop("Usage: Rscript tools/check-package.R [--bioc]", call. = FALSE)
}
script <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE))
if (length(script) != 1L) {
    stop("Run this file with Rscript.", call. = FALSE)
}
root <- dirname(dirname(normalizePath(script, mustWork = TRUE)))
description <- read.dcf(file.path(root, "DESCRIPTION"))
package <- description[1L, "Package"]
version <- description[1L, "Version"]
output <- file.path(root, ".checks")
dir.create(output, showWarnings = FALSE)
work <- tempfile(paste0(package, "-"), tmpdir = output)
dir.create(work)
work <- normalizePath(work, mustWork = TRUE)
message("Check output: ", work)
if (!nzchar(Sys.getenv("TZ"))) {
    Sys.setenv(TZ = "UTC")
}
writeLines(capture.output(sessionInfo()), file.path(work, "session-info.log"))
r <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "R.exe" else "R")

run <- function(arguments, logfile) {
    logfile <- file.path(work, logfile)
    status <- system2(r, arguments, stdout = logfile, stderr = logfile)
    if (!identical(status, 0L)) {
        cat(tail(readLines(logfile, warn = FALSE), 60L), sep = "\n")
        stop("Command failed; inspect ", logfile, call. = FALSE)
    }
    invisible(logfile)
}

setwd(work)
run(c("CMD", "build", shQuote(root)), "build.log")
archive <- file.path(work, paste0(package, "_", version, ".tar.gz"))
if (!file.exists(archive)) {
    stop("The source build did not create the expected archive.", call. = FALSE)
}
run(c("CMD", "check", "--no-manual", shQuote(archive)), "check.log")
checklog <- file.path(work, paste0(package, ".Rcheck"), "00check.log")
lines <- readLines(checklog, warn = FALSE)
if (!any(grepl("^Status: OK[[:space:]]*$", lines))) {
    cat(tail(lines, 60L), sep = "\n")
    stop("R CMD check did not report Status: OK.", call. = FALSE)
}

if ("--bioc" %in% args) {
    if (!requireNamespace("BiocCheck", quietly = TRUE)) {
        stop("Install BiocCheck before using --bioc.", call. = FALSE)
    }
    expression <- paste(
        "a <- commandArgs(TRUE);",
        "result <- BiocCheck::BiocCheck(a[[1L]], checkDir = a[[2L]],",
        "`new-package` = TRUE);",
        "if (length(result$error) || length(result$warning)) quit(status = 1L)"
    )
    run(c("--vanilla", "--slave", "-e", shQuote(expression), "--args",
        shQuote(archive), shQuote(work)), "bioccheck.log")
}
message(package, " ", version, ": all requested checks passed.")
