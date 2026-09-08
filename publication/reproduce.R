# Source package → isolated install → full analysis → independent replay → article.
script <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE))
root <- dirname(dirname(normalizePath(script, mustWork = TRUE)))
setwd(root)
args <- commandArgs(trailingOnly = TRUE)
stopifnot(all(args %in% c("--no-replay", "--no-render", "--replay-only")),
    !anyDuplicated(args), !("--no-replay" %in% args && "--replay-only" %in% args))
.libPaths(c(file.path(root, ".checks/publication-library"),
    file.path(dirname(root), ".scratch/publication/R-library"),
    file.path(dirname(root), ".R-library"), .libPaths()))
Sys.setenv(R_LIBS = paste(.libPaths(), collapse = .Platform$path.sep), TZ = "UTC")
required <- c("Rcpp", "RcppArmadillo", "Matrix", "BiocParallel", "ggplot2",
    "airway", "SummarizedExperiment", "singscore", "GSVA", "rmarkdown", "knitr")
missing <- required[!vapply(required, function(p) nzchar(system.file(package = p)), logical(1))]
if (length(missing)) stop("Install publication dependencies: ", paste(missing, collapse = ", "))
work <- file.path(root, ".checks/publication")
dir.create(work, recursive = TRUE, showWarnings = FALSE)
library <- file.path(root, ".checks/publication-library")
dir.create(library, recursive = TRUE, showWarnings = FALSE)
run <- function(program, arguments, name) {
    logfile <- file.path(work, paste0(name, ".log"))
    status <- system2(program, arguments, stdout = logfile, stderr = logfile)
    if (status != 0L) {
        cat(tail(readLines(logfile, warn = FALSE), 40), sep = "\n")
        stop(name, " failed; log=", logfile)
    }
    message(name, ": passed")
}
mapping_files <- c("publication/inputs/wang2004-activated-hgnc.tsv",
    "publication/inputs/hgnc-provenance.json")
mapping_before <- tools::md5sum(mapping_files)
python <- Sys.which("python3")
if (!nzchar(python)) stop("Install Python 3 to replay the offline HGNC mapping.")
run(python, "publication/inputs/map-wang2004.py", "mapping")
stopifnot(identical(mapping_before, tools::md5sum(mapping_files)))
rscript <- file.path(R.home("bin"), "Rscript")
if (!"--replay-only" %in% args) {
    run(file.path(R.home("bin"), "R"), c("CMD", "INSTALL",
        paste0("--library=", shQuote(library)), shQuote(root)), "install")
    run(rscript, "publication/analyze.R", "analysis")
} else {
    stopifnot(file.exists(file.path(library, "genefunnel", "DESCRIPTION")))
}
source("publication/replay.R")
run(rscript, "publication/check-replay.R", "replay-regression")
if (!"--no-replay" %in% args) {
    replay <- new_replay_directory(work)
    message("Fresh replay directory: ", replay)
    Sys.setenv(GENEFUNNEL_SKIP_TIMING = "1")
    run(rscript, c("publication/analyze.R", shQuote(file.path(replay, "results")),
        shQuote(file.path(replay, "figures"))), "replay")
    Sys.unsetenv("GENEFUNNEL_SKIP_TIMING")
    comparison <- compare_replay("publication", replay)
    utils::write.csv(comparison, "publication/results/replay.csv", row.names = FALSE)
    stopifnot(all(comparison$replay_identical))
    message("Replay: ", nrow(comparison), " byte-identical tables/figures")
    unlink(replay, recursive = TRUE)
}
if (!"--no-render" %in% args) {
    rmarkdown::render("publication/article.Rmd", output_format = "all", quiet = TRUE)
    message("Article: HTML + DOCX")
}
