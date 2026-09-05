# Run with the package and ggplot2 installed; no input data are downloaded.
script <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE))
if (length(script) != 1L || length(commandArgs(trailingOnly = TRUE))) {
    stop("Usage: Rscript tools/render-quickstart.R", call. = FALSE)
}
root <- dirname(dirname(normalizePath(script, mustWork = TRUE)))

render <- function(root) {
    package <- read.dcf(file.path(root, "DESCRIPTION"))[1L, "Package"]
    work <- tempfile(paste0(package, "-readme-"))
    dir.create(work)
    previous <- setwd(work)
    on.exit({
        setwd(previous)
        unlink(work, recursive = TRUE)
    }, add = TRUE)
    sys.source(file.path(root, "inst", "examples", "quickstart.R"),
        envir = new.env(parent = globalenv()))
    image <- list.files(file.path(work, paste0(package, "-results")),
        pattern = "\\.png$", full.names = TRUE)
    stopifnot(length(image) == 1L)
    destination <- file.path(root, "man", "figures")
    dir.create(destination, showWarnings = FALSE, recursive = TRUE)
    stopifnot(file.copy(image, file.path(destination, "quickstart.png"),
        overwrite = TRUE))
    message("README figure: ", file.path(destination, "quickstart.png"))
}
render(root)
