# Assisted-by: OpenAI Codex.

markdown_files <- c(
    "README.md",
    list.files("inst", pattern = "[.]md$", recursive = TRUE, full.names = TRUE),
    list.files(
        "vignettes",
        pattern = "[.](Rmd|md)$",
        recursive = TRUE,
        full.names = TRUE
    )
)

missing_links <- character()
for (path in markdown_files) {
    text <- paste(readLines(path, warn = FALSE), collapse = "\n")
    matches <- regmatches(
        text,
        gregexpr("]\\([^)]+\\)", text, perl = TRUE)
    )[[1L]]
    if (identical(matches, character(0)) || identical(matches, "")) {
        next
    }
    targets <- sub("^]\\(([^)]+)\\)$", "\\1", matches)
    targets <- sub("[[:space:]].*$", "", targets)
    local <- !grepl("^[A-Za-z][A-Za-z0-9+.-]*:", targets) &
        !startsWith(targets, "#")
    targets <- targets[local]
    for (target in targets) {
        target <- utils::URLdecode(sub("[?#].*$", "", target))
        resolved <- file.path(dirname(path), target)
        if (!file.exists(resolved)) {
            missing_links <- c(missing_links, paste0(path, " -> ", target))
        }
    }
}
if (length(missing_links) > 0L) {
    stop("Missing local documentation links:\n", paste(missing_links, collapse = "\n"))
}

rd_files <- list.files("man", pattern = "[.]Rd$", full.names = TRUE)
rd_problems <- lapply(rd_files, tools::checkRd)
names(rd_problems) <- rd_files
rd_problems <- rd_problems[lengths(rd_problems) > 0L]
if (length(rd_problems) > 0L) {
    for (path in names(rd_problems)) {
        cat(path, ":\n", sep = "")
        print(rd_problems[[path]])
    }
    stop("Rd validation reported problems.")
}

rendered <- tempfile(fileext = ".html")
on.exit(unlink(rendered, force = TRUE), add = TRUE)
rmarkdown::render(
    "vignettes/genefunnel.Rmd",
    output_file = basename(rendered),
    output_dir = dirname(rendered),
    envir = new.env(parent = globalenv()),
    quiet = TRUE,
    clean = TRUE
)
if (!file.exists(rendered) || file.info(rendered)$size < 1000) {
    stop("The rendered vignette is missing or unexpectedly small.")
}

html <- paste(readLines(rendered, warn = FALSE), collapse = "\n")
remote_resource <- paste0(
    "<(script|img|iframe)[^>]+src=[\"']https?://|",
    "<link[^>]+href=[\"']https?://"
)
if (grepl(remote_resource, html, ignore.case = TRUE, perl = TRUE)) {
    stop("The standalone vignette references an external rendering resource.")
}

cat(
    "Documentation checks passed: Rd, local links, and standalone vignette.\n"
)
