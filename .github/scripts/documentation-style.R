# Assisted-by: OpenAI Codex.

.documentation_words <- function(text) {
    text <- gsub("`[^`]+`", " TOKEN ", text, perl = TRUE)
    text <- gsub("\\[([^]]+)\\]\\([^)]+\\)", "\\1", text, perl = TRUE)
    text <- gsub("[$][^$]+[$]", " TOKEN ", text, perl = TRUE)
    text <- gsub("[*][^*]+[*]", " TOKEN ", text, perl = TRUE)
    text <- gsub("<[^>]+>", " ", text, perl = TRUE)
    words <- regmatches(
        text,
        gregexpr("[[:alnum:]][[:alnum:]_/'’-]*", text, perl = TRUE)
    )[[1L]]
    if (identical(words, "")) character() else words
}

.documentation_instruction <- function(text) {
    verbs <- paste0(
        "Build|Calculate|Check|Choose|Commit|Declare|Deduplicate|Define|",
        "Execute|Exercise|Include|Inspect|Install|Keep|Label|Let|List|Load|",
        "Match|Open|Pair|Place|Prepare|Preserve|Record|Reject|Remove|Report|",
        "Return|Run|Score|Select|Start|Treat|Use|Validate|Vary"
    )
    grepl("\\b(must|shall|should)\\b", text, ignore.case = TRUE, perl = TRUE) ||
        grepl(paste0("^(?:", verbs, ")\\b"), text, perl = TRUE) ||
        grepl(
            paste0(
                "^(?:After|Before|If|Once|When)\\b[^,]{0,160},[[:space:]]*",
                "(?:", verbs, ")\\b"
            ),
            text,
            perl = TRUE
        )
}

.documentation_emit <- function(path, line, text) {
    text <- trimws(paste(text, collapse = " "))
    if (!nzchar(text)) return(NULL)
    sentences <- strsplit(text, "(?<=[.!?])\\s+", perl = TRUE)[[1L]]
    data.frame(
        path = rep(path, length(sentences)),
        line = rep(line, length(sentences)),
        text = trimws(sentences),
        stringsAsFactors = FALSE
    )
}

.documentation_blocks <- function(path) {
    if (identical(path, "DESCRIPTION")) {
        return(.documentation_emit(
            path,
            1L,
            read.dcf(path, fields = "Description")[[1L]]
        ))
    }

    lines <- readLines(path, warn = FALSE)
    is_r <- grepl("[.]R$", path)
    yaml <- grepl("[.]Rmd$", path) && length(lines) > 0L &&
        identical(lines[[1L]], "---")
    fenced <- math <- html <- roxygen_skip <- FALSE
    block <- character()
    block_line <- 1L
    result <- list()
    flush <- function() {
        if (length(block)) {
            result[[length(result) + 1L]] <<-
                .documentation_emit(path, block_line, block)
            block <<- character()
        }
    }

    for (index in seq_along(lines)) {
        line <- lines[[index]]
        trimmed <- trimws(line)

        if (yaml) {
            if (index > 1L && identical(trimmed, "---")) yaml <- FALSE
            next
        }
        if (grepl("^[[:space:]]*(```|~~~)", line)) {
            flush()
            fenced <- !fenced
            next
        }
        if (fenced) next
        if (!is_r &&
            (startsWith(trimmed, "$$") || startsWith(trimmed, "\\["))) {
            flush()
            math <- !(nchar(trimmed) > 4L && endsWith(trimmed, "$$")) &&
                !(nchar(trimmed) > 4L && endsWith(trimmed, "\\]"))
            next
        }
        if (!is_r && math) {
            if (identical(trimmed, "\\]") || endsWith(trimmed, "$$")) {
                math <- FALSE
            }
            next
        }
        if (html || grepl("^[[:space:]]*<!--", line)) {
            flush()
            html <- !grepl("-->", line)
            next
        }

        if (is_r) {
            if (!grepl("^[[:space:]]*#'", line)) {
                flush()
                roxygen_skip <- FALSE
                next
            }
            content <- sub("^[[:space:]]*#' ?", "", line)
            if (math) {
                if (identical(trimws(content), "}")) math <- FALSE
                next
            }
            if (startsWith(trimws(content), "\\deqn{")) {
                flush()
                math <- TRUE
                next
            }
            if (grepl("^@(examples|references|seealso)\\b", content)) {
                flush()
                roxygen_skip <- TRUE
                next
            }
            if (grepl("^@", content)) {
                roxygen_skip <- FALSE
                if (grepl(
                    "^@(aliases|export|inheritParams|import|keywords|useDynLib)\\b",
                    content
                )) {
                    flush()
                    next
                }
                content <- sub("^@[[:alnum:]_]+[[:space:]]*", "", content)
            }
            if (roxygen_skip) next
            line <- content
            trimmed <- trimws(line)
        }

        if (!nzchar(trimmed) ||
            grepl("^[[:space:]]*#+[[:space:]]", line) ||
            grepl("^[[:space:]]*\\|", line) ||
            grepl("^[[:space:]]*[-*_]{3,}[[:space:]]*$", line)) {
            flush()
            next
        }
        if (grepl("^[[:space:]]*[-*+] ", line) ||
            grepl("^[[:space:]]*[0-9]+[.] ", line)) {
            flush()
            block_line <- index
            line <- sub("^[[:space:]]*([-*+] |[0-9]+[.] )", "", line)
        } else if (!length(block)) {
            block_line <- index
        }
        block <- c(block, line)
    }
    flush()
    do.call(rbind, result)
}

documentation_style_check <- function() {
    paths <- c(
        "DESCRIPTION",
        "NEWS.md",
        "README.md",
        "benchmark/README.md",
        list.files("benchmark", "-result[.]md$", full.names = TRUE),
        "formal/README.md",
        list.files("inst", "[.]md$", full.names = TRUE),
        list.files("R", "[.]R$", full.names = TRUE),
        "vignettes/genefunnel.Rmd"
    )
    paths <- sort(unique(paths))
    missing <- paths[!file.exists(paths)]
    if (length(missing)) {
        stop("Human-facing documentation is missing: ", paste(missing, collapse = ", "))
    }
    sentences <- do.call(rbind, lapply(paths, .documentation_blocks))
    counts <- vapply(sentences$text, function(text) {
        length(.documentation_words(text))
    }, integer(1L))
    limits <- ifelse(
        vapply(sentences$text, .documentation_instruction, logical(1L)),
        20L,
        25L
    )
    too_long <- counts > limits
    filler <- grepl(
        "\\b(simply|robust|seamlessly|leverage)\\b",
        sentences$text,
        ignore.case = TRUE,
        perl = TRUE
    )
    invalid <- too_long | filler
    if (any(invalid)) {
        details <- sprintf(
            "%s:%d [%d/%d words]%s %s",
            sentences$path[invalid],
            sentences$line[invalid],
            counts[invalid],
            limits[invalid],
            ifelse(filler[invalid], " [filler]", ""),
            sentences$text[invalid]
        )
        stop("Human-facing documentation style violations:\n", paste(details, collapse = "\n"))
    }
    invisible(length(sentences$text))
}
