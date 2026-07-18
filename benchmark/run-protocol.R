# Assisted-by: OpenAI Codex.

if (!nzchar(Sys.getenv("TZ", unset = ""))) {
    Sys.setenv(TZ = "UTC")
}

harness_file <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(harness_file) != 1L) {
    stop("Cannot locate benchmark/run-protocol.R.", call. = FALSE)
}
benchmark_dir <- dirname(normalizePath(sub("^--file=", "", harness_file)))
repo_root <- dirname(benchmark_dir)
source(file.path(benchmark_dir, "protocol-index.R"))

benchmark_protocol_usage <- function() {
    cat(paste0(
        "Usage: Rscript benchmark/run-protocol.R [harness options] [runner options]\n",
        "  --protocol=VERSION  Exact protocol version; no implicit latest\n",
        "  --suite=SUITE       Exact executable suite\n",
        "  --list              List registered protocol/suite pairs\n",
        "  --help              Show this text\n",
        "All other options pass unchanged to the selected frozen runner.\n"
    ))
}

benchmark_protocol_parse_options <- function(args) {
    options <- list(protocol = NULL, suite = NULL, list = FALSE, runner = character())
    assign_once <- function(current, value, label) {
        if (!is.null(current)) stop("Duplicate harness option: ", label, call. = FALSE)
        value
    }
    for (arg in args) {
        if (identical(arg, "--help")) {
            benchmark_protocol_usage()
            quit(save = "no", status = 0L)
        } else if (identical(arg, "--list")) {
            if (options$list) stop("Duplicate harness option: --list", call. = FALSE)
            options$list <- TRUE
        } else if (grepl("^--protocol=", arg)) {
            options$protocol <- assign_once(
                options$protocol, sub("^--protocol=", "", arg), "--protocol"
            )
        } else if (grepl("^--suite=", arg)) {
            options$suite <- assign_once(
                options$suite, sub("^--suite=", "", arg), "--suite"
            )
        } else {
            options$runner <- c(options$runner, arg)
        }
    }
    options
}

benchmark_protocol_print_available <- function(suites) {
    view <- suites
    names(view) <- c("protocol", "suite")
    utils::write.table(
        view, stdout(), sep = "\t", row.names = FALSE, quote = FALSE
    )
}

options <- benchmark_protocol_parse_options(commandArgs(trailingOnly = TRUE))
validated <- benchmark_protocol_validate_index(repo_root)
if (options$list) {
    if (!is.null(options$protocol) || !is.null(options$suite) ||
        length(options$runner) > 0L) {
        stop("--list cannot be combined with execution options.", call. = FALSE)
    }
    benchmark_protocol_print_available(validated$suites)
    quit(save = "no", status = 0L)
}
if (is.null(options$protocol) || is.null(options$suite)) {
    stop("Explicit --protocol and --suite are required.", call. = FALSE)
}
selected <- benchmark_protocol_select(
    validated$index, options$protocol, options$suite
)
runner <- benchmark_protocol_runner(selected, repo_root)
rscript <- Sys.which("Rscript")
if (!nzchar(rscript)) stop("Rscript is unavailable.", call. = FALSE)
cat("Validated protocol ", options$protocol, "/", options$suite, ".\n", sep = "")
status <- system2(
    rscript,
    c(
        "--vanilla", shQuote(runner),
        vapply(options$runner, shQuote, character(1L))
    )
)
if (!identical(status, 0L)) {
    quit(save = "no", status = as.integer(status))
}
