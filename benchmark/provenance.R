# Assisted-by: OpenAI Codex.

benchmark_command_output <- function(command, args = character(), wd = NULL) {
    old_wd <- getwd()
    if (!is.null(wd)) {
        setwd(wd)
        on.exit(setwd(old_wd), add = TRUE)
    }
    tryCatch(
        system2(command, args, stdout = TRUE, stderr = TRUE),
        error = function(error) character()
    )
}

benchmark_linux_field <- function(path, pattern) {
    lines <- tryCatch(readLines(path, warn = FALSE), error = function(error) character())
    value <- grep(pattern, lines, value = TRUE)
    if (length(value) == 0L) NA_character_ else trimws(sub(pattern, "", value[[1L]]))
}

benchmark_package_version <- function(package) {
    if (!requireNamespace(package, quietly = TRUE)) {
        return(NA_character_)
    }
    as.character(utils::packageVersion(package))
}

benchmark_metadata <- function(
    repo_root,
    runner,
    preset,
    repeats,
    workers,
    time_bin = ""
) {
    git_head <- benchmark_command_output("git", c("rev-parse", "HEAD"), repo_root)
    git_status <- benchmark_command_output("git", c("status", "--porcelain"), repo_root)
    cpu_model <- benchmark_linux_field(
        "/proc/cpuinfo",
        "^model name[[:space:]]*:[[:space:]]*"
    )
    memory_kib <- benchmark_linux_field(
        "/proc/meminfo",
        "^MemTotal:[[:space:]]*"
    )
    time_version <- if (nzchar(time_bin)) {
        version <- benchmark_command_output(time_bin, "--version")
        if (length(version)) version[[1L]] else NA_character_
    } else {
        NA_character_
    }
    data.frame(
        key = c(
            "generated_utc", "protocol_version", "runner", "preset",
            "repeats", "requested_workers", "git_head", "git_dirty",
            "hostname", "sysname", "release", "machine", "cpu_model",
            "logical_cores", "memory_total", "r_version", "r_platform",
            "genefunnel_version", "Matrix_version", "BiocParallel_version",
            "Rcpp_version", "RcppArmadillo_version", "rscript",
            "external_time"
        ),
        value = c(
            format(Sys.time(), tz = "UTC", usetz = TRUE),
            GENEFUNNEL_BENCHMARK_PROTOCOL,
            runner,
            preset,
            repeats,
            workers,
            if (length(git_head)) git_head[[1L]] else NA_character_,
            length(git_status) > 0L,
            unname(Sys.info()[["nodename"]]),
            unname(Sys.info()[["sysname"]]),
            unname(Sys.info()[["release"]]),
            unname(Sys.info()[["machine"]]),
            cpu_model,
            parallel::detectCores(logical = TRUE),
            memory_kib,
            paste(R.version$major, R.version$minor, sep = "."),
            R.version$platform,
            benchmark_package_version("genefunnel"),
            benchmark_package_version("Matrix"),
            benchmark_package_version("BiocParallel"),
            benchmark_package_version("Rcpp"),
            benchmark_package_version("RcppArmadillo"),
            Sys.which("Rscript"),
            time_version
        ),
        stringsAsFactors = FALSE,
        check.names = FALSE
    )
}

benchmark_write_session_info <- function(path, packages) {
    writeLines(
        capture.output({
            print(sessionInfo())
            cat("\nInstalled benchmark packages:\n")
            for (package in packages) {
                cat(package, benchmark_package_version(package), "\n")
            }
        }),
        path
    )
}
