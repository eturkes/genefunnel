# Assisted-by: OpenAI Codex.

if (!nzchar(Sys.getenv("TZ", unset = ""))) {
    Sys.setenv(TZ = "UTC")
}

runner_file <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(runner_file) != 1L) {
    stop("Cannot locate run-sensitivity-profile.R.", call. = FALSE)
}
benchmark_dir <- dirname(normalizePath(sub("^--file=", "", runner_file)))
repo_root <- dirname(benchmark_dir)
source(file.path(benchmark_dir, "protocol.R"))
source(file.path(benchmark_dir, "provenance.R"))
source(file.path(benchmark_dir, "components-installation.R"))
source(file.path(benchmark_dir, "sensitivity-protocol.R"))
source(file.path(benchmark_dir, "sensitivity-profile.R"))

sensitivity_profile_usage <- function() {
    cat(paste0(
        "Usage: Rscript benchmark/run-sensitivity-profile.R [options]\n",
        "  --mode=smoke|gate       Downscaled smoke or frozen profile\n",
        "  --candidate-id=SHA      Full clean Git SHA required for gate\n",
        "  --output=DIR            New result directory inside repository\n",
        "  --help                  Show this text\n"
    ))
}

sensitivity_profile_options <- function(args) {
    options <- list(mode = "smoke", candidate_id = NULL, output = NULL)
    for (arg in args) {
        if (identical(arg, "--help")) {
            sensitivity_profile_usage()
            quit(save = "no", status = 0L)
        }
        matched <- FALSE
        for (key in names(options)) {
            flag <- paste0("--", gsub("_", "-", key), "=")
            if (startsWith(arg, flag)) {
                options[[key]] <- substring(arg, nchar(flag) + 1L)
                matched <- TRUE
                break
            }
        }
        if (!matched) stop("Unknown sensitivity profile option: ", arg)
    }
    if (!options$mode %in% c("smoke", "gate")) {
        stop("`--mode` must be smoke or gate.", call. = FALSE)
    }
    options
}

sensitivity_profile_output <- function(path, mode) {
    if (is.null(path)) {
        stamp <- format(Sys.time(), "%Y%m%d-%H%M%S", tz = "UTC")
        path <- file.path("benchmark", "results", paste(stamp, "sensitivity", mode, sep = "-"))
    }
    root <- normalizePath(repo_root, mustWork = TRUE)
    if (!startsWith(path, .Platform$file.sep)) {
        path <- file.path(repo_root, path)
    }
    path <- normalizePath(path, mustWork = FALSE)
    if (!startsWith(path, paste0(root, .Platform$file.sep))) {
        stop("Sensitivity profile output must be inside the repository.")
    }
    if (file.exists(path)) {
        stop("Sensitivity profile output already exists: ", path)
    }
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    dir.create(path)
    normalizePath(path, mustWork = TRUE)
}

sensitivity_profile_git <- function(options) {
    head <- benchmark_command_output("git", c("rev-parse", "HEAD"), repo_root)
    status <- benchmark_command_output("git", c("status", "--porcelain"), repo_root)
    if (length(head) != 1L || !grepl("^[0-9a-f]{40}$", head)) {
        stop("Cannot resolve sensitivity profile Git SHA.", call. = FALSE)
    }
    candidate <- if (is.null(options$candidate_id)) head[[1L]] else options$candidate_id
    if (options$mode == "gate" &&
        (is.null(options$candidate_id) || !grepl("^[0-9a-f]{40}$", candidate) ||
            !identical(candidate, head[[1L]]) || length(status) > 0L)) {
        stop("Gate mode requires the explicit clean current Git SHA.", call. = FALSE)
    }
    list(head = head[[1L]], candidate = candidate, dirty = length(status) > 0L)
}

sensitivity_profile_install <- function(output_dir, git_sha) {
    archive <- file.path(output_dir, "source.tar")
    source_dir <- file.path(output_dir, "source")
    library <- file.path(output_dir, "library")
    log <- file.path(output_dir, "install.log")
    dir.create(source_dir)
    dir.create(library)
    archive_status <- system2(
        "git",
        c("-C", shQuote(repo_root), "archive", "--format=tar",
            shQuote(paste0("--output=", archive)), shQuote(git_sha)),
        stdout = log,
        stderr = log
    )
    if (!identical(archive_status, 0L)) stop("Cannot create profile source archive.")
    utils::untar(archive, exdir = source_dir)
    install_status <- system2(
        file.path(R.home("bin"), "R"),
        c("CMD", "INSTALL", "--no-multiarch", "--with-keep.source",
            shQuote(paste0("--library=", library)), shQuote(source_dir)),
        stdout = log,
        stderr = log
    )
    package <- file.path(library, "genefunnel")
    if (!identical(install_status, 0L) || !dir.exists(package)) {
        stop("Cannot install profile package; see ", log, ".")
    }
    list(
        library = normalizePath(library),
        archive_md5 = unname(tools::md5sum(archive)),
        archive_sha256 = unname(tools::sha256sum(archive)),
        manifest_md5 = components_manifest_digest(package)
    )
}

sensitivity_profile_metadata <- function(git, install, options) {
    system <- Sys.info()
    data.frame(
        key = c(
            "generated_utc", "profile_protocol", "parent_protocol", "mode",
            "candidate_id", "git_head", "git_dirty", "hostname", "sysname",
            "release", "machine", "cpu_model", "logical_cores", "memory_total",
            "r_version", "r_platform", "source_archive_md5",
            "source_archive_sha256", "installed_manifest_md5", "rscript"
        ),
        value = c(
            format(Sys.time(), tz = "UTC", usetz = TRUE),
            SENSITIVITY_PROFILE_VERSION, "E-1.0.0", options$mode,
            git$candidate, git$head, git$dirty, unname(system[["nodename"]]),
            unname(system[["sysname"]]), unname(system[["release"]]),
            unname(system[["machine"]]), benchmark_linux_field(
                "/proc/cpuinfo", "^model name[[:space:]]*:[[:space:]]*"
            ), parallel::detectCores(logical = TRUE), benchmark_linux_field(
                "/proc/meminfo", "^MemTotal:[[:space:]]*"
            ), paste(R.version$major, R.version$minor, sep = "."),
            R.version$platform, install$archive_md5, install$archive_sha256,
            install$manifest_md5, Sys.which("Rscript")
        ),
        stringsAsFactors = FALSE,
        check.names = FALSE
    )
}

sensitivity_profile_report <- function(output_dir, git) {
    runs <- utils::read.delim(file.path(output_dir, "runs.tsv"), check.names = FALSE)
    profile <- utils::read.delim(file.path(output_dir, "profile.tsv"), check.names = FALSE)
    allocation <- utils::read.delim(
        file.path(output_dir, "allocation.tsv"), check.names = FALSE
    )
    decision <- utils::read.delim(file.path(output_dir, "decision.tsv"), check.names = FALSE)
    run_rows <- paste0(
        "| ", runs$repeat_id, " | ", formatC(runs$elapsed_sec, digits = 6, format = "f"),
        " | `", runs$output_md5, "` |"
    )
    lines <- c(
        "# Exact sensitivity profile", "",
        paste0("Candidate: `", git$candidate, "`; protocol `", SENSITIVITY_PROFILE_VERSION,
            "` supplementing `E-1.0.0`."), "",
        "| Repeat | Elapsed seconds | Output MD5 |", "|---:|---:|---|", run_rows, "",
        paste0("Rprof: ", profile$total_samples, " samples; exact-arithmetic share = ",
            formatC(profile$exact_sample_share, digits = 6, format = "f"), "."),
        paste0("Rprofmem: ", allocation$manager_r_alloc_bytes, " bytes across ",
            allocation$manager_r_alloc_events, " events."), "",
        paste0("Optimization eligible: **", decision$optimization_eligible, "** (elapsed >60 s: ",
            decision$elapsed_trigger, "; exact share >0.50: ", decision$exact_trigger, ")."),
        "Performance is descriptive; this decision makes no reliability or public-API claim."
    )
    writeLines(lines, file.path(output_dir, "report.md"))
}

profile_write_tsv <- function(value, path) {
    utils::write.table(
        value, path, sep = "\t", row.names = FALSE, col.names = TRUE,
        quote = TRUE, na = "NA"
    )
}

options <- sensitivity_profile_options(commandArgs(trailingOnly = TRUE))
git <- sensitivity_profile_git(options)
output_dir <- sensitivity_profile_output(options$output, options$mode)
protocol <- sensitivity_profile_validate(repo_root)
parent <- sensitivity_validate_protocol(repo_root)
registry <- sensitivity_read_registry(file.path(benchmark_dir, "sensitivity-protocol.tsv"))
config <- sensitivity_profile_config(registry, options$mode)
install <- sensitivity_profile_install(output_dir, git$candidate)
saveRDS(config, file.path(output_dir, "config.rds"), version = 3L)
invisible(file.copy(
    file.path(benchmark_dir, "sensitivity-profile-protocol.tsv"),
    file.path(output_dir, "protocol.tsv")
))
invisible(file.copy(
    file.path(benchmark_dir, "sensitivity-protocol.tsv"),
    file.path(output_dir, "parent-protocol.tsv")
))
profile_write_tsv(
    sensitivity_profile_metadata(git, install, options),
    file.path(output_dir, "metadata.tsv")
)

worker <- file.path(benchmark_dir, "sensitivity-profile-worker.R")
stdout <- file.path(output_dir, "worker-stdout.log")
stderr <- file.path(output_dir, "worker-stderr.log")
paths <- unique(c(install$library, normalizePath(.libPaths(), mustWork = TRUE)))
status <- system2(
    Sys.which("Rscript"),
    c("--vanilla", shQuote(worker), shQuote(file.path(output_dir, "config.rds")),
        shQuote(output_dir), shQuote(install$library)),
    stdout = stdout,
    stderr = stderr,
    env = c(
        paste0("R_LIBS=", paste(paths, collapse = .Platform$path.sep)),
        "GENEFUNNEL_BENCHMARK_TOKEN=sensitivity-profile"
    )
)
if (!identical(status, 0L)) {
    logs <- c(readLines(stdout, warn = FALSE), readLines(stderr, warn = FALSE))
    stop("Sensitivity profile worker failed:\n", paste(tail(logs, 40L), collapse = "\n"))
}
sensitivity_profile_report(output_dir, git)
core_artifacts <- c(
    "protocol.tsv", "parent-protocol.tsv", "metadata.tsv", "manifest.tsv",
    "runs.tsv", "profile.tsv", "allocation.tsv", "decision.tsv",
    "worker-environment.tsv", "session-info.txt", "report.md"
)
profile_write_tsv(
    data.frame(
        path = core_artifacts,
        sha256 = unname(tools::sha256sum(file.path(output_dir, core_artifacts))),
        stringsAsFactors = FALSE,
        check.names = FALSE
    ),
    file.path(output_dir, "artifact-hashes.tsv")
)
cat(readLines(stdout, warn = FALSE), sep = "\n")
cat("Sensitivity profile complete: ", output_dir, "\n", sep = "")
