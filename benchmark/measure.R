benchmark_rss_kib <- function(pid) {
    path <- sprintf("/proc/%d/status", as.integer(pid))
    lines <- tryCatch(readLines(path, warn = FALSE), error = function(error) character())
    field <- grep("^VmRSS:[[:space:]]", lines, value = TRUE)
    if (length(field) != 1L) {
        return(NA_real_)
    }
    as.numeric(sub("^VmRSS:[[:space:]]*([0-9]+).*$", "\\1", field))
}

benchmark_child_pids <- function(pid) {
    path <- sprintf(
        "/proc/%d/task/%d/children",
        as.integer(pid),
        as.integer(pid)
    )
    children <- tryCatch(readLines(path, warn = FALSE), error = function(error) "")
    children <- trimws(paste(children, collapse = " "))
    if (!nzchar(children)) {
        return(integer())
    }
    suppressWarnings(as.integer(strsplit(children, "[[:space:]]+")[[1L]]))
}

benchmark_process_tree <- function(root, exclude = integer()) {
    queued <- as.integer(root)
    seen <- integer()
    while (length(queued) > 0L) {
        pid <- queued[[1L]]
        queued <- queued[-1L]
        if (pid %in% seen || pid %in% exclude) {
            next
        }
        seen <- c(seen, pid)
        queued <- c(queued, benchmark_child_pids(pid))
    }
    seen
}

benchmark_proc_text <- function(pid, field) {
    path <- sprintf("/proc/%d/%s", as.integer(pid), field)
    value <- tryCatch(
        suppressWarnings(readBin(path, what = "raw", n = 1048576L)),
        error = function(error) raw()
    )
    if (length(value) == 0L) {
        return("")
    }
    value[value == as.raw(0)] <- charToRaw("\n")
    rawToChar(value)
}

benchmark_tagged_processes <- function(token, exclude = integer()) {
    pids <- suppressWarnings(as.integer(list.files(
        "/proc",
        pattern = "^[0-9]+$"
    )))
    pids <- pids[is.finite(pids) & !pids %in% exclude]
    marker <- paste0("GENEFUNNEL_BENCHMARK_TOKEN=", token)
    tagged <- vapply(pids, function(pid) {
        grepl(marker, benchmark_proc_text(pid, "environ"), fixed = TRUE)
    }, logical(1))
    pids <- pids[tagged]
    commands <- vapply(pids, benchmark_proc_text, character(1), field = "cmdline")
    pids[!grepl("(^|/)time([[:space:]\n]|$)", commands)]
}

benchmark_start_memory_monitor <- function(interval_sec = 0.01) {
    root <- Sys.getpid()
    baseline <- benchmark_rss_kib(root)
    token <- Sys.getenv("GENEFUNNEL_BENCHMARK_TOKEN", unset = "")
    supported <- .Platform$OS.type == "unix" &&
        dir.exists("/proc/self") &&
        is.finite(baseline)
    if (!supported) {
        return(list(
            supported = FALSE,
            baseline_rss_kib = baseline,
            scope = "unavailable"
        ))
    }

    done <- tempfile("genefunnel-memory-monitor-")
    job <- parallel::mcparallel({
        monitor_pid <- Sys.getpid()
        peak_rss_kib <- baseline
        peak_processes <- 1L
        samples <- 0L
        known_pids <- root
        repeat {
            if (nzchar(token) && samples %% 5L == 0L) {
                known_pids <- benchmark_tagged_processes(
                    token,
                    exclude = monitor_pid
                )
            }
            pids <- if (nzchar(token)) {
                unique(c(root, known_pids))
            } else {
                benchmark_process_tree(root, exclude = monitor_pid)
            }
            rss <- vapply(pids, benchmark_rss_kib, numeric(1))
            aggregate_rss <- sum(rss, na.rm = TRUE)
            if (aggregate_rss > peak_rss_kib) {
                peak_rss_kib <- aggregate_rss
                peak_processes <- sum(is.finite(rss))
            }
            samples <- samples + 1L
            if (file.exists(done) || !dir.exists(sprintf("/proc/%d", root))) {
                break
            }
            Sys.sleep(interval_sec)
        }
        c(
            score_peak_tree_rss_kib = peak_rss_kib,
            score_peak_processes = peak_processes,
            memory_samples = samples
        )
    }, silent = TRUE, mc.set.seed = FALSE)

    list(
        supported = TRUE,
        baseline_rss_kib = baseline,
        scope = if (nzchar(token)) "tagged-process aggregate" else "descendant aggregate",
        interval_sec = interval_sec,
        done = done,
        job = job
    )
}

benchmark_stop_memory_monitor <- function(monitor) {
    if (!isTRUE(monitor$supported)) {
        return(c(
            score_baseline_rss_kib = monitor$baseline_rss_kib,
            score_peak_tree_rss_kib = NA_real_,
            score_peak_increment_kib = NA_real_,
            score_peak_processes = NA_real_,
            memory_samples = NA_real_
        ))
    }

    file.create(monitor$done)
    result <- parallel::mccollect(monitor$job, wait = TRUE)[[1L]]
    unlink(monitor$done)
    if (inherits(result, "try-error")) {
        warning("Process-tree memory monitor failed: ", result, call. = FALSE)
        result <- c(
            score_peak_tree_rss_kib = NA_real_,
            score_peak_processes = NA_real_,
            memory_samples = NA_real_
        )
    }
    peak <- unname(result[["score_peak_tree_rss_kib"]])
    c(
        score_baseline_rss_kib = monitor$baseline_rss_kib,
        score_peak_tree_rss_kib = peak,
        score_peak_increment_kib = if (is.finite(peak)) {
            max(0, peak - monitor$baseline_rss_kib)
        } else {
            NA_real_
        },
        score_peak_processes = unname(result[["score_peak_processes"]]),
        memory_samples = unname(result[["memory_samples"]])
    )
}

benchmark_allocation_stats <- function(path) {
    lines <- if (file.exists(path)) readLines(path, warn = FALSE) else character()
    allocation_lines <- grep("^[0-9]+[[:space:]]*:", lines, value = TRUE)
    bytes <- suppressWarnings(as.numeric(sub(":.*$", "", allocation_lines)))
    bytes <- bytes[is.finite(bytes)]
    c(
        manager_r_alloc_bytes = sum(bytes),
        manager_r_alloc_events = length(bytes)
    )
}

benchmark_output_md5 <- function(value) {
    path <- tempfile("genefunnel-output-", fileext = ".rds")
    on.exit(unlink(path), add = TRUE)
    saveRDS(value, path, compress = FALSE, version = 3L)
    unname(tools::md5sum(path))
}
