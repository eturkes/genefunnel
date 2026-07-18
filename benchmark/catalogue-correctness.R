# Assisted-by: OpenAI Codex.

if (!nzchar(Sys.getenv("TZ", unset = ""))) {
    Sys.setenv(TZ = "UTC")
}

runner_file <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(runner_file) != 1L) {
    stop("Cannot locate benchmark/catalogue-correctness.R.", call. = FALSE)
}
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
    stop(
        "Usage: catalogue-correctness.R CANDIDATE_LIBRARY RESULT_TSV",
        call. = FALSE
    )
}
candidate_library <- normalizePath(args[[1L]], mustWork = TRUE)
result_path <- args[[2L]]
.libPaths(c(candidate_library, .libPaths()))
Sys.setenv(R_LIBS = paste(.libPaths(), collapse = .Platform$path.sep))
for (package in c("genefunnel", "Matrix", "BiocParallel")) {
    if (!requireNamespace(package, quietly = TRUE)) {
        stop("Catalogue correctness dependency is unavailable: ", package, call. = FALSE)
    }
}
expected_package <- normalizePath(
    file.path(candidate_library, "genefunnel"),
    mustWork = TRUE
)
if (!identical(normalizePath(find.package("genefunnel")), expected_package)) {
    stop("Catalogue correctness loaded the wrong package.", call. = FALSE)
}

compile_catalogue <- getFromNamespace(".compile_gene_sets", "genefunnel")
score_compiled <- getFromNamespace(".genefunnel_compiled", "genefunnel")
serialize_catalogue <- getFromNamespace(
    ".serialize_compiled_catalogue",
    "genefunnel"
)
unserialize_catalogue <- getFromNamespace(
    ".unserialize_compiled_catalogue",
    "genefunnel"
)
validate_catalogue <- getFromNamespace(
    ".validate_compiled_catalogue",
    "genefunnel"
)
serial <- BiocParallel::SerialParam(progressbar = FALSE)
mat <- matrix(
    c(4, 4, 2, 0, 0, 0, 2, 8, 4, NA_real_, 2, 1, 1, 3, NaN, 5),
    nrow = 4L,
    dimnames = list(c("A", "B", "C", "D"), paste0("sample_", seq_len(4L)))
)
gene_sets <- list(
    full = c("A", "B", "A"),
    partial = c("C", "absent"),
    overlap = c("B", "C", "D"),
    empty = character()
)
compiled <- compile_catalogue(gene_sets, rownames(mat))
bytes <- serialize_catalogue(compiled)

capture_call <- function(expression) {
    warnings <- character()
    value <- withCallingHandlers(
        expression,
        warning = function(condition) {
            warnings <<- c(warnings, conditionMessage(condition))
            invokeRestart("muffleWarning")
        }
    )
    list(value = value, warnings = warnings)
}

capture_error <- function(expression) {
    tryCatch(expression, error = identity)
}

error_message <- function(value) {
    if (inherits(value, "error")) conditionMessage(value) else "<no error>"
}

rows <- list()
record <- function(assertion, passed, observed) {
    rows[[length(rows) + 1L]] <<- data.frame(
        assertion = assertion,
        passed = isTRUE(passed),
        observed = as.character(observed),
        stringsAsFactors = FALSE,
        check.names = FALSE
    )
}

list_call <- capture_call(genefunnel::genefunnel(mat, gene_sets, BPPARAM = serial))
compiled_call <- capture_call(score_compiled(mat, compiled, serial))
record(
    "serial_list_identity",
    identical(compiled_call$value, list_call$value),
    "compiled dense = named-list dense"
)
record(
    "warning_identity",
    identical(compiled_call$warnings, list_call$warnings) &&
        length(list_call$warnings) == 1L,
    paste(compiled_call$warnings, collapse = " | ")
)
sparse_call <- capture_call(score_compiled(
    Matrix::Matrix(mat, sparse = TRUE),
    compiled,
    serial
))
record(
    "sparse_identity",
    identical(sparse_call$value, list_call$value) &&
        identical(sparse_call$warnings, list_call$warnings),
    "compiled sparse = named-list dense"
)
restored <- unserialize_catalogue(bytes)
record(
    "wire_roundtrip",
    identical(restored, compiled) &&
        identical(serialize_catalogue(restored), bytes),
    unname(tools::sha256sum(bytes = bytes))
)
append_digest <- getFromNamespace(".catalogue_append_digest", "genefunnel")
body <- bytes[seq_len(length(bytes) - 32L)]
corrupt <- bytes
corrupt[[20L]] <- as.raw(bitwXor(as.integer(corrupt[[20L]]), 1L))
schema_body <- body
schema_body[[9L]] <- as.raw(2L)
formula_body <- body
formula_body[[15L]] <- charToRaw("X")
raw_variants <- list(
    corrupt = corrupt,
    truncated = bytes[seq_len(32L)],
    trailing = append_digest(c(body, as.raw(0L))),
    schema = append_digest(schema_body),
    formula = append_digest(formula_body)
)
raw_errors <- lapply(raw_variants, function(payload) {
    capture_error(unserialize_catalogue(payload))
})
raw_messages <- vapply(raw_errors, error_message, character(1))
raw_patterns <- c("fingerprint", "truncated", "trailing", "schema", "formula")
record(
    "malformed_wire_rejection",
    all(vapply(raw_errors, inherits, logical(1), "error")) &&
        all(mapply(
            grepl,
            raw_patterns,
            raw_messages,
            MoreArgs = list(ignore.case = TRUE),
            USE.NAMES = FALSE
        )),
    paste(raw_messages, collapse = " | ")
)

stale_mat <- mat[c(2L, 1L, 3L, 4L), , drop = FALSE]
stale_error <- capture_error(score_compiled(stale_mat, compiled, serial))
record(
    "stale_order_rejection",
    inherits(stale_error, "error") &&
        grepl("feature identities and order", error_message(stale_error)),
    error_message(stale_error)
)
stale_value <- mat
rownames(stale_value)[[4L]] <- "changed"
stale_value_error <- capture_error(
    score_compiled(stale_value, compiled, serial)
)
record(
    "stale_identity_rejection",
    inherits(stale_value_error, "error") &&
        grepl(
            "feature identities and order",
            error_message(stale_value_error)
        ),
    error_message(stale_value_error)
)
duplicate_error <- capture_error(compile_catalogue(
    gene_sets,
    c("A", "A", "C", "D")
))
record(
    "duplicate_feature_rejection",
    inherits(duplicate_error, "error") &&
        grepl("duplicat|unique", error_message(duplicate_error)),
    error_message(duplicate_error)
)
tampered <- compiled
tampered$member_index$full[[1L]] <- 2L
tamper_error <- capture_error(validate_catalogue(tampered))
record(
    "tamper_rejection",
    inherits(tamper_error, "error") &&
        grepl("malformed", error_message(tamper_error)),
    error_message(tamper_error)
)
schema_object <- compiled
schema_object$schema_version <- 2L
formula_object <- compiled
formula_object$formula_version <- "changed"
version_errors <- lapply(
    list(schema_object, formula_object),
    function(object) capture_error(validate_catalogue(object))
)
version_messages <- vapply(version_errors, error_message, character(1))
record(
    "object_version_rejection",
    all(vapply(version_errors, inherits, logical(1), "error")),
    paste(version_messages, collapse = " | ")
)
invalid <- mat
invalid[[1L]] <- Inf
list_error <- capture_error(
    genefunnel::genefunnel(invalid, gene_sets, BPPARAM = serial)
)
compiled_error <- capture_error(score_compiled(invalid, compiled, serial))
record(
    "matrix_error_identity",
    identical(class(compiled_error), class(list_error)) &&
        identical(error_message(compiled_error), error_message(list_error)),
    error_message(compiled_error)
)

fresh <- BiocParallel::SnowParam(workers = 2L, type = "SOCK", progressbar = FALSE)
fresh_start <- unname(proc.time()[["elapsed"]])
fresh <- BiocParallel::bpstart(fresh)
fresh_startup_elapsed <- unname(proc.time()[["elapsed"]]) - fresh_start
fresh_score_start <- unname(proc.time()[["elapsed"]])
fresh_scores <- capture_call(score_compiled(mat, compiled, fresh))
fresh_score_elapsed <- unname(proc.time()[["elapsed"]]) - fresh_score_start
fresh_paths <- BiocParallel::bplapply(
    seq_len(2L),
    function(index) normalizePath(find.package("genefunnel")),
    BPPARAM = fresh
)
BiocParallel::bpstop(fresh)
record(
    "fresh_sock_identity",
    identical(fresh_scores$value, list_call$value) &&
        identical(fresh_scores$warnings, list_call$warnings) &&
        all(vapply(fresh_paths, identical, logical(1), expected_package)),
    sprintf(
        "startup %.6f sec; score %.6f sec",
        fresh_startup_elapsed,
        fresh_score_elapsed
    )
)

reused <- BiocParallel::SnowParam(
    workers = 2L,
    type = "SOCK",
    progressbar = FALSE
)
reused_start <- unname(proc.time()[["elapsed"]])
reused <- BiocParallel::bpstart(reused)
reused_startup_elapsed <- unname(proc.time()[["elapsed"]]) - reused_start
on.exit(BiocParallel::bpstop(reused), add = TRUE)
first_reused_start <- unname(proc.time()[["elapsed"]])
first_reused <- capture_call(score_compiled(mat, compiled, reused))
first_reused_elapsed <- unname(proc.time()[["elapsed"]]) - first_reused_start
second_reused_start <- unname(proc.time()[["elapsed"]])
second_reused <- capture_call(score_compiled(
    Matrix::Matrix(mat, sparse = TRUE),
    compiled,
    reused
))
second_reused_elapsed <- unname(proc.time()[["elapsed"]]) - second_reused_start
worker_bytes <- BiocParallel::bplapply(
    list(bytes, bytes),
    function(payload) {
        object <- genefunnel:::.unserialize_compiled_catalogue(payload)
        genefunnel:::.serialize_compiled_catalogue(object)
    },
    BPPARAM = reused
)
record(
    "reused_sock_identity",
    identical(first_reused$value, list_call$value) &&
        identical(second_reused$value, list_call$value) &&
        identical(worker_bytes, list(bytes, bytes)) &&
        BiocParallel::bpisup(reused),
    sprintf(
        "startup %.6f sec; first %.6f sec; reused %.6f sec",
        reused_startup_elapsed,
        first_reused_elapsed,
        second_reused_elapsed
    )
)

results <- do.call(rbind, rows)
results$fresh_sock_startup_sec <- fresh_startup_elapsed
results$fresh_sock_score_sec <- fresh_score_elapsed
results$reused_sock_startup_sec <- reused_startup_elapsed
results$reused_sock_first_score_sec <- first_reused_elapsed
results$reused_sock_warm_score_sec <- second_reused_elapsed
utils::write.table(
    results,
    result_path,
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE,
    quote = TRUE,
    qmethod = "double",
    na = "NA"
)
if (!all(results$passed)) {
    quit(save = "no", status = 1L)
}
cat("Catalogue correctness: ", sum(results$passed), "/", nrow(results), " passed.\n", sep = "")
