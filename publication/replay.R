new_replay_directory <- function(work) {
    path <- tempfile("replay-", tmpdir = work)
    stopifnot(dir.create(path, recursive = FALSE))
    path
}

compare_replay <- function(reference, replay) {
    paths <- c(list.files(file.path(reference, "results"), pattern = "[.]csv$", full.names = TRUE),
        list.files(file.path(reference, "figures"), pattern = "[.]png$", full.names = TRUE))
    paths <- paths[!basename(paths) %in% c("runtime.csv", "replay.csv")]
    stopifnot(length(paths) > 0L)
    relative <- substring(paths, nchar(reference) + 2L)
    targets <- file.path(replay, relative)
    identical_bytes <- file.exists(targets) &
        unname(tools::md5sum(paths)) == unname(tools::md5sum(targets))
    data.frame(file = paths, md5 = unname(tools::md5sum(paths)),
        replay_identical = identical_bytes)
}
