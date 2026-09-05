.run_quickstart_example <- function() {
    script <- system.file("examples", "quickstart.R", package = "genefunnel",
        mustWork = TRUE)
    work <- tempfile("genefunnel-example-")
    dir.create(work)
    previous <- setwd(work)
    on.exit({
        setwd(previous)
        unlink(work, recursive = TRUE)
    }, add = TRUE)
    example <- new.env(parent = globalenv())
    capture.output(suppressMessages(sys.source(script, envir = example)))
    files <- list.files(example$output_dir, full.names = TRUE)
    checksums <- tools::md5sum(files)
    capture.output(suppressMessages(sys.source(script, envir = example)))
    expect_identical(tools::md5sum(files), checksums)
    expect_identical(list.files(work), "genefunnel-results")
    list(
        example = example,
        files = basename(files),
        scores = utils::read.csv(file.path(example$output_dir, "scores.csv")),
        coverage = utils::read.csv(file.path(example$output_dir,
            "gene-set-coverage.csv")),
        saved = readRDS(file.path(example$output_dir, "result.rds")),
        png = readBin(file.path(example$output_dir, "comparison.png"),
            what = "raw", n = 24L)
    )
}

test_that("the installed quickstart exports reproducible, interpretable scores", {
    skip_if_not_installed("ggplot2")
    run <- .run_quickstart_example()
    expect_setequal(run$files, c("scores.csv", "gene-set-coverage.csv",
        "result.rds", "comparison.png"))
    expect_identical(run$scores$sample, c("coherent", "split", "concentrated"))
    expect_identical(run$scores$gene_set, rep("pathway", 3L))
    expect_equal(run$scores$member_total, c(12, 12, 12))
    expect_equal(run$scores$score, c(12, 4, 0))
    expect_equal(run$scores$observed_members, rep(4L, 3L))
    expect_equal(run$scores$coverage, rep(1, 3L))
    expect_equal(run$scores, run$example$score_table)
    expect_equal(run$coverage, run$example$result$set_coverage)
    expect_identical(run$saved, run$example$result)
    expect_equal(unname(run$saved$scores["pathway", ]), c(12, 4, 0))
    expect_equal(run$example$score_plot$data$value, c(12, 12, 12, 12, 4, 0))
    expect_identical(run$png[1:8], as.raw(c(137, 80, 78, 71, 13, 10, 26, 10)))
    expect_equal(sum(as.integer(run$png[17:20]) * 256^(3:0)), 1440)
    expect_equal(sum(as.integer(run$png[21:24]) * 256^(3:0)), 864)
})
