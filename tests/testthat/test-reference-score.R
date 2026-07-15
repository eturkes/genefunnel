# Assisted-by: OpenAI Codex.

test_that("the reference scorer implements the normative formula", {
    expect_identical(reference_score(c(4, 4, 4)), 12)
    expect_identical(reference_score(c(0, 0, 0)), 0)
    expect_identical(reference_score(c(4, 0)), 0)
    expect_identical(reference_score(c(4, 0, 0)), 0)
    expect_identical(reference_score(c(1, 2, NA_real_)), 2)
    expect_true(is.na(reference_score(c(4, NA_real_))))
    expect_identical(reference_score(c(1e150, 1, 0)), 1.5)
})
