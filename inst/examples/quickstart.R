# Run from a writable folder. Re-running replaces these example output files.
# All inputs are synthetic. No data are downloaded.
if (!requireNamespace("genefunnel", quietly = TRUE)) {
    stop("Install genefunnel before running this example.", call. = FALSE)
}
if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Install ggplot2 to draw the example figure.", call. = FALSE)
}

x <- rbind(
    A = c(coherent = 3, split = 6, concentrated = 12),
    B = c(coherent = 3, split = 6, concentrated = 0),
    C = c(coherent = 3, split = 0, concentrated = 0),
    D = c(coherent = 3, split = 0, concentrated = 0)
)
gene_sets <- list(pathway = rownames(x))
result <- genefunnel::genefunnel(
    x, gene_sets,
    BPPARAM = BiocParallel::SerialParam(),
    details = TRUE
)
score_table <- data.frame(
    gene_set = "pathway",
    sample = colnames(x),
    member_total = unname(colSums(x)),
    score = unname(result$scores["pathway", ]),
    observed_members = unname(result$observed["pathway", ]),
    coverage = unname(result$coverage["pathway", ])
)

comparison <- data.frame(
    sample = factor(rep(score_table$sample, 2L), levels = colnames(x)),
    measure = factor(rep(c("Member total", "GeneFunnel"), each = ncol(x)),
        levels = c("Member total", "GeneFunnel")),
    value = c(score_table$member_total, score_table$score)
)
score_plot <- ggplot2::ggplot(
    comparison, ggplot2::aes(x = sample, y = value, fill = measure)
) +
    ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.7),
        width = 0.6) +
    ggplot2::geom_text(ggplot2::aes(label = value),
        position = ggplot2::position_dodge(width = 0.7), vjust = -0.5) +
    ggplot2::scale_fill_manual(values = c(
        "Member total" = "#A6ADB5", "GeneFunnel" = "#0F5E62")) +
    ggplot2::scale_x_discrete(labels = c(
        coherent = "Coherent\n(3, 3, 3, 3)",
        split = "Split\n(6, 6, 0, 0)",
        concentrated = "Concentrated\n(12, 0, 0, 0)")) +
    ggplot2::scale_y_continuous(limits = c(0, 14), breaks = seq(0, 12, 4),
        expand = ggplot2::expansion(mult = c(0, 0))) +
    ggplot2::labs(
        title = "Same total, different GeneFunnel scores",
        subtitle = "One four-gene set; all members observed in every sample",
        x = NULL, y = "Total or score (input units)", fill = NULL,
        caption = paste(
            "Synthetic inputs. GeneFunnel penalizes disagreement within the set.",
            "A zero score is not a claim of biological inactivity.", sep = "\n")
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
        panel.grid.major.x = ggplot2::element_blank(),
        panel.grid.minor = ggplot2::element_blank(),
        legend.position = "top",
        plot.title.position = "plot",
        plot.caption = ggplot2::element_text(hjust = 0)
    )

output_dir <- "genefunnel-results"
if (!dir.exists(output_dir) && !dir.create(output_dir)) {
    stop("Cannot create genefunnel-results in the working folder.", call. = FALSE)
}
utils::write.csv(score_table, file.path(output_dir, "scores.csv"),
    row.names = FALSE)
utils::write.csv(result$set_coverage,
    file.path(output_dir, "gene-set-coverage.csv"), row.names = FALSE)
saveRDS(result, file.path(output_dir, "result.rds"))
ggplot2::ggsave(file.path(output_dir, "comparison.png"), score_plot,
    width = 8, height = 4.8, units = "in", dpi = 180, bg = "white")
print(score_table, row.names = FALSE)
message("Example files: ", normalizePath(output_dir))
if (interactive()) {
    print(score_plot)
}
