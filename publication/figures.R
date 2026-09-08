render_figures <- function(results, output) {
    dir.create(output, recursive = TRUE, showWarnings = FALSE)
    read <- function(name) read.csv(file.path(results, paste0(name, ".csv")), check.names = FALSE)
    labels <- function(x) gsub("_", " ", x, fixed = TRUE)
    palette <- c(GeneFunnel = "#007F82", Mean = "#C85B3F", Half_penalty = "#8979A6",
        Geometric_log1p = "#C49A27", Minimum = "#777777", singscore = "#395C99",
        ssGSEA = "#944F75", GSVA = "#648348")
    shapes <- setNames(c(16, 17, 15, 18, 8, 0, 1, 2), names(palette))
    theme <- ggplot2::theme_minimal(base_size = 11) + ggplot2::theme(
        panel.grid.minor = ggplot2::element_blank(), legend.position = "bottom")
    save <- function(p, filename, width, height) {
        ggplot2::ggsave(file.path(output, filename), p, width = width, height = height,
            dpi = 300, type = "cairo")
    }
    summary <- read("simulation-summary")
    scenario_order <- c("uniform_activation", "sparse_activation", "same_total_equality", "null",
        "library_size_nuisance", "mixed_direction", "heterogeneous_activation", "correlated_activation")
    summary$scenario <- factor(summary$scenario, levels = scenario_order,
        labels = c("Uniform activation", "Sparse activation", "Same-total equality", "Null",
            "Library-size nuisance", "Mixed direction", "Heterogeneous activation", "Correlated activation"))
    summary$size <- factor(summary$size, labels = paste("n =", sort(unique(summary$size))))
    summary$method <- factor(summary$method, levels = names(palette))
    p <- ggplot2::ggplot(summary, ggplot2::aes(method, mean_auc, colour = method, shape = method)) +
        ggplot2::geom_hline(yintercept = .5, linetype = 2, colour = "grey60") +
        ggplot2::geom_point(size = 2) + ggplot2::geom_errorbar(ggplot2::aes(
            ymin = mean_auc - 1.96 * mc_se, ymax = mean_auc + 1.96 * mc_se), width = .15) +
        ggplot2::facet_wrap(ggplot2::vars(size, scenario), ncol = 4) +
        ggplot2::scale_colour_manual(values = palette, labels = labels) +
        ggplot2::scale_shape_manual(values = shapes, labels = labels) +
        ggplot2::coord_cartesian(ylim = c(0, 1)) +
        ggplot2::labs(x = NULL, y = "Mean directional AUROC ± 1.96 Monte Carlo SE", colour = "Method", shape = "Method") +
        theme + ggplot2::theme(axis.text.x = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank(),
            strip.text = ggplot2::element_text(size = 10), legend.text = ggplot2::element_text(size = 9)) +
        ggplot2::guides(colour = ggplot2::guide_legend(nrow = 2), shape = ggplot2::guide_legend(nrow = 2))
    save(p, "simulation.png", 11, 9)
    ad <- read("ad-inference"); a <- ad[ad$metric == "rho", ]
    p <- ggplot2::ggplot(a, ggplot2::aes(estimate, reorder(method, estimate), colour = method)) +
        ggplot2::geom_errorbar(ggplot2::aes(xmin = low, xmax = high), orientation = "y", width = .15) +
        ggplot2::geom_point(size = 2) + ggplot2::scale_colour_manual(values = palette) +
        ggplot2::scale_y_discrete(labels = labels) +
        ggplot2::labs(x = "Spearman correlation with Braak stage; 95% subject bootstrap interval", y = NULL) +
        theme + ggplot2::theme(legend.position = "none")
    save(p, "ad-association.png", 8, 3.4)
    pairs <- read("airway-pairs")
    pairs <- pairs[pairs$preprocessing == "CPM", ]
    p <- ggplot2::ggplot(pairs, ggplot2::aes(donor, treated_minus_control, colour = method, shape = donor)) +
        ggplot2::geom_hline(yintercept = 0, colour = "grey65", linetype = 2) +
        ggplot2::geom_point(size = 2) +
        ggplot2::facet_wrap(ggplot2::vars(method), scales = "free_y", ncol = 4,
            labeller = ggplot2::as_labeller(labels)) +
        ggplot2::scale_colour_manual(values = palette) +
        ggplot2::scale_y_continuous(labels = function(values)
            format(signif(values, 3), trim = TRUE, scientific = FALSE, digits = 3)) +
        ggplot2::labs(x = NULL, y = "Treated − control; method-specific units", shape = "Donor") +
        theme + ggplot2::guides(colour = "none") +
        ggplot2::theme(axis.text.x = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank())
    save(p, "airway-pairs.png", 10, 5)
}
