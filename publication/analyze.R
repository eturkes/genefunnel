# Deterministic evidence; invoked by reproduce.R after a clean local install.
args <- commandArgs(trailingOnly = TRUE)
out <- if (length(args)) args[[1L]] else "publication/results"
fig <- if (length(args) > 1L) args[[2L]] else "publication/figures"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
dir.create(fig, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(".checks/publication-library", "../.scratch/publication/R-library",
    "../.R-library", .libPaths()))
stopifnot(requireNamespace("genefunnel"), requireNamespace("ggplot2"))
options(stringsAsFactors = FALSE, digits = 16)
Sys.setenv(TZ = "UTC")
RNGkind("Mersenne-Twister", "Inversion", "Rejection")
set.seed(7447)
serial <- BiocParallel::SerialParam(progressbar = FALSE)
write_table <- function(x, name) {
    utils::write.csv(x, file.path(out, paste0(name, ".csv")), row.names = FALSE,
        na = "NA")
}
score <- function(x, sets = list(target = rownames(x))) {
    genefunnel::genefunnel(x, sets, BPPARAM = serial)
}
auc <- function(x, y) {
    keep <- is.finite(x) & !is.na(y)
    x <- signif(x[keep], 12L); y <- y[keep]
    n1 <- sum(y); n0 <- sum(!y)
    if (!n1 || !n0) return(NA_real_)
    (sum(rank(x, ties.method = "average")[y]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}
rho <- function(x, y) suppressWarnings(cor(x, y, method = "spearman", use = "complete.obs"))
interval <- function(x) unname(quantile(x, c(.025, .975), na.rm = TRUE))
core_scores <- function(x, members = rownames(x)) {
    z <- x[members, , drop = FALSE]
    n <- colSums(!is.na(z))
    gf <- as.numeric(score(x, list(target = members))) / n
    means <- colMeans(z, na.rm = TRUE)
    rbind(GeneFunnel = gf, Mean = means,
        Half_penalty = (gf + means) / 2,
        Geometric_log1p = expm1(colMeans(log1p(z), na.rm = TRUE)),
        Minimum = apply(z, 2L, min, na.rm = TRUE))
}
canonical <- function(x, members) {
    stopifnot(requireNamespace("singscore", quietly = TRUE),
        requireNamespace("GSVA", quietly = TRUE))
    ranks <- singscore::rankGenes(x, tiesMethod = "average")
    sing <- singscore::simpleScore(ranks, upSet = members, centerScore = TRUE,
        knownDirection = TRUE)$TotalScore
    par <- GSVA::ssgseaParam(log1p(x), list(target = members),
        minSize = 2, maxSize = Inf, alpha = .25, normalize = FALSE,
        checkNA = "no", verbose = FALSE, ondisk = "no")
    ss <- suppressMessages(GSVA::gsva(par, verbose = FALSE, BPPARAM = serial))
    rbind(singscore = sing, ssGSEA = as.numeric(ss))
}
gsva_score <- function(x, members) {
    par <- GSVA::gsvaParam(log1p(x), list(target = members),
        minSize = 2, maxSize = Inf, kcdf = "Gaussian",
        maxDiff = TRUE, absRanking = FALSE, checkNA = "no", verbose = FALSE, ondisk = "no")
    as.numeric(suppressMessages(GSVA::gsva(par, verbose = FALSE, BPPARAM = serial)))
}

# Independent algebra oracle; integer grid contains ties, zeros, and plateaus.
algebra <- list(); a <- 0L
for (n in 2:5) {
    grid <- as.matrix(expand.grid(rep(list(0:3), n)))
    x <- t(grid); rownames(x) <- paste0("g", seq_len(n))
    actual <- as.numeric(score(x))
    total <- colSums(x); center <- colMeans(x)
    deviation <- colSums(abs(sweep(x, 2L, center)))
    pietra <- ifelse(total == 0, 0, deviation / (2 * total))
    oracle <- total * (1 - n * pietra / (n - 1))
    inc <- x; inc[1L, ] <- inc[1L, ] + 1
    perm <- x[n:1L, , drop = FALSE]
    mid <- (x + perm) / 2
    checks <- c(identity = max(abs(actual - oracle)),
        lower_bound = max(pmax(-actual, 0)),
        upper_bound = max(pmax(actual - total, 0)),
        monotonicity = max(pmax(actual - as.numeric(score(inc)), 0)),
        concavity = max(pmax(actual - as.numeric(score(mid)), 0)),
        homogeneity = max(abs(as.numeric(score(x * 3)) - 3 * actual)),
        translation = max(abs(as.numeric(score(x + 2)) - actual - 2 * n)))
    a <- a + 1L
    algebra[[a]] <- data.frame(n = n, vectors = ncol(x), check = names(checks),
        maximum_error = unname(checks))
}
algebra <- do.call(rbind, algebra)
stopifnot(all(algebra$maximum_error < 1e-10))
write_table(algebra, "algebra")

# Protocol: two set sizes × eight mechanisms × 50 independent replicates.
# Same-total equality is a construction probe, not biological accuracy.
scenarios <- c("uniform_activation", "sparse_activation", "same_total_equality",
    "null", "library_size_nuisance", "mixed_direction",
    "heterogeneous_activation", "correlated_activation")
sim <- list(); s <- 0L
for (size in c(10L, 50L)) for (scenario in scenarios) {
    message("Simulation: n=", size, "; ", scenario)
    matrices <- vector("list", 50L)
    for (replicate in seq_len(50L)) {
        x <- matrix(rgamma(500L * 80L, shape = 2, rate = 2), 500L, 80L)
        members <- seq_len(size); cases <- 41:80
        if (scenario == "uniform_activation") x[members, cases] <- 1.5 * x[members, cases]
        if (scenario == "sparse_activation") {
            active <- seq_len(max(1L, as.integer(size * .1)))
            x[active, cases] <- 6 * x[active, cases]
        }
        if (scenario == "heterogeneous_activation") {
            x <- x * rlnorm(nrow(x), meanlog = 0, sdlog = 1.5)
            x[members, cases] <- 1.5 * x[members, cases]
        }
        if (scenario == "correlated_activation") {
            x[members, ] <- sweep(x[members, , drop = FALSE], 2L,
                rlnorm(ncol(x), meanlog = -.08, sdlog = .4), "*")
            x[members, cases] <- 1.5 * x[members, cases]
        }
        if (scenario == "same_total_equality") {
            x[members, 1:40] <- matrix(rgamma(size * 40L, .5), size)
            x[members, cases] <- matrix(rgamma(size * 40L, 8), size)
            z <- x[members, , drop = FALSE]
            x[members, ] <- sweep(z, 2L, colSums(z), "/") * size
        }
        if (scenario == "library_size_nuisance") x[, cases] <- 1.5 * x[, cases]
        if (scenario == "mixed_direction") {
            x[seq_len(size / 2), cases] <- 1.5 * x[seq_len(size / 2), cases]
            x[(size / 2 + 1):size, cases] <- x[(size / 2 + 1):size, cases] / 1.5
        }
        rownames(x) <- paste0("g", seq_len(nrow(x)))
        matrices[[replicate]] <- x
    }
    x <- do.call(cbind, matrices); colnames(x) <- paste0("s", seq_len(ncol(x)))
    members <- rownames(x)[seq_len(size)]
    z <- rbind(core_scores(x, members), canonical(x, members))
    for (replicate in seq_len(50L)) {
        idx <- ((replicate - 1L) * 80L + 1L):(replicate * 80L)
        block <- rbind(z[, idx, drop = FALSE],
            GSVA = gsva_score(x[, idx, drop = FALSE], members))
        vals <- apply(block, 1L, auc, y = rep(c(FALSE, TRUE), each = 40))
        s <- s + 1L
        sim[[s]] <- data.frame(size = size, scenario = scenario,
            replicate = replicate, method = names(vals), auc = unname(vals))
    }
}
sim <- do.call(rbind, sim); write_table(sim, "simulation-replicates")
keys <- unique(sim[c("size", "scenario", "method")])
summary <- do.call(rbind, lapply(seq_len(nrow(keys)), function(i) {
    k <- keys[i, ]; z <- sim[sim$size == k$size & sim$scenario == k$scenario, ]
    m <- z$auc[z$method == k$method]; gf <- z$auc[z$method == "GeneFunnel"]
    delta <- gf - m
    cbind(k, replicates = length(m), mean_auc = mean(m),
        mc_se = sd(m) / sqrt(length(m)),
        delta_gf = mean(delta),
        delta_low = mean(delta) - qt(.975, length(m) - 1) * sd(delta) / sqrt(length(m)),
        delta_high = mean(delta) + qt(.975, length(m) - 1) * sd(delta) / sqrt(length(m)))
}))
write_table(summary, "simulation-summary")

# Missingness stress; same complete profiles and exact masked counts across mechanisms.
x <- matrix(rgamma(20 * 300, 2, 2), 20, dimnames = list(paste0("g", 1:20), NULL))
complete <- as.numeric(score(x)); coverage <- list(); k <- 0L
for (fraction in c(.1, .3, .5)) for (mechanism in c("random", "lowest", "highest")) {
    z <- x
    for (j in seq_len(ncol(z))) {
        idx <- switch(mechanism, random = sample.int(20, 20 * fraction),
            lowest = head(order(x[, j]), 20 * fraction),
            highest = head(order(x[, j], decreasing = TRUE), 20 * fraction))
        z[idx, j] <- NA_real_
    }
    actual <- as.numeric(score(z)); normalized <- actual * 20 / colSums(!is.na(z))
    for (policy in c(0, .8, 1)) {
        accepted <- colMeans(!is.na(z)) >= policy
        k <- k + 1L
        coverage[[k]] <- data.frame(missing_fraction = fraction, mechanism = mechanism,
            min_coverage = policy, retained_profiles = sum(accepted), profiles = ncol(z),
            median_raw_relative_error = if (any(accepted)) median((actual / complete - 1)[accepted]) else NA,
            median_per_gene_relative_error = if (any(accepted)) median((normalized / complete - 1)[accepted]) else NA,
            rank_correlation = if (any(accepted)) rho(actual[accepted], complete[accepted]) else NA)
    }
}
write_table(do.call(rbind, coverage), "coverage")

# Fixture integrity and subject-level inference; resample preprocessing with subjects.
inputs <- c("inst/extdata/GSE125583_presynaptic_release.csv", "inst/extdata/GSE60361_brain_cell_types.rds")
expected <- c("4f3038f701f2015f35f46c938e3b0cdd", "93e3022388eefd47dc06dfcad9c194aa")
stopifnot(identical(unname(tools::md5sum(inputs)), expected))
write_table(data.frame(file = inputs, md5 = expected), "input-manifest")
ad <- read.csv(inputs[[1L]], check.names = FALSE)
genes <- c("SYN1", "SYN2", "SYP", "SNAP25", "STX1A", "VAMP2")
raw <- t(as.matrix(ad[genes])); control <- ad$diagnosis == "control"
stopifnot(!anyDuplicated(ad$case_id), !anyDuplicated(ad$sample))
ad_scores <- function(idx, transform = "control_median") {
    x <- raw[, idx, drop = FALSE]
    if (transform == "control_median") {
        ref <- apply(x[, control[idx], drop = FALSE], 1L, median)
        stopifnot(all(ref > 0)); x <- x / ref
    }
    if (transform == "log1p") x <- log1p(x)
    core_scores(x)
}
point <- ad_scores(seq_len(nrow(ad)))
partial_rank <- function(x, y, age, sex) {
    keep <- complete.cases(x, y, age, sex)
    d <- data.frame(x = rank(x[keep]), y = rank(y[keep]),
        age = rank(age[keep]), sex = factor(sex[keep]))
    cor(resid(lm(x ~ age + sex, d)), resid(lm(y ~ age + sex, d)))
}
ad_metrics <- function(idx, transform = "control_median") {
    z <- ad_scores(idx, transform)
    rbind(rho = apply(z, 1L, rho, y = ad$braak[idx]),
        partial_rank = apply(z, 1L, partial_rank, y = ad$braak[idx],
            age = ad$age[idx], sex = ad$sex[idx]),
        apparent_auc = apply(-z, 1L, auc, y = !control[idx]))
}
point_metrics <- ad_metrics(seq_len(nrow(ad)))
boot <- array(NA_real_, c(3L, nrow(point), 2000L), dimnames = list(
    rownames(point_metrics), colnames(point_metrics), NULL))
for (b in seq_len(dim(boot)[[3L]])) {
    idx <- c(sample(which(control), sum(control), replace = TRUE),
        sample(which(!control), sum(!control), replace = TRUE))
    boot[, , b] <- ad_metrics(idx)
}
ad_results <- do.call(rbind, lapply(rownames(point_metrics), function(metric) {
    do.call(rbind, lapply(colnames(point_metrics), function(method) {
        ci <- interval(boot[metric, method, ])
        delta <- boot[metric, "GeneFunnel", ] - boot[metric, method, ]
        delta_ci <- interval(delta)
        data.frame(metric = metric, method = method, estimate = point_metrics[metric, method],
            low = ci[[1]], high = ci[[2]], delta_gf = point_metrics[metric, "GeneFunnel"] - point_metrics[metric, method],
            delta_low = delta_ci[[1]], delta_high = delta_ci[[2]], bootstrap_replicates = 2000L)
    }))
}))
write_table(ad_results, "ad-inference")
write_table(data.frame(ad[c("sample", "case_id", "diagnosis", "age", "sex", "braak")],
    t(point), check.names = FALSE), "ad-scores")
sensitivity <- do.call(rbind, lapply(c("control_median", "raw", "log1p"), function(preprocessing) {
    m <- ad_metrics(seq_len(nrow(ad)), preprocessing)
    data.frame(preprocessing = preprocessing, method = colnames(m), rho = m["rho", ],
        partial_rank = m["partial_rank", ], apparent_auc = m["apparent_auc", ])
}))
write_table(sensitivity, "ad-preprocessing")

# Large public matrix → descriptive stress, no subtype-as-donor inference.
mouse <- readRDS(inputs[[2L]])
counts <- mouse$subtype_counts
cpm <- sweep(counts, 2, colSums(counts), "/") * 1e6
broad <- sweep(mouse$broad_class_counts, 2, colSums(mouse$broad_class_counts), "/") * 1e6
reference <- apply(broad, 1L, median)
ordered <- order(-reference, names(reference)); ranked <- names(reference)[ordered][reference[ordered] > 0]
mouse_results <- list(); m <- 0L
for (size in c(500, 5000)) for (preprocessing in c("CPM", "reference_scaled")) {
    x <- if (preprocessing == "CPM") cpm else cpm / pmax(reference, .Machine$double.eps)
    members <- head(ranked, size)
    gf <- as.numeric(score(x, list(panel = members)))
    total <- colSums(x[members, , drop = FALSE])
    m <- m + 1L
    mouse_results[[m]] <- data.frame(panel_size = size, preprocessing = preprocessing,
        mouse$subtype_metadata, total = total, GeneFunnel = gf, retained_fraction = gf / total)
}
write_table(do.call(rbind, mouse_results), "mouse-scale")

# Independent intervention labels and an external, complete positive signature.
data("airway", package = "airway", envir = environment())
aw_counts <- as.matrix(SummarizedExperiment::assay(airway))
aw_meta <- as.data.frame(SummarizedExperiment::colData(airway))
aw <- sweep(aw_counts, 2L, colSums(aw_counts), "/") * 1e6
# Shared expressed universe; constant/all-zero rows supply no between-sample information.
variable <- apply(aw, 1L, function(v) diff(range(v)) > 0)
aw <- aw[variable, , drop = FALSE]
panel <- read.delim("publication/inputs/wang2004-activated-hgnc.tsv", check.names = FALSE)
stopifnot(nrow(panel) == 50L, !anyDuplicated(panel$symbol),
    !anyDuplicated(panel$ensembl_gene_id), all(nzchar(panel$ensembl_gene_id)))
panel$present_in_assay <- panel$ensembl_gene_id %in% rownames(aw_counts)
panel$in_scoring_universe <- panel$ensembl_gene_id %in% rownames(aw)
write_table(panel, "airway-coverage")
aw_members <- panel$ensembl_gene_id[panel$in_scoring_universe]
stopifnot(length(aw_members) >= 2L)
aw_scores <- list(); aw_pairs <- list(); aw_summary <- list(); h <- 0L
for (preprocessing in c("CPM", "log1p_CPM")) {
    x <- if (preprocessing == "CPM") aw else log1p(aw)
    values <- rbind(core_scores(x, aw_members), canonical(x, aw_members),
        GSVA = gsva_score(aw, aw_members))
    aw_scores[[preprocessing]] <- data.frame(preprocessing = preprocessing,
        sample = colnames(x), donor = aw_meta$cell, treatment = aw_meta$dex,
        t(values), check.names = FALSE)
    for (method in rownames(values)) {
        donors <- unique(as.character(aw_meta$cell))
        differences <- vapply(donors, function(donor) {
            selected <- as.character(aw_meta$cell) == donor
            unname(values[method, selected & aw_meta$dex == "trt"] -
                values[method, selected & aw_meta$dex == "untrt"])
        }, numeric(1))
        nonzero <- signif(differences, 12) != 0
        test <- if (any(nonzero)) {
            binom.test(sum(differences[nonzero] > 0), sum(nonzero), p = .5)
        } else list(conf.int = c(NA_real_, NA_real_), p.value = NA_real_)
        h <- h + 1L
        aw_pairs[[h]] <- data.frame(preprocessing = preprocessing, method = method,
            donor = donors, treated_minus_control = differences)
        aw_summary[[h]] <- data.frame(preprocessing = preprocessing, method = method,
            declared_genes = nrow(panel), observed_genes = length(aw_members),
            donor_pairs = length(differences), increased = sum(differences > 0),
            tied = sum(!nonzero), median_difference = median(differences),
            direction_fraction = if (any(nonzero)) mean(differences[nonzero] > 0) else NA_real_,
            fraction_low = test$conf.int[[1]], fraction_high = test$conf.int[[2]],
            two_sided_sign_p = test$p.value)
    }
}
write_table(do.call(rbind, aw_scores), "airway-scores")
write_table(do.call(rbind, aw_pairs), "airway-pairs")
write_table(do.call(rbind, aw_summary), "airway-summary")
manifest <- read.csv(file.path(out, "input-manifest.csv"))
extra_inputs <- c(system.file("data", "airway.RData", package = "airway"),
    "publication/inputs/wang2004-activated-hgnc.tsv", "publication/inputs/hgnc_complete_set.txt.gz")
manifest <- rbind(manifest, data.frame(file = c("airway/data/airway.RData", extra_inputs[-1L]),
    md5 = unname(tools::md5sum(extra_inputs))))
write_table(manifest, "input-manifest")

# Timing is an observation, not part of deterministic replay or a superiority claim.
if (!identical(Sys.getenv("GENEFUNNEL_SKIP_TIMING"), "1")) {
    set.seed(9101L)
    timing <- list(); t <- 0L
    for (density in c(.05, .5, 1)) {
        x <- matrix(rexp(5000 * 200), 5000, dimnames = list(paste0("g", 1:5000), NULL))
        x[matrix(runif(length(x)) > density, nrow(x))] <- 0
        sets <- setNames(replicate(200, sample(rownames(x), 50), simplify = FALSE), paste0("set", 1:200))
        sparse <- Matrix::Matrix(x, sparse = TRUE)
        d <- score(x, sets); sp <- score(sparse, sets)
        stopifnot(isTRUE(all.equal(d, sp, tolerance = 1e-12)))
        for (storage in c("dense", "sparse")) {
            input <- if (storage == "dense") x else sparse
            invisible(score(input, sets))
            times <- replicate(7, system.time(score(input, sets))[["elapsed"]])
            t <- t + 1L
            timing[[t]] <- data.frame(density = density, storage = storage,
                features = nrow(x), samples = ncol(x), sets = length(sets), members = 50,
                repetitions = length(times), median_seconds = median(times),
                minimum_seconds = min(times), maximum_seconds = max(times),
                input_bytes = as.numeric(object.size(input)), max_absolute_difference = max(abs(d - sp)))
        }
    }
    write_table(do.call(rbind, timing), "runtime")
    hardware <- c(capture.output(Sys.info()),
        paste(names(Sys.getenv(c("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"))),
            Sys.getenv(c("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS")), sep = "="))
    if (nzchar(Sys.which("lscpu"))) hardware <- c(hardware, system2("lscpu", stdout = TRUE))
    writeLines(sub("[[:blank:]]+$", "", hardware),
        file.path(out, "runtime-environment.txt"))
}

source("publication/figures.R")
render_figures(out, fig)
packages <- sort(loadedNamespaces())
write_table(data.frame(package = packages, version = vapply(packages,
    function(p) as.character(utils::packageVersion(p)), character(1))), "dependencies")
writeLines(sub("[[:blank:]]+$", "", capture.output(sessionInfo())),
    file.path(out, "environment.txt"))
message("Deterministic analysis complete: ", out)
