# Assisted-by: OpenAI Codex.

SENSITIVITY_PROTOCOL_VERSION <- "E-1.0.0"
SENSITIVITY_PROTOCOL_MD5 <- "09d9429ae18f64ba00b3bd84955ad71d"

sensitivity_read_registry <- function(path) {
    registry <- utils::read.delim(
        path,
        stringsAsFactors = FALSE,
        check.names = FALSE,
        colClasses = "character",
        na.strings = character()
    )
    required <- c("protocol_version", "section", "key", "value")
    if (!identical(names(registry), required) || nrow(registry) == 0L) {
        stop("Sensitivity registry schema is invalid.", call. = FALSE)
    }
    if (any(!nzchar(as.matrix(registry))) ||
        !identical(unique(registry$protocol_version),
            SENSITIVITY_PROTOCOL_VERSION) ||
        anyDuplicated(paste(registry$section, registry$key, sep = "::"))) {
        stop("Sensitivity registry contents are invalid.", call. = FALSE)
    }
    registry
}

sensitivity_registry_value <- function(registry, section, key) {
    selected <- registry$section == section & registry$key == key
    if (sum(selected) != 1L) {
        stop(
            "Sensitivity registry key is missing or duplicated: ",
            section,
            "::",
            key,
            call. = FALSE
        )
    }
    registry$value[selected]
}

sensitivity_registry_levels <- function(registry, key) {
    value <- sensitivity_registry_value(registry, "controlled", key)
    strsplit(value, ",", fixed = TRUE)[[1L]]
}

sensitivity_scenario_keys <- function() {
    c(
        "member_count", "archetype", "dynamic_range", "log_sd",
        "expected_total", "dispersion", "dropout_fraction",
        "profile_replicate"
    )
}

sensitivity_controlled_design <- function(registry) {
    keys <- sensitivity_scenario_keys()
    levels <- lapply(keys, sensitivity_registry_levels, registry = registry)
    design <- do.call(expand.grid, c(
        levels,
        list(KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
    ))
    names(design) <- keys
    design$scenario_id <- seq_len(nrow(design))
    design$fold <- as.integer(design$profile_replicate)
    design
}

sensitivity_mask_design <- function(registry) {
    masks <- expand.grid(
        removed_fraction = as.numeric(sensitivity_registry_levels(
            registry, "removed_fraction"
        )),
        mask_mechanism = sensitivity_registry_levels(
            registry, "mask_mechanism"
        ),
        mask_repeat = seq_len(as.integer(sensitivity_registry_value(
            registry, "controlled", "mask_repeats"
        ))),
        KEEP.OUT.ATTRS = FALSE,
        stringsAsFactors = FALSE
    )
    masks$fraction_index <- match(
        masks$removed_fraction,
        unique(masks$removed_fraction)
    )
    masks$mechanism_index <- match(
        masks$mask_mechanism,
        unique(masks$mask_mechanism)
    )
    masks
}

sensitivity_design_counts <- function(registry) {
    scenario_keys <- sensitivity_scenario_keys()
    mask_keys <- c(
        "removed_fraction", "mask_mechanism", "mask_repeats",
        "absence_mode"
    )
    level_count <- function(key) {
        if (key == "mask_repeats") {
            return(as.integer(sensitivity_registry_value(
                registry, "controlled", key
            )))
        }
        length(sensitivity_registry_levels(registry, key))
    }
    scenarios <- as.integer(prod(vapply(
        scenario_keys, level_count, integer(1L)
    )))
    feature_rows <- as.integer(scenarios * prod(vapply(
        mask_keys, level_count, integer(1L)
    )))
    c(scenarios = scenarios, feature_rows = feature_rows)
}

sensitivity_validate_design <- function(registry, design) {
    masks <- sensitivity_mask_design(registry)
    member_counts <- as.integer(sensitivity_registry_levels(
        registry, "member_count"
    ))
    remaining <- outer(
        member_counts,
        masks$removed_fraction,
        function(size, fraction) {
            size - pmax(1, pmin(size - 3, floor(fraction * size + 0.5)))
        }
    )
    seed_base <- as.integer(sensitivity_registry_value(
        registry, "controlled", "mask_seed_base"
    ))
    seeds <- unlist(lapply(design$scenario_id, function(scenario_id) {
        seed_base + 1000L * (scenario_id - 1L) +
            100L * (masks$fraction_index - 1L) +
            10L * (masks$mechanism_index - 1L) + masks$mask_repeat - 1L
    }), use.names = FALSE)
    stopifnot(
        all(remaining >= 3L),
        !anyDuplicated(seeds),
        max(seeds) <= .Machine$integer.max,
        identical(as.integer(table(design$fold)), rep.int(576L, 10L))
    )
    invisible(masks)
}

sensitivity_validate_protocol <- function(root = ".") {
    registry_path <- file.path(root, "benchmark", "sensitivity-protocol.tsv")
    observed_md5 <- unname(tools::md5sum(registry_path))
    if (!identical(observed_md5, SENSITIVITY_PROTOCOL_MD5)) {
        stop("Sensitivity protocol identity is invalid.", call. = FALSE)
    }
    registry <- sensitivity_read_registry(registry_path)
    counts <- sensitivity_design_counts(registry)
    design <- sensitivity_controlled_design(registry)
    expected_scenarios <- as.integer(sensitivity_registry_value(
        registry, "controlled", "scenario_count"
    ))
    expected_feature <- as.integer(sensitivity_registry_value(
        registry, "controlled", "feature_rows"
    ))
    expected_technical <- as.integer(sensitivity_registry_value(
        registry, "controlled", "technical_rows"
    ))
    folds <- as.integer(sensitivity_registry_levels(
        registry, "profile_replicate"
    ))
    stopifnot(
        nrow(registry) == 90L,
        nrow(design) == expected_scenarios,
        identical(unname(counts[["scenarios"]]), expected_scenarios),
        identical(unname(counts[["feature_rows"]]), expected_feature),
        identical(expected_technical, expected_scenarios),
        identical(folds, seq_len(10L)),
        expected_scenarios %% length(folds) == 0L
    )
    sensitivity_validate_design(registry, design)
    list(
        registry_rows = nrow(registry),
        scenarios = expected_scenarios,
        feature_rows = expected_feature,
        technical_rows = expected_technical,
        scenarios_per_fold = expected_scenarios %/% length(folds)
    )
}
