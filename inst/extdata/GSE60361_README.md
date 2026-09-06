# GSE60361 mouse brain cell-type fixture

`GSE60361_brain_cell_types.rds` aggregates Zeisel et al.'s mouse cortex and
hippocampus single-cell RNA-seq data. It contains molecule counts for 19,972
endogenous genes across 47 subtypes and seven broad classes, representing
2,816 cells with named subtypes.

Run `data-raw/build-gse60361-example.R` from the package root to regenerate the
fixture. The script retrieves the `zeisel-brain-2015` dataset through
`scRNAseq::ZeiselBrainData()`, validates the source dimensions and annotations,
aggregates counts, and checks the generated fixture against MD5
`93e3022388eefd47dc06dfcad9c194aa`. The scRNAseq resource is versioned
`2023-12-14`.

The vignette selects genes by median CPM across the seven broad classes,
without subtype labels. The top 500 and 5,000 form nested abundance panels,
demonstrating larger-set scoring rather than pathways or cell-type signatures.
