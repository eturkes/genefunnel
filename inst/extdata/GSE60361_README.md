# GSE60361 mouse brain cell-type fixture

`GSE60361_brain_cell_types.rds` is a compact aggregation of the Zeisel et al.
mouse cortex and hippocampus single-cell RNA-seq experiment. The fixture
contains molecule counts for 19,972 endogenous genes across 47 named subtypes
and seven broad cell classes. These profiles aggregate 2,816 cells with a named
subtype annotation.

Run `data-raw/build-gse60361-example.R` from the package root to regenerate the
fixture. The script retrieves the `zeisel-brain-2015` dataset through
`scRNAseq::ZeiselBrainData()`, validates the source dimensions and annotations,
aggregates counts, and checks the generated fixture against MD5
`93e3022388eefd47dc06dfcad9c194aa`. The scRNAseq resource is versioned
`2023-12-14`.

The vignette converts the aggregated counts to CPM. It selects genes by their
median CPM across the seven broad classes, without using subtype labels. The
top 5,000 genes form the larger scale panel, and the top 500 form the nested
smaller panel. These abundance-defined panels demonstrate score behavior at
larger set sizes; they are not biological pathways or cell-type signatures.
