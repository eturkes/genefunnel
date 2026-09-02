# GSE125583 presynaptic-release fixture

`GSE125583_presynaptic_release.csv` is a compact extract from the independent
fusiform-gyrus bulk RNA-seq cohort in GEO accession GSE125583. It contains 289
subjects, selected sample metadata, and processed nRPKM values for six
presynaptic genes used by the vignette.

Run `data-raw/build-gse125583-example.R` from the package root to regenerate the
fixture. The script downloads the GEO series matrix and per-sample processed
tables, verifies their published byte-level MD5 checksums, extracts the fixed
gene panel, and checks the generated fixture against MD5
`4f3038f701f2015f35f46c938e3b0cdd`.

The vignette scales each gene to its median among control subjects before
scoring. The selected panel and downstream summaries are illustrative and are
not part of the source study's reported analysis.
