# GSE125583 presynaptic-release fixture

`GSE125583_presynaptic_release.csv` contains metadata and six presynaptic genes'
nRPKM values for 289 subjects from the GSE125583 fusiform-gyrus bulk RNA-seq
cohort.

Run `data-raw/build-gse125583-example.R` from the package root to regenerate the
fixture. The script downloads the GEO series matrix and processed tables,
verifies published MD5 checksums, extracts the fixed panel, and checks fixture MD5
`4f3038f701f2015f35f46c938e3b0cdd`.

The vignette scales genes to their control medians before scoring. Its panel
and summaries are illustrative, not part of the source study's analysis.
