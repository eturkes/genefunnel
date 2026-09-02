# GSE129308 vignette fixture

`GSE129308_AD_pseudobulk_counts.csv` contains selected-gene UMI totals from the
16 paired Alzheimer’s disease samples used in the GeneFunnel vignette. The
samples comprise MAP2-sorted NFT-free neurons and AT8-sorted NFT-bearing neurons
from eight donors.

The fixture was generated from the Cell Ranger-filtered HDF5 matrices in the GEO
archive `GSE129308_RAW.tar`. Counts were summed across deposited barcodes for
each donor-population sample. The fixture also records barcode counts and total
library sizes. It does not reproduce the thesis quality-control pipeline or its
single-soma-first scoring strategy.

Source: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE129308

Archive SHA-256:
`dcd8ffee06d71c50f30e04b1aa7a6658dfc9a359e3a53b4e59f0f8afeabe9e0a`

Regenerate the file from the repository root:

```bash
python data-raw/build-gse129308-example.py
```

The regeneration script requires Python, NumPy, and h5py.
