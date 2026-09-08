# External glucocorticoid signature

The panel contains every activated gene in Wang et al.'s Table 1.
The selection is independent of the airway expression data.

- Source article: [Wang et al., PNAS](https://doi.org/10.1073/pnas.0407008101).
- Original labels: `wang2004-activated-original.tsv`.
- Approved symbols and Ensembl identifiers: `wang2004-activated-hgnc.tsv`.
- Mapping rules and checks: `map-wang2004.py`.
- Complete frozen HGNC source: `hgnc_complete_set.txt.gz`.
- Source and output SHA-256 identities: `hgnc-provenance.json`.

Run the offline mapping recipe with Python 3:

```sh
python3 publication/inputs/map-wang2004.py
```

The recipe verifies the decompressed source SHA-256 before mapping.
It writes the same mapping and provenance bytes on repeated runs.
HGNC releases these data under [CC0](https://www.genenames.org/about/license/).
The source article's duplicate GenBank accession does not determine the mapping.
The original compound label `Kip2/p57` requires both aliases to resolve its target uniquely.
