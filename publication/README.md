# GeneFunnel article and evidence

This directory contains a separate software-method article and its reproducible analyses.
The article evaluates a Pietra-penalized total. It does not claim a new inequality index or universal biological accuracy.

## Reproduce

Install the package dependencies from `DESCRIPTION`.
Also install `airway`, `singscore`, and `GSVA` from Bioconductor.
Install `rmarkdown`, `knitr`, and Pandoc for article rendering.
Install Python 3 for the offline HGNC mapping replay.
GSVA dependencies can require the ImageMagick C++ development libraries on Linux.

Run this command from the repository root:

```sh
Rscript publication/reproduce.R
```

The script installs the current source into `.checks/publication-library`.
It generates the evidence, repeats the deterministic analysis, compares bytes, and renders HTML and DOCX articles.
It uses the bundled, checksum-pinned public-data fixtures and the installed `airway` data.
It does not download data during analysis.

Use `--no-replay` only for exploratory iterations.
Use `--no-render` to run the evidence gate without Pandoc.
After a full run, use `--replay-only` to verify existing outputs in a fresh directory.
The replay regression rejects stale orphan files.
Run the package gate separately:

```sh
Rscript tools/check-package.R --bioc
```

## Evidence

- `article.Rmd`, `article.html`, `article.docx`: editable source and readable article formats.
- `analyze.R`: complete experiment definitions, seeds, estimators, and figure generation.
- `results/`: subject-level outputs, simulation replicates, uncertainty, coverage, runtime, and environment.
- `figures/`: standalone publication figures.
- `ACCEPTANCE.md`, `REVIEW.md`: fixed acceptance criteria and adjudications.

Timing values describe this machine. They are excluded from byte-level replay.
Object sizes describe stored inputs, not peak process memory.
The article's declarations must be confirmed before submission.
