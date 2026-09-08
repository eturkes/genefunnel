# GeneFunnel

- Package gate = `Rscript tools/check-package.R --bioc`; dependencies → applicable R library.
- Publication gate = `Rscript publication/reproduce.R`; source install → `.checks/publication-library`; deterministic analysis tables/PNGs replay byte-identically.
- Publication inputs = bundled checksum-pinned GEO fixtures + installed `airway`; canonical comparators = `singscore`, ssGSEA via `GSVA`, and GSVA. Preserve preprocessing/estimand differences.
- Score = finite-size Pietra-penalized total; within-sample equality ≠ coexpression, activation, causal effect, or validated biological accuracy.
- Public article = `publication/article.Rmd`; experiment protocol/code = `publication/analyze.R`; fixed review rows = `publication/ACCEPTANCE.md`.
