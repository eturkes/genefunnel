<!-- Assisted-by: OpenAI Codex. -->

# Kang/Reactome aggregation result B-1.0.3

**Decision:** PASS for frozen technical stability; FAIL for the combined
held-out replication gate and both biological-effect endpoints. This completes
the donor-level perturbation characterization without supporting perturbation
sensitivity or a biological effect claim. Every endpoint is retained in the
[`machine result`](aggregation-kang-result.tsv).

## Provenance and execution

- Protocol SHA-256:
  `c95a2eb91c9e3f8027b461e85ee532a621550038b27f20167b0838dd95c2f7ad`.
- Data-manifest SHA-256:
  `c8743b696f3e05fa623d7f3adc19e27b379303e7f0544a88e054aedf528d5a28`.
- Runner/package source: clean commit
  `62d42dc3652d810b3c7ea58c70d29986b42504ce`; generated 2026-07-18
  08:48:00 UTC with one worker and C collation.
- Inputs: four exact pinned files, 77,527,715 bytes. The two sparse matrices
  had dimensions 35,635 x 14,619 and 35,635 x 14,446 with 8,732,747 and
  8,838,098 coordinate entries; neither contained a duplicate/explicit-zero
  entry after canonicalization.
- Environment: Linux x86_64, R 4.6.1, genefunnel 0.99.0, Matrix 1.7.5,
  BiocParallel 1.46.0, Rcpp 1.1.2, and RcppArmadillo 15.4.0.1.
- The artifact manifest covers 19 files and 14,292,063 bytes;
  `artifacts.tsv` SHA-256 is
  `5247a8c71abfca58d929f0ce17e56c25fd8dbbcdc893fcec788529c033b21e2f`.
  All listed hashes match. Independent barcode/unit and Reactome membership
  reconstructions were exact; endpoint signs and tests reconstructed from the
  retained tables. Maximum relative score-identity residual was `9.30e-15`.
  Rounded audit-table rank ties changed reconstructed correlations by at most
  `4.34e-5`, far from either cutoff. A second isolated execution reproduced all
  12 scientific tables byte-for-byte; SHA-256 of their ordered hash stream is
  `040d2beca3c0bf41314c6c152fdffd666601bd037e023cab9dfd09aaa837168b`.

## Preprocessing and support

- The 29,065 raw barcodes included 313 cross-file duplicates. Frozen
  `make.unique(..., sep = "")` output exactly matched metadata order. Nine
  records lacked author cell-type labels and were excluded by the registered
  cell-type filter, leaving 23,981 registered singlet candidates.
- Of the fixed 96 donor/condition/cell-type units, 84 met the 40-total/20-half
  rule. Their 23,649 cells entered all views. Twelve units containing 332
  candidate cells failed eligibility and entered none; eligible units ranged
  from 40 to 1,484 cells.
- All 35,635 gene rows had symbols; first-occurrence collapse produced 32,938
  symbols and summed 2,697 duplicate-symbol rows. Of 2,868 human Reactome
  pathways, 1,728 retained 8-128 measured symbols (59,969 memberships). The
  two primary pathways retained 76 and 97 symbols.
- The fixed design produced 82,944 pathway/view rows: 1,728 pathways x 16
  donor-conditions x three views. Every row was audit-eligible. Normalized gap
  was defined for 80,633 rows; 2,311 exact-zero aggregates correctly retained
  `NA`. Across defined rows, median/type-8 90th-percentile normalized gaps were
  0.08096/0.28788.

## Frozen decisions

All 16 donor-condition split correlations were defined across 1,621-1,690
common pathways. Correlations ranged 0.68118-0.93642. Their median was 0.91661
and type-8 10th percentile 0.83992, passing thresholds 0.70 and 0.50.

- `R-HSA-909733` (interferon alpha/beta signaling): training median contrast
  was +0.02340. The positive direction repeated in 3/4 held-out donors, so the
  descriptive held-out endpoint passed. Six of eight donor contrasts were
  positive; exact two-sided sign-test `p = 0.2890625`, Holm-adjusted
  `p = 0.578125`, so the biological endpoint failed.
- `R-HSA-877300` (interferon gamma signaling): training median contrast was
  -0.001880. The negative direction repeated in only 2/4 held-out donors, so
  replication failed. Four of eight contrasts were positive; exact and
  Holm-adjusted `p = 1`, so the biological endpoint also failed.

Technical split repeatability therefore does not establish donor-replicated
perturbation sensitivity. Together with the prior CellBench failure, this
rejects public audit promotion under B-1.0.3 and closes Workstream B with the
theorem-plus-warning fallback. Cells and pathways were never substituted for
the eight donor units.
