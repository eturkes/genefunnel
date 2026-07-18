<!-- Assisted-by: OpenAI Codex. -->

# CellBench aggregation result B-1.0.3

**Decision:** FAIL for both frozen CellBench co-primary gates. This is a
known-mixture process-control result, not a biological effect test. Every exact
endpoint and its fixed-grid completeness count is retained in the
[`machine result`](aggregation-cellbench-result.tsv).

## Provenance and execution

- Protocol SHA-256:
  `c95a2eb91c9e3f8027b461e85ee532a621550038b27f20167b0838dd95c2f7ad`.
- Data-manifest SHA-256:
  `c8743b696f3e05fa623d7f3adc19e27b379303e7f0544a88e054aedf528d5a28`.
- Runner/package source: clean commit
  `9b3cab65f626200d28a9ce494b79700247eaf0be`; generated 2026-07-18
  08:02:47 UTC with one worker and C collation.
- Inputs: four exact pinned files, 4,072,960 bytes; 13,906 common Ensembl
  identifiers; 12 CEL-seq2-pure-profile-derived sets; 4,848 mixed
  library/set rows. No mixed outcome or SORT-seq profile selected a member.
- Environment: Linux x86_64, R 4.6.1, genefunnel 0.99.0, Matrix 1.7.5,
  BiocParallel 1.46.0, Rcpp 1.1.2, and RcppArmadillo 15.4.0.1.
- Ignored full bundle: 17 files and 5,314,921 bytes. `artifacts.tsv` SHA-256
  is `beff1ccd5c86cbd52608c4fcaf4aaecbe4d9d7b6964fe1c91b15956ac71ed11b`.
  Independent table reconstruction matched all 16 listed artifact hashes;
  selection metrics were exact and maximum decision-table transport delta was
  `5.41e-13`. A second isolated execution reproduced the 12 scientific tables
  byte-for-byte; SHA-256 of their ordered hash stream is
  `bcc21fc6c2d8f1dcaa82105fa623bec2e126e6aba5f970be23d8127212346015`.

## Failure anatomy

All 96 platform/composition/set reference audits were eligible with defined
normalized gaps. Failure arose after measuring individual mixtures:

- 1,434 of 4,848 library/set scores were exactly zero, making the observed
  ratio undefined. Only 191 of 384 fixed groups were complete. The frozen
  maximum-group curve endpoints are therefore `NA` and fail.
- Of 384 groups, 299 failed: every one of the 288 pair-set groups plus 11 of
  96 complex-control groups. Only 85 complex groups passed.
- Failure was not rescued by defined rows. Across pair sets of size 8/32/128,
  defined observations numbered 182/827/1,193 of 1,212, while their pooled
  median absolute errors were 0.361/3.700/5.224 and type-8 90th percentiles
  were 1.565/16.007/27.714. The largest finite group median/q90 was
  75.910/168.355.
- Complex controls were materially better: pooled median/q90 errors for sizes
  8/32/128 were 0.0809/0.2218, 0.0647/0.1683, and 0.0254/0.0654.
- Complete-case stability estimates exceeded their numerical cutoffs, but the
  fixed grid did not: CEL-seq2 had 98/192 complete conditions, SORT-seq 93/192,
  and the cross-platform comparison 88/192. Their correlations were
  0.97798, 0.94450, and 0.90995, respectively; all three endpoints fail the
  pre-specified completeness rule.
- Finite observed ratios ranged from -248.794 to 0.834, versus reference ratios
  from 0.00120 to 0.590. This is allowed by the protocol and exposes the
  instability of dividing by a small measured library score; it is not clamped.

The result rejects public promotion under B-1.0.3 and reinforces PLAN's
theorem-plus-warning fallback. Kang/Reactome remains useful for completing the
pre-specified technical and donor characterization, but cannot rescue the
failed CellBench export gate.
