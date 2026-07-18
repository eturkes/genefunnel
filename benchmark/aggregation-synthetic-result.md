<!-- Assisted-by: OpenAI Codex. -->

# Aggregation synthetic result B-1.0.2

**Decision:** PASS for all frozen controlled co-primary gates. This is
synthetic evidence for a narrow recovery task, not biological validation or
evidence that the normalized aggregation gap is robust to severe dropout.
Every exact endpoint is retained in the
[`machine result`](aggregation-synthetic-result.tsv).

## Provenance and execution

- Frozen protocol SHA-256:
  `aec8fd4e3e49b953e5ca75e0c2059c1d68436409def0dfd3d351b4ad8a49356f`.
- Frozen data-manifest SHA-256:
  `ca6dedb4957fa6ef3c21649d342996d1d482cf6f3203e4f2f9bab0bbf1a823e8`.
- Runner/package source: clean commit
  `d1cbf1510af3bc27d6baf14e3433866449645a77`; generated 2026-07-18
  07:12:03 UTC with four workers and chunk size 128.
- Environment: Linux x86_64, R 4.6.1, genefunnel 0.99.0, Matrix 1.7.5,
  BiocParallel 1.46.0, Rcpp 1.1.2, and RcppArmadillo 15.4.0.1.
- Design: 62,208 shared latent scenarios, two independently measured
  replicates, 124,416 rows, and 2,000 paired prediction-error bootstraps.
  Protocol-defined count sampling, post-count dropout, and per-unit rescaling
  to 1,000,000 were applied; no external biological data were used.
- Eligibility: 124,331 measurements eligible. Eighty-five zero-observed-unit
  measurements affected 84 pairs. One further fully eligible pair had exact
  zero aggregate score and therefore undefined normalized gap. The model
  excluded those 85 pairs as frozen.
- Ignored full bundle: 501 files and 68,034,582 bytes. `artifacts.tsv` SHA-256
  is `d413d40f5f1b2160ccadc4fa2a552d3344df5bb39d071587c69485ba463f8c4f`;
  SHA-256 of the sorted `sha256sum checkpoints/*.rds` stream is
  `de9c1075214145dc560d7883f7734a3c2e365dc31912b9acbc48c91f1b1f2cf2`.
  A binary replay verified all 486 checkpoint fingerprints and recomputed
  endpoints, folds, predictions, coefficients, and bootstrap values with
  maximum text-transport delta below `5e-16`; all 14 listed artifact hashes
  matched.

## Frozen endpoints

| Gate | Estimates | Frozen requirement | Result |
|---|---|---|---|
| Curve recovery | median error 0.0102821; 90th percentile 0.182020 | <=0.10; <=0.25 | pass |
| A/B stability | overall 0.929844; archetypes 0.919627, 0.930966, 0.936522 | >=0.80; each >=0.65 | pass |
| Null leakage | median 0.000149946; 95th percentile 0.0182094 | <=0.05; <=0.20 | pass |
| Incremental prediction | median fold RMSE reduction 0.218438; bootstrap lower bound 0.210229 | >=0.10; >=0.05 | pass |

## Failure envelope

The co-primary curve gate pools its protocol-defined non-complex domain. Its
factor-stratified diagnostics reveal a material severe-dropout limitation:

- At dropout 0, median/90th-percentile error was 0.00147/0.01228 and median
  observed/latent gap was 0.05584/0.05286. At dropout 0.5, these were
  0.04918/0.30464 and 0.15346/0.05333. The stratified 90th percentile therefore
  exceeds the overall gate's 0.25 cutoff, although it was not a separate frozen
  endpoint.
- The marginal median observed-gap rise from dropout 0 to 0.5 was 0.09762. The
  corresponding rise from planted complementarity 0 to 1 was 0.10728, so the
  qualitative PLAN fallback concern - technical dropout comparable with the
  planted effect - is present.
- One pooled subunit produced median/90th-percentile error 0.04945/0.27079;
  16 subunits reduced it to 0.00521/0.02659.

Consequently, the result supports the exact theorem and a warning about
aggregate-then-score behavior, but not a general robustness or biological
claim. The audit remains internal. Pinned cross-platform RNA mixtures and
donor-replicated perturbation data are the next independent gates.
