<!-- Assisted-by: OpenAI Codex. -->

# Controlled sensitivity result E-1.0.0

**Decision:** FAIL both frozen co-primary prediction gates. Retain the exact
diagnostic as an internal adversarial-test oracle and omit a public sensitivity
API. This controlled result cannot support an external technical-repeat,
biological-replicate, causal, or assay-independent reliability claim.

## Provenance and completeness

- Parent `E-1.0.0` registry MD5:
  `09d9429ae18f64ba00b3bd84955ad71d`; execution supplement `E-C-1.0.0`
  MD5: `1c119004ee749b4242f1b73ca6fd8c4c`.
- Clean runner/package candidate:
  `5920ea920e0e0e85659485c26f8096a61f08ea40`; source archive SHA-256:
  `ba49307a94b8cb1fd347db910be71975a5e39ed92b7d7a1ff65557f509663013`;
  installed manifest MD5: `fe6476f760cf6aedc9343b8299ad5088`.
- Generated 2026-07-18 12:25:03 UTC with eight workers and 128-scenario
  chunks: 5,760 latent scenarios, 345,600 feature-loss rows, 5,760 repeat rows,
  ten folds of 576 scenarios, and 2,000 bootstraps per target. No fixed row or
  status was excluded; 66 duplicated encoding rows had valid `zero_total`
  partial inputs.
- The ignored full bundle has 20 core artifacts totaling 251,274,443 bytes
  plus 45 atomic checkpoints. `artifacts.tsv` SHA-256 is
  `73872ce7576643fef4db48f9994e1b17180647117887090191f6c19423ec2038`;
  the sorted checkpoint-hash stream SHA-256 is
  `0a6ebe60a6d5c22c263eecd17f53de3be2c5b36804533d533de386af7e87296e`.
  Every listed artifact hash re-verifies. Independent TSV endpoint
  recomputation differs by at most `6.6e-17`, below one binary64 epsilon.

## Frozen endpoints

| Target | Median fold RMSE reduction | Bootstrap 95% lower | Required | Result |
|---|---:|---:|---:|---|
| Feature loss | 0.00110097 | 0.000407703 | >=0.10; >=0.05 | fail |
| Controlled repeat | 0.0193130 | 0.0112902 | >=0.10; >=0.05 | fail |

Feature-loss fold reductions ranged from -0.000674 to 0.001463; controlled-
repeat reductions ranged from -0.002184 to 0.029976. Every feature stratum lay
between -0.003506 and 0.003538; every repeat stratum lay between -0.002411 and
0.029976. No factor, fraction, mechanism, encoding, or fold approaches the
frozen primary thresholds, and descriptive strata cannot rescue them. Exact
machine endpoints are in
[`sensitivity-controlled-result.tsv`](sensitivity-controlled-result.tsv).

## Artificial thinning curves

Feature-loss median/90th-percentile error across all mechanisms was
0.134718/0.338044. It rose monotonically with removed fraction. At 50% removal,
low-abundance masks had median/90th-percentile loss 0.165583/0.344220, while
high-abundance masks had 0.361659/0.442211; low- and high-detection masks lay
between them. Paired global-absence and sample-missing rows have identical
observed values, scores, targets, and sensitivities but retain their distinct
coverage facts. The 30 fixed curve rows are retained in
[`sensitivity-controlled-curves.tsv`](sensitivity-controlled-curves.tsv).

These masks depend on the simulated profile and model-implied detection
probability. They are study-composition-dependent artificial deletions, not
evidence that real unmeasured genes are zero, missing at random, or governed by
the same detection model. The curves describe the failure envelope only; they
are not rescue endpoints.

The repeat target itself had median/90th percentile 0.052547/0.207585, but the
compact sensitivities added only the failed incremental prediction above the
frozen score/coverage/sum/support baseline. Algebraic correctness therefore
survives while empirical reliability and public promotion fail.
