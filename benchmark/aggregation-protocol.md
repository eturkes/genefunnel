<!-- Assisted-by: OpenAI Codex. -->

# Aggregation-audit protocol B-1.0.0

**Frozen:** 2026-07-18, before a group-audit implementation or validation
result. [`aggregation-protocol.tsv`](aggregation-protocol.tsv) is the machine
registry and [`aggregation-data.tsv`](aggregation-data.tsv) pins every external
byte. An endpoint, threshold, split, factor, source, preprocessing rule, or
schema change requires a new version committed before replacement results.

This protocol separates four questions:

1. Does an additive prototype implement the theorem and eligibility contract?
2. Does the observed gap recover latent cancellation under controlled noise?
3. Does it track known RNA mixtures across two laboratory protocols?
4. Is it technically stable and donor-replicated in an IFN-beta experiment?

The theorem can remain true when every empirical gate fails. Synthetic and
observational results support only the claims assigned below.

## Prototype schema

The first implementation is internal `.aggregation_audit()` with fixed
signature

```r
.aggregation_audit(mat, gene_sets, groups, weights, missing = "reject")
```

It remains unexported until all public-promotion gates pass. `mat` uses the
existing feature-by-unit base/`Matrix` storage contract. It additionally needs
unique non-empty column names. `groups` is an unclassed character vector and
`weights` an unclassed numeric/integer vector; both are named exactly and in
the same order as `colnames(mat)`. Group order is first occurrence. Input
weights are finite non-negative masses, normalized separately within each
group and always reported. Equal weights therefore require explicit equal
masses.

Matrix structure, identifiers, group labels, and all weights are validated
first. A zero-weight unit is then excluded before its matrix values are
validated. Every positive-weight column must satisfy the scorer's finite,
non-negative-or-missing value contract across all features. Every declared set
and group receives one summary row; sets below two matched members are
reported rather than silently omitted.

`missing = "reject"` makes a group/set ineligible when any matched member is
`NA`/`NaN` in an active unit. `missing = "intersection"` removes a member when
any active unit is missing it, preserves declared first-match order, recomputes
all terms on that support, reports the member, and rejects fewer than two.
Unmatched declared members are also reported, separately from missing ones.

The result is a fixed-order unclassed list of four unclassed data frames:

1. `summary`: columns fixed in the machine registry. One row per group then
   declared set.
2. `weights`: one row per input unit, preserving matrix order. All-zero groups
   have `effective_weight = NA`; otherwise inactive units have zero.
3. `removed_members`: zero or more rows in group/set/declared-member order;
   reason is `"unmatched"` or `"missing"`.
4. `unit_scores`: one row per eligible group/set/active unit, preserving unit
   order.

An eligible row has reason `"eligible"`. Fixed ineligible reasons are
`"no_positive_weight"`, `"too_few_matched_members"`,
`"active_missing_values"`, `"too_few_common_members"`, and
`"numerically_unavailable"`. Ineligible scalar metrics and
`identity_residual` are `NA_real_`; factual counts remain populated.
`normalized_status` is `"defined"`, `"zero_aggregate"`, or `"ineligible"`.

### Numerical contract

Unit and aggregate scores use the authoritative native score on the exact
retained block. The gap is independently evaluated from the non-negative
opposing-centered-mass formula in `inst/AGGREGATION_SPEC.md`, after division by
the block's largest value and final homogeneous rescaling. The returned
identity residual is

$$
F(\bar{x})-\sum_kw_kF(x_k)-J.
$$

Let `s` be the unscaled maximum value, `m` active units, and `n` retained
members. The fixed binary64 allowance is

```text
tau = 256 * eps * max(1, m, n) *
      max(F(mean), weighted_unit_score, J, s)
```

where `eps = .Machine$double.eps`. Defined finite metrics require
`abs(identity_residual) <= tau`, `J >= 0`, and `J <= F(mean) + tau`.
Unrepresentable rescaling or a failed identity produces
`"numerically_unavailable"`; the implementation neither clamps a negative
direct subtraction nor substitutes a rounded identity. Exact zero aggregate
gives scores/gap zero, normalized gap `NA`, and status `"zero_aggregate"`.

API correctness is an intersection gate: schema/types/order, authoritative
scores, independent oracle agreement, missingness modes, every reason,
zero-weight exclusion, dense/sparse identity, scaling, physical sums, and
unit/member/group/set permutations all pass. Performance is descriptive for
this internal prototype and cannot change the frozen scorer.

## Controlled synthetic experiment

### Design

The machine generator crosses seven core factors:

- archetype: `complex_like`, `cascade_like`, `regulatory_like`, or
  `mixed_direction`;
- retained members: 8, 32, or 128;
- active units: 2 or 4;
- weight profile: equal, moderate, or dominant;
- per-unit library depth: 1,000, 10,000, or 100,000;
- independent post-count dropout: 0, 0.2, or 0.5; and
- complementarity fraction $\lambda$: 0, 0.5, or 1.

For two units, weight profiles are `(0.5, 0.5)`, `(0.25, 0.75)`, and
`(0.1, 0.9)`. For four they are `(0.25, 0.25, 0.25, 0.25)`,
`(0.1, 0.2, 0.3, 0.4)`, and `(0.05, 0.05, 0.1, 0.8)`.

Each core row crosses a 32-run resolution-IV fractional factorial. Base signs
`A:E` enumerate $2^5$ in lexicographic order; `F = A*B*C` and `G = B*C*D`.
Signs select, respectively:

| Sign | low | high |
|---|---|---|
| A - baseline | monotone | U-shaped |
| B - member dynamic range | 4-fold | 64-fold |
| C - log-scale SD | 0.1 | 0.5 |
| D - member common-factor fraction | 0 | 0.6 |
| E - nuisance-program overlap | 0 | 0.5 |
| F - one-unit/member outlier multiplier | 1 | 20 |
| G - measured subunits pooled per unit | 1 | 16 |

This gives 62,208 latent scenarios. Independent measurement replicates `A`
and `B` give exactly 124,416 observations. Seeds are `18072027 + scenario row`
for replicate A and plus 1,000,000 for B. Scenario generation uses base R only.

### Latent profiles

For $t_i=(i-1)/(n-1)$ and dynamic range $r$, the baseline is either
$q_i=r^{t_i-1/2}$ or $q_i=r^{2|t_i-1/2|-1/2}$, then rescaled to mean one.
Let $z_{0i}=\sin(2\pi t_i)$. Unit patterns are:

- complex-like: $z_{ki}^{A}=z_{0i}$;
- cascade-like: a mean-centered Gaussian peak
  $\exp[-16(t_i-(k-1/2)/m)^2]$;
- regulatory-like: $\cos[2\pi(t_i-(k-1)/m)]$; and
- mixed-direction: $(-1)^{k-1}$ times `+1` for the first member half and `-1`
  for the second.

Each archetype row is divided by its maximum absolute entry. Define
$z_{ki}=(1-\lambda)z_{0i}+\lambda z_{ki}^{A}$ and initial
$x_{ki}=q_i(1+0.8z_{ki})$, which remains positive.

Seeded log-normal perturbation uses SD $\sigma$, correlation $\rho$, member
loading $l_i=0.5+t_i$, unit normal $u_k$, and independent $e_{ki}$:

$$
x_{ki}\leftarrow x_{ki}
\exp\{\sigma[\sqrt{\rho}\,u_kl_i+\sqrt{1-\rho}\,e_{ki}]
-\sigma^2[\rho l_i^2+(1-\rho)]/2\}.
$$

Complex-like rows reuse one perturbation vector across units, preserving an
exact latent no-cancellation control before an outlier. When overlap is 0.5,
the first half of members is multiplied by a bounded second alternating
unit/member programme in `[0.75, 1.25]`; complex-like rows reuse that programme
across units. The outlier, when active, multiplies `x[1, 1]` by 20. Finally each
unit is rescaled to total 1,000,000. The latent target $R^*$ is the exact
reference audit of these compatible profiles.

For each measured unit, multinomial draws allocate the fixed library depth.
The 16-subunit condition divides the integer depth as evenly as possible,
applies dropout independently within each subunit, and sums subunits. The
one-subunit condition applies dropout once. Nonzero observed totals are
rescaled to 1,000,000; total-zero rows are recorded ineligible. Raw total and
detected-member fraction remain covariates.

### Primary endpoints

- Curve error: `abs(observed R - latent R*)` over defined non-complex rows.
  Median must be `<= 0.10` and the 90th percentile `<= 0.25`.
- Replicate stability: Spearman correlation of A/B observed R over matched
  rows must be `>= 0.80` overall and `>= 0.65` separately for each non-null
  archetype.
- Null leakage: among complex-like, no-outlier, zero-dropout rows, median R
  must be `<= 0.05` and the 95th percentile `<= 0.20`.
- Incremental effect: predict latent $R^*$ with fixed 10-fold folds. Baseline
  predictors are aggregate score, raw total, detected fraction, declared
  depth/dropout, member/unit count, maximum weight, and all design-factor main
  effects. The augmented model adds absolute and normalized aggregation gap.
  Train on replicate A outside a fold and predict replicate B inside it. Median
  fold RMSE reduction must be `>= 10%`; a 2,000-replicate scenario-stratified
  bootstrap 95% lower bound must be `>= 5%`.

All four are co-primary intersection gates. Factor-stratified errors and
coefficients are reported without rescuing a failed gate. The 0.5-dropout
stratum is retained even if it rejects the method: PLAN's fallback explicitly
asks whether missingness can dominate planted cancellation.

## Known RNA-mixture validation

[Tian et al.](https://doi.org/10.1038/s41592-019-0425-8) created RNA mixtures
of H2228, H1975, and HCC827 at three pure, three 68/16/16, and one equal
composition, across 3.75, 7.5, 15, and 30 pg inputs. The original repository is
pinned at one commit in the data manifest. CEL-seq2 is training; SORT-seq is a
held-out laboratory protocol. Rows marked `outliers` are excluded by the
authors' supplied metadata. Counts and metadata must align exactly by sample.

Within each platform, every library becomes counts per million. Pure-sample
means define the three unit profiles. Gene selection uses only CEL-seq2 pure
profiles and genes present in both platforms:

- for each cell-line pair and size 8/32/128, take half the largest log2 ratio
  in each direction, advancing through the stable rank to avoid duplicates;
- for each size, a complex-like control takes genes with the largest minimum
  pure-profile abundance across all three lines.

Ties resolve by Ensembl identifier. This fixes 12 data-derived sets without
using mixed-sample outcomes.

For each mixed library $y$, known proportions $w$ and platform pure profiles
$q_k$ give reference $R_{ref}$ from the audit. The observed process-control
quantity

$$
R_{obs}=\{F(y)-\sum_kw_kF(q_k)\}/F(y)
$$

is allowed outside `[0, 1]`: a measured library need not equal the reference
weighted mean. It is never presented as an API result. Median absolute
`R_obs - R_ref` by platform/composition/mRNA amount/set must be `<= 0.15`; its
90th percentile must be `<= 0.30`. Sorted sample identities split technical
replicates odd/even. The split-half condition/set medians must correlate
`>= 0.75` within each platform, and complete-condition medians must correlate
`>= 0.60` between platforms. Both error and stability gates are co-primary.

## Donor-replicated perturbation validation

[Kang et al.](https://doi.org/10.1038/nbt.4042) profiled PBMCs from eight
donors with paired control/IFN-beta conditions; the exact GEO
[GSE96583](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE96583) files
are pinned. Only author-labelled singlets from B cells, CD14+ monocytes, CD4 T
cells, CD8 T cells, FCGR3A+ monocytes, and NK cells enter. No learned cell
annotation or cell-level inference is added.

From the raw tar, the runner extracts only `GSM2560248_2.1.mtx.gz`,
`GSM2560248_barcodes.tsv.gz`, `GSM2560249_2.2.mtx.gz`, and
`GSM2560249_barcodes.tsv.gz`. Matrix rows align exactly to the pinned batch-2
gene table; matrix columns align exactly to the respective barcodes and joined
author metadata. Any mismatch rejects the run.

Within donor, condition, and cell type, raw UMI counts are summed and divided
by retained cell count. Cell types are the audit units; cell-count proportions
are their weights. This exactly reconstructs the donor-condition mean UMI per
cell and performs no log transform or unit-specific library normalization.
Duplicate non-empty gene symbols are summed. A cell type needs at least 40
cells, with at least 20 in each deterministic sorted-barcode odd/even half.

[Reactome v97](https://reactome.org/download/97/ReactomePathways.gmt.zip) is
pinned rather than `current`. Unique human-symbol sets retaining 8-128 measured
members are audited. The two pre-specified perturbation pathways are
`R-HSA-909733` (interferon alpha/beta signaling) and `R-HSA-877300`
(interferon gamma signaling).

Technical stability is Spearman correlation of half-1/half-2 R across eligible
pathways within each donor-condition. Its median across 16 donor-conditions
must be `>= 0.70` and 10th percentile `>= 0.50`.

Training donors are 101, 1015, 1039, and 1256; held-out donors are 107, 1016,
1244, and 1488. For each primary pathway, training fixes the sign of median
paired `stim - ctrl` R. Replication requires the same nonzero sign in at least
three of four held-out donors. This is a held-out descriptive gate, not an
inferential claim. A biological perturbation claim additionally requires a
two-sided exact paired sign-test Holm-adjusted `p <= 0.05` across the two
pathways. Cells and pathways are never substituted for the eight donor units.

## Decision and evidence

Prototype correctness, all synthetic primary gates, both CellBench gates, and
Kang technical stability must pass before an exported audit is proposed. The
held-out donor rule must also pass before documentation suggests perturbation
sensitivity. The adjusted sign-test gate is mandatory for any biological
effect claim.

Failure invokes PLAN's fallback: retain the theorem and aggregate-then-score
warning; keep or remove the internal prototype according to maintenance value.
No threshold is relaxed after seeing results. Generated downloads, expanded
manifests, raw observations, reports, and decisions stay under ignored
`benchmark/results/` or `.agent/tmp/`; tracked evidence is compact and names
exact protocol/data hashes, Git state, R/session state, preprocessing, excluded
facts, and every failed as well as passed endpoint.
