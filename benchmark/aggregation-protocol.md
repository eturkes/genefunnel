<!-- Assisted-by: OpenAI Codex. -->

# Aggregation-audit protocol B-1.0.3

**Frozen:** 2026-07-18. `B-1.0.0` preceded implementation; `B-1.0.2` preceded
every synthetic result; `B-1.0.3` precedes every external-data endpoint. The
[`machine registry`](aggregation-protocol.tsv) and
[`data manifest`](aggregation-data.tsv) are authoritative. An endpoint,
threshold, split, factor, source, preprocessing rule, or schema change requires
a new version committed before replacement results.

This protocol separates four questions:

1. Does an additive prototype implement the theorem and eligibility contract?
2. Does the observed gap recover latent cancellation under controlled noise?
3. Does it track known RNA mixtures across two laboratory protocols?
4. Is it technically stable and donor-replicated in an IFN-beta experiment?

The theorem can remain true when every empirical gate fails. Synthetic and
observational results support only the claims assigned below.

`B-1.0.1` is a schema-only implementation-review correction. `B-1.0.0`
promised the affected units required by `AGGREGATION_SPEC.md`, but its
`removed_members` table had no unit column. The first local red/green prototype
fixtures exposed that guarantee gap before any synthetic, downloaded-data, or
threshold result. The table now records `unit`; unmatched members use `NA`, and
missing members get one row per affected unit. Data bytes, designs, splits,
endpoints, thresholds, decision rules, and every other field remain unchanged.

`B-1.0.2` closes execution degrees of freedom found while designing the first
synthetic runner: it separates one shared latent seed from two measurement
seeds; defines the overlap, subunit/dropout, zero-total, and covariate
operators; and fixes fold assignment, factor encoding, QR alias handling,
quantiles, and a paired prediction-error bootstrap. It changes no API field,
factor level, external byte/URL, endpoint, threshold, or decision rule. This
amendment was committed before generating a synthetic measurement, fitting a
validation model, or reading downloaded data.

`B-1.0.3` closes external execution degrees of freedom before computing a pure
profile, selecting a data-derived set, running an external-data audit, or
calculating any CellBench/Kang endpoint. It fixes CellBench zero-abundance
ranking, equal-mixture normalization, gate units, condition grids, and split
ordering; and Kang's duplicated cross-file barcode join, sparse-count
validation, cell-type eligibility, pathway parsing, direction, and exact
sign-test rules. Only checksums, CSV/Matrix Market/GMT schemas, dimensions,
identifier alignment, and condition counts were inspected first to establish
feasibility. All source URLs/payload bytes, synthetic fields/results, API
fields, thresholds, donor splits, pathways, and decision rules are unchanged.

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
3. `removed_members`: zero or more rows in group/set/declared-member/unit
   order; reason is `"unmatched"` or `"missing"`. Unmatched rows have
   `unit = NA_character_`; missing rows identify each affected active unit.
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
and `B` give exactly 124,416 observations. For latent integer ID `s`, their
measurement seeds are `18072027 + s` and `19072027 + s`; both share latent
seed `20072027 + s`. Generation resets base R's
`Mersenne-Twister`/`Inversion`/`Rejection` RNG for each seed, so worker count
cannot change a row.

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

The latent RNG draws all `u` values, then `e` in unit-major order. Complex-like
rows reuse the first `u` and first `e` row across units, preserving an exact
latent no-cancellation control before an outlier. When overlap is 0.5, let
$a_i=(-1)^{i-1}$ and $b_k=(-1)^{k-1}$; the first member half is multiplied by
$1+0.25a_ib_k$. Complex-like rows set every $b_k=1$, reusing the programme
across units. The outlier, when active, multiplies `x[1, 1]` by 20. Finally each
unit is rescaled to total 1,000,000. The latent target $R^*$ is the exact
independent reference audit of these compatible profiles.

For each measured unit, multinomial draws allocate the fixed library depth.
Subunit sizes are the integer quotient plus one for the first remainder
subunits. In unit then subunit order, each multinomial draw is followed by
independent member-level Bernoulli dropout that zeroes the whole post-count
member value; subunits are then summed. Any zero-total unit makes the
measurement ineligible. Otherwise each unit is rescaled to 1,000,000.
`raw_total` is the post-dropout count sum over units/members;
`detected_fraction` is the fraction of post-dropout unit-member cells above
zero.

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

All descriptive quantiles use R type 8. For one-indexed 32-run factorial and
1,944-row core indices `f` and `c`, the fold is
`1 + (((f - 1) + 3 * (c - 1)) %% 10)`; paired A/B measurements share it. Each
declared factor level differs by at most seven scenarios across folds. Modeling
retains only pairs with defined finite targets and predictors. The four
continuous baseline columns are aggregate score, raw total, detected fraction,
and maximum weight.
Every declared design field is a factor in registry level order:
`archetype`, member/unit count, weight profile, depth, dropout,
complementarity, baseline, dynamic range, log SD, correlation, overlap,
outlier, and subunit count. A global treatment matrix fixes columns across
folds. Base R `lm.fit` supplies QR least squares; aliased coefficients are zero
for prediction. Fold reduction is `(baseline RMSE - augmented RMSE) / baseline
RMSE`.

The 2,000 bootstrap replicates use seed `21072027`. Each replicate resamples
held-out B scenarios with replacement inside each fold, reuses the frozen
cross-validation predictions, recomputes ten reductions, and takes their
median. The one-sided 95% lower bound is the type-8 0.025 quantile. Models are
not refit inside this prediction-error bootstrap.

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

Each count table must have unique non-empty Ensembl row identifiers, and count
columns must be identical in order to metadata rows. Logical `outliers` rows
are removed. Within each platform, every retained library becomes counts per
million using its total over all count rows; cross-platform gene filtering
occurs afterward. The three pure unit profiles are arithmetic mean CPM across
all retained pure libraries for that platform and cell line, pooling mRNA
amounts. Common genes follow CEL-seq2 row order and must also occur in
SORT-seq. Gene selection uses only CEL-seq2 pure profiles in that universe:

- for each cell-line pair and size 8/32/128, take half the largest log2 ratio
  in each direction, advancing through the stable rank to avoid duplicates;
- for each size, a complex-like control takes genes with the largest minimum
  pure-profile abundance across all three lines.

For a pair `(first, second)`, positive/zero has ratio `Inf`, zero/positive has
`-Inf`, and zero/zero is the neutral ratio zero. Rank by decreasing signed
ratio then Ensembl identifier; choose the first-direction half, then advance
through the reverse-direction rank while skipping selected identifiers. Pair
order is H2228/H1975, H2228/HCC827, H1975/HCC827; sizes follow registry order.
Controls rank decreasing minimum pure CPM then identifier and follow the nine
pair sets. This fixes 12 sets without using mixed-sample outcomes.

For each mixed library $y$, the raw metadata triple must match a registered
composition exactly. It is normalized to sum one - in particular, the three
reported `0.33` values become equal audit weights. Those weights and platform
pure profiles $q_k$ give reference $R_{ref}$ from the audit. Pure libraries
train profiles only and never enter an endpoint. The observed process-control
quantity

$$
R_{obs}=\{F(y)-\sum_kw_kF(q_k)\}/F(y)
$$

is undefined when `F(y) <= 0` and otherwise may leave `[0, 1]`: a measured
library need not equal the reference weighted mean. It is never an API result.

The error unit is one mixed library/set absolute `R_obs - R_ref`. For every
fixed platform/composition/mRNA-amount/set group, its median must be `<= 0.15`
and type-8 90th percentile `<= 0.30`; the error gate requires every group.
There are four non-pure compositions x four amounts x 12 sets = 192 conditions
per platform. Within platform/composition/amount, C-locale sorted sample
identities split odd/even. The 192 matched odd/even condition/set medians must
have Spearman correlation `>= 0.75` separately in CEL-seq2 and SORT-seq. The
192 matched complete-condition medians must correlate `>= 0.60` between
platforms. A missing or undefined fixed-grid metric records scientific `FAIL`;
malformed input aborts. Error and stability gates are co-primary.

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
gene table and columns to their respective barcode files. The two raw barcode
files reuse 313 identifiers. Concatenate them in matrix order and apply base R
`make.unique(..., sep = "")`; that sequence must identically equal the author
metadata identifiers. Coordinate entries must be finite, non-negative integer
counts; duplicate coordinates are summed. Any mismatch rejects the run.

Empty gene symbols are discarded; duplicate symbols are summed in first-symbol
occurrence order. Registered donors, `ctrl`/`stim`, six registered cell types,
and author-labelled singlets are retained. Within each donor/condition/cell
type, barcodes sort in C locale and split by odd/even position. A cell type
needs at least 40 cells and 20 in each half; failure removes that unit from the
full and both half views. Within each view, raw UMI counts are summed and
divided by its retained cell count. Cell types are audit units and their cell
count proportions among retained group cells are weights. This reconstructs
the donor-condition/view mean UMI per cell without log transformation or
unit-specific library normalization.

[Reactome v97](https://reactome.org/download/97/ReactomePathways.gmt.zip) is
pinned rather than `current`. The archive must contain exact member
`ReactomePathways.gmt`. `R-HSA-` pathways retain unique first-occurrence member
symbols intersected with collapsed measured symbols; sets of size 8-128 remain
in GMT order. Full, odd, and even audits use `missing = "reject"`. The two
pre-specified perturbation pathways are
`R-HSA-909733` (interferon alpha/beta signaling) and `R-HSA-877300`
(interferon gamma signaling).

Technical stability is Spearman correlation of odd/even normalized gaps across
common defined pathways within each of the fixed 16 donor-conditions. All 16
must be defined. Their median must be `>= 0.70` and type-8 10th percentile
`>= 0.50`.

Training donors are 101, 1015, 1039, and 1256; held-out donors are 107, 1016,
1244, and 1488. Full-view donor contrasts are `stim - ctrl` normalized gap.
For each primary pathway, training fixes the sign of the four-donor median;
undefined or exact-zero direction fails. Replication requires the same nonzero
sign in at least three of four held-out donors. This is a held-out descriptive
gate, not inference. The separate two-sided exact sign test discards exact-zero
contrasts, uses a binomial null probability 0.5 (`p = 1` if none remain), and
Holm-adjusts across the two pathways. A biological effect claim requires the
held-out rule plus adjusted `p <= 0.05` for both pathways. Cells and pathways
are never substituted for the eight donor units.

## Decision and evidence

Prototype correctness, all synthetic primary gates, both CellBench gates, and
Kang technical stability must pass before an exported audit is proposed. The
held-out donor rule must also pass before documentation suggests perturbation
sensitivity. The adjusted sign-test gate is mandatory for any biological
effect claim.

The completed B-1.0.2 synthetic result is inherited unchanged: B-1.0.3 alters
no synthetic field, endpoint, threshold, byte, or implementation used by that
run. It does not reinterpret the severe-dropout failure envelope.

`run-aggregation-synthetic.R` requires a clean committed tree, installs that
tree into an isolated temporary library, checks deterministic generator/audit
and fabricated model-orchestration controls, and retains all observations,
folds, coefficients, bootstrap estimates, endpoints, failures, environment
facts, and artifact hashes. Checkpoints bind to protocol hash, Git commit, and
ordered scenario IDs. Scientific gate failure completes successfully and is
reported as `FAIL`; provenance/schema/arithmetic failure aborts.

The CellBench and Kang runners likewise require clean commits, verify only
manifest-pinned payloads, install the exact source snapshot, and retain fixed
grid failures. The Kang evidence includes exact matrix/barcode preprocessing,
all cell-type eligibility and split assignments, retained Reactome membership,
every pathway/view audit, 16 split correlations, donor contrasts, Holm-adjusted
decisions, environment facts, and artifact hashes. Missing or undefined
scientific endpoints record `FAIL`; malformed bytes, schema, alignment, or
arithmetic aborts.

The complete pinned execution passed Kang technical stability (median/type-8
10th-percentile split Spearman 0.91661/0.83992) but failed the combined held-out
gate: interferon gamma matched its training direction in only two of four
held-out donors. Holm-adjusted exact sign-test p-values were 0.578125 and 1, so
neither biological-effect endpoint passed. The tracked
[`result`](aggregation-kang-result.md) retains exact provenance and all six
endpoints. Together with CellBench failure, this completes negative external
validation and selects PLAN's fallback without exporting the prototype.

Failure invokes PLAN's fallback: retain the theorem and aggregate-then-score
warning; keep or remove the internal prototype according to maintenance value.
No threshold is relaxed after seeing results. Generated downloads, expanded
manifests, raw observations, reports, and decisions stay under ignored
`benchmark/results/` or `.agent/tmp/`; tracked evidence is compact and names
exact protocol/data hashes, Git state, R/session state, preprocessing, excluded
facts, and every failed as well as passed endpoint.
