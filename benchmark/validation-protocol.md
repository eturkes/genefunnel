<!-- Assisted-by: OpenAI Codex. -->

# Prospective scientific-validation protocol F-2.0.0

**Frozen:** 2026-07-18, before comparative implementation, dataset selection,
molecular-value access, or scoring. The
[`machine registry`](validation-protocol.tsv) is authoritative. This document
explains its estimands and decision rules. It is a design, not a result.

Protocol `1.0.0` remains the immutable analytic/resource control. `F-2.0.0`
introduces methods, assays, and inferential assertions, so it starts a new major
benchmark protocol. It is intentionally not executable or present in the
version-selecting index yet. Source, task, and role selection must be committed
before retrieval. `F-2.0.1` then pins exact studies, tasks, roles, licenses,
URLs, bytes, and checksums after opaque retrieval but before an archive or
molecular table is parsed or opened.
`F-2.0.2` must pin the complete generator, matching, schema, and runner closure
after development-data plumbing but before held-out access or the first
comparative score. Any scientifically material change gets a new version and a
prospective rationale.

## Domain and units

The primary domain is human bulk RNA-seq, donor-level/source-labelled
pseudobulk RNA-seq, and bulk proteomics. Accepted RNA inputs start as raw
gene-level counts. Pseudobulk counts are summed within biological unit,
condition, and source-supplied cell type, with at least 20 cells in every such
stratum. Accepted proteomics inputs start as linear non-negative protein
abundances or intensities whose zero and missing sentinels are documented.
Missing proteomics measurements become `NA`, never zero. A recorded zero is
accepted only when the source defines it as a quantified zero.

The inference unit is a donor, organism, independent culture, or independent
biospecimen - never a cell or technical repeat. Technical repeats are aggregated
before primary scoring and retained only for the secondary stability endpoint.
Each biological unit is scored before a condition contrast is formed; pooling
units and then scoring would change the GeneFunnel estimand.

Direct single-cell use is a different hypothesis. Technical zeros at cell level
can violate the observed-zero contract, and bulk methods have already shown
context-dependent cell-level behavior. It therefore requires a later protocol
rather than a favourable subgroup inside this one ([Noureen et al.,
2022](https://elifesciences.org/articles/71994)).

Each task is one study, assay, intervention contrast, and prespecified Reactome
target. Eligibility requires at least four independent units per condition or
four complete pairs, a concurrent control, an externally supported mechanistic
target and expected direction, resolvable metadata/covariates, and no perfect
condition-batch confounding. Source-declared QC rules may exclude samples only
when recorded before scores. Score- or endpoint-driven exclusion is invalid.

## Prospective split

At least one development study per assay stratum is used only to prove plumbing
and is absent from every claim endpoint. Each stratum then needs at least three
independent held-out studies, two tasks per study, six tasks total, and two
mechanistically distinct targets. Shared samples or alternate processing of the
same raw material count as one study. Source, task, and role selection is
committed before retrieval. Opaque archives may then be downloaded to establish
exact byte hashes; the data supplement is committed before parsing or
molecular-value access. Studies are preprocessed and scored independently; no
cross-study normalization or batch correction is allowed.
All candidate studies and their order are fixed up front. A source that becomes
ineligible after access is retained as an exclusion and is not replaced.

Development can expose implementation defects, but cannot move the scientific
thresholds below. After the implementation supplement is frozen, held-out data
are verified and the gate runs once. An unavailable eligible task, failed
coverage match, or undefined method output remains visible and counts against
the minimum replication requirements.

## Preprocessing and meaningful zero

RNA Ensembl version suffixes are stripped; only human `ENSG` identifiers are
accepted and duplicate rows are summed. Genes with zero total count across the
study are removed. [edgeR 4.10.1](https://bioconductor.org/packages/edgeR/)
performs TMM normalization through the exact calls in the registry. The primary
matrix is unlogged TMM CPM, which is non-negative and preserves a raw zero. It
does not turn a technical non-detection into evidence of biological absence.

Proteomics identifiers reduce to one canonical human UniProt accession.
Ambiguous protein groups are excluded; duplicate quantifications need a
source-defined, non-overlapping linear aggregation rule in the data manifest.
Each sample's linear values are scaled to a total of one million over all its
non-missing retained proteins before common-support filtering. Primary
comparisons use proteins observed in every unit and no imputation.

After assay preprocessing, rows constant across all primary units are removed.
Every primary method receives the identical matrix and identical ordered set
members. This same-compatible-input regime isolates the scoring rules. A
secondary method-native regime lets GSVA/ssGSEA use TMM log-CPM for RNA and a
source-provided normalized log2 proteomics assay where one exists; GeneFunnel,
sum, and mean remain on the compatible non-negative matrix. Native-pipeline
results cannot rescue a primary failure.

## Fixed catalogue and coverage

Reactome release 97 supplies assay-specific human sets from its all-level
Ensembl and UniProt mappings plus its pathway list. Exact URLs and SHA-256
fingerprints are in the registry. Reactome data are
[CC0](https://reactome.org/license), and the official download page provides
versioned mappings and quarterly archives ([Reactome downloads](https://reactome.org/download-data)).

RNA members are unique `ENSG` identifiers. Proteomics members are unique
canonical UniProt accessions. All human hierarchy levels and disease statuses
remain eligible; no outcome-dependent pathway class is removed. Declared sets
contain 10-200 members. A primary target additionally needs at least 70%
dataset coverage, at least ten matched members, and common observation of every
matched member in every primary unit. Coverage failures occur before scoring
and are reported with member identities.

## Methods and exact APIs

The ordered comparison is GeneFunnel, sum, mean, singscore, GSVA, and ssGSEA.
Sum and mean are the direct magnitude ablations. The three current methods span
sample-independent ranks, a within-sample random walk, and cohort-level CDF
normalization without expanding into a method zoo. This is motivated by the
large current single-sample benchmark showing substantial method-family,
missing-gene, set-size, and sample-size effects ([Toro-Domínguez et al.,
2025](https://academic.oup.com/bib/article/26/6/bbaf684/8382580)).

The protocol pins the Bioconductor 3.23 source archives, Git commits, SHA-256
hashes, constructor arguments, and serial execution calls for
[singscore 1.32.0](https://bioconductor.org/packages/singscore/) and
[GSVA 2.6.3](https://bioconductor.org/packages/GSVA/). GSVA's current
parameter-object API supplies both GSVA and ssGSEA; the old `method=` API is not
used. Singscore uses `tiesMethod="min"`, no stable-gene panel, and its
non-directional one-set score. GSVA uses an explicit Gaussian kernel for the
continuous matrices; ssGSEA uses `alpha=0.25` and normalization. Constant rows
have already been removed commonly, so GSVA row filtering is disabled to keep
membership identical.

PLAGE is excluded because its singular-vector sign is arbitrary, which is
incompatible with a prospectively oriented perturbation endpoint. The combined
z-score is excluded as a cohort-standardized mean already represented by the
mean ablation. ORA is group-level rather than a single-sample score. Direct-cell
methods belong to the separate cell-level hypothesis.

## Matched controls and task statistic

For each task, the data manifest fixes one Reactome stable ID, expected
activation (`+1`) or inhibition (`-1`), external mechanistic citation,
pair/block variables, and covariates. For every method and set, base R
`lm.fit(tol=1e-7)` estimates the condition coefficient from an identically
ordered design: intercept, condition, pair/block fixed effects, then manifest
covariates. The condition column must be estimable. Multiplication by the
expected direction orients larger values toward the expected response.

One hundred random feature sets are matched to the target separately within
each study. One hundred thousand seeded, exact-size, without-replacement
proposals exclude every target member. Matching uses only control-condition
units and six set summaries: mean and SD of member log-median abundance, mean
detection, median absolute within-set Spearman correlation, mean log catalogue
degree, and maximum Jaccard overlap with another retained Reactome pathway.
Proposal MADs scale the distances; all six coordinates must lie within one
scale unit. The lowest squared distances win with lexical ties. Fewer than 100
matches makes the task ineligible before any method score. This directly
controls size, abundance, detection, co-expression, and catalogue overlap/
degree without treating competitor agreement as truth.

For method `m`, task `t`, target `g`, and matched controls `c_j`, the endpoint is

\[
C_{mt}=\frac1{100}\sum_j
\{I(d_t\hat\beta_{mtg}>d_t\hat\beta_{mtc_j})+
\tfrac12I(d_t\hat\beta_{mtg}=d_t\hat\beta_{mtc_j})\}.
\]

It is the fraction of matched controls beaten by the target, with half credit
for ties. It ranges from zero to one; 0.5 is random target-control ordering.
It is a perturbation-ranking endpoint, not proof that a pathway is active.

## Co-primary held-out gate

Tasks receive equal weight within study and studies receive equal weight within
assay stratum. A 10,000-replicate hierarchical bootstrap resamples studies,
tasks within selected study, and biological units within task condition or
pair, then recomputes contrasts and matched-control statistics. Control-set
memberships stay fixed: resampling changes condition coefficients and ranks,
not the prospectively selected controls. All fifteen hypotheses - five for each
of bulk RNA, pseudobulk RNA, and proteomics - use one-sided Bonferroni
familywise alpha 0.05. Thus each lower bound is the type-8 bootstrap quantile at
`0.05/15`.

For each assay, every condition must hold:

1. GeneFunnel target rank lower bound is at least 0.55.
2. Sum and mean task ranks are identical: target and every control have the
   same member count, so the two score matrices differ only by one positive
   task-wide factor. Any observed disagreement is an execution error.
3. GeneFunnel minus that common sum/mean ablation has lower bound at least
   0.05.
4. GeneFunnel minus each of singscore, GSVA, and ssGSEA has lower bound at
   least -0.03.
5. Every held-out study has GeneFunnel point rank above 0.5, positive
   GeneFunnel-minus-ablation differences, and no current-method deficit below
   -0.03.

Five rank points means the target beats five additional matched controls per
hundred; three points is the maximum accepted loss to a current method. These
are prospective minimum material effects, not values reverse-engineered from a
result. Bulk RNA, pseudobulk RNA, and proteomics pass separately. A
stratum-specific claim requires that stratum's complete gate; a general claim
requires all three. Cross-stratum pooling, subgroups, native pipelines,
simulations, and secondary endpoints cannot rescue a failure.

## Secondary and null endpoints

Secondary outputs report method-native target rank; GeneFunnel magnitude,
balance, penalty, and effective-size contrasts; eligible pseudobulk aggregation
gaps; technical-repeat pathway-rank Spearman correlation; feature-thinning rank
loss; leave-one-unit/pair-out estimates; and five-run isolated elapsed/RSS
measurements. Feature thinning crosses 10%, 25%, and 50% removal with uniform,
low/high abundance, and low/high detection mechanisms. The earlier controlled
sensitivity gate failed, so these are failure-envelope descriptions, never a
reliability rescue. Named assay/endpoint families use BH at 0.05 and retain all
raw and adjusted values.

At least 20 independently declared null contrasts per assay compare
exchangeable control groups. The full model and covariates are the positive-task
model; the null is a zero condition coefficient. A Freedman-Lane procedure
permutes reduced-model residuals within manifest exchangeability blocks (or
swaps within pairs), using 9,999 draws or all distinct permutations. Two-sided
plus-one p-values receive BH correction within assay, study, and method across
eligible pathways. Because every tested contrast is known null, the report
retains discoveries/tests and any-discovery frequency with binomial intervals.
This evaluates a specified downstream pipeline; a GeneFunnel score alone has no
false-positive rate.

## Controlled simulation and seeds

The registry freezes the simulation factors before construction. RNA negative
binomial and proteomics log-normal measurement models cross target sizes
16/64/192, 4/8 units per condition, baseline balance 0.4/0.8, magnitude fold
1/1.25/2, balance change -0.2/0/0.2, and four set archetypes. Magnitude and
balance are fully crossed and set independently so the latent GeneFunnel score
is exactly `total * balance` before measurement.

Eleven nuisance factors cover depth, dropout, partial coverage, subunits,
catalogue overlap, contamination, outliers, complementary mixtures, member
dynamic range, variance, and correlation. A registered 128-run resolution-IV
fraction uses seven base signs and four three-factor products; the
implementation supplement must independently verify its resolution before a
value is generated. Two measurement repeats support latent-score RMSE,
target-control rank, repeat stability, and factor-stratified failure envelopes.
Simulation explains mechanism and failure; held-out external evidence remains
mandatory.

All RNG namespaces are fixed in the registry and reset using base R's
`Mersenne-Twister`/`Inversion`/`Rejection` combination. The implementation
supplement must map scenario/task indices to these namespaces without changing
them.

## Evidence and decisions

Every execution records clean source identity, package archive fingerprints,
installed versions, calls, seeds, environment, data manifests, exclusions, and
artifact hashes. Compact decisions/endpoints/failure envelopes are tracked, or
an immutable external result bundle receives a tracked hash manifest. Ignored
`benchmark/results/` files are working storage and never sole evidence.

Malformed provenance, schema, alignment, or arithmetic aborts. A scientific
gate failure is a completed result and retains its evidence. Passing only a
subset narrows the supported domain to those strata. If none passes, GeneFunnel
remains the faithful thesis implementation with its exact mathematical
diagnostics and no claim of general utility or superiority.
