<!-- Assisted-by: OpenAI Codex. -->

# Prospective validation amendment F-2.1.0

**Frozen before candidate-source selection, planned-object retrieval,
molecular-value access, comparative implementation, or results.** The machine
[`amendment registry`](validation-amendment.tsv) is authoritative. It applies
to the exact SHA-256-pinned `F-2.0.0` parent without modifying its bytes.
Amendment `replace` rows supersede the corresponding parent fields and `add`
rows extend them; that merged contract is the effective future design. This
document is explanation, not evidence or a result. No candidate source has
been selected and no candidate archive or molecular table has been accessed in
this work; that process record is not proof against external access.

Registry SHA-256:
`7bda0577d077a5204165bd570af6bd6eaeeb2d581571e17ff54e832b4f2cf2a5`.
Selection validator SHA-256:
`2caf5070d1d02180f3082d69d82026935384a52f99ec4a3fb100b72753e96888`.

## Why an amendment precedes source selection

Metadata-only feasibility review exposed six design defects:

1. pseudobulk tasks did not identify the source cell type that defines their
   aggregation and observations;
2. tasks sharing donors, controls, contrasts, or raw material were resampled as
   if independent;
3. a three-study, non-studentized `1/300` percentile bootstrap cannot justify a
   confidence-bound or familywise-error claim;
4. twenty supposedly independent null contrasts were generally infeasible,
   while 9,999 permutations cannot reach BH's first threshold when more than
   500 pathways are tested;
5. the proteomics matrix contract did not identify source quantity,
   normalization, inference, zero/missing, or archive-member semantics; and
6. the parent supplement mixed prospective source selection with facts available
   only after retrieval and overstated what timestamps and hashes can prove.

The amendment repairs those contracts prospectively. It deliberately narrows
the primary interpretation to a fixed held-out panel instead of manufacturing
population inference from purposively selected studies.

## Analysis contrast and task

An **analysis contrast** is one study, assay, biological population, reference
condition, intervention condition, exposure regimen, collection time, and
biological-unit frame. A positive **target task** adds one prespecified Reactome
stable ID, expected direction, and mechanistic citation. A held-out positive
study must contribute at least two distinct analysis contrasts; adding nested
targets to one contrast does not satisfy that requirement.

For pseudobulk RNA, population means exactly one source-supplied cell-type
field and one exact value fixed from metadata. Cell types are never pooled,
averaged, wildcard-selected, or chosen from molecular values. Raw counts are
summed within biological unit, condition, and that population. Every retained
unit-condition-population row needs at least 20 cells, and the contrast still
needs four independent units per arm or four complete pairs. Cells never add
replication.

Any tasks sharing a biological unit or raw material are dependent and carry
one transitive dependence-group ID. A unit-frame ID defines biological-unit
construction, not an exact realized observation set. F-S pool IDs plus exact
sample/cell predicates close planned membership; the post-parse manifest closes
the exact realized members. Identical source, unit frame, arm, and predicate
sets reuse one pool ID regardless of nuisance terms. Paired designs use joint
pools, unpaired arms use separate pools, and different pools may overlap. Every
post-parse row maps to one transitive dependence group and biological-cluster
ID so every pool and overlap is explicit. Source independence groups and task
dependence-group IDs are bijective: one raw-material group cannot be split
across apparently fresh dependence IDs or merged with another group.

Controls match population, exposure regimen, collection time, and every
declared factor except intervention. Constant provenance labels remain filters,
not aliased model terms. Every declared model has full rank, an estimable
condition coefficient, and positive residual degrees of freedom.

Every task-design row declares `field_type`. Unit fields are identifiers;
pair and block fields are categorical; covariates are categorical or numeric.
String model fields decode as non-null UTF-8 scalars without normalization or
coercion. Categorical levels sort by decoded bytes under the C locale, use the
first level as treatment reference, and ignore ambient contrast options.
Numeric covariates use finite, non-null, unscaled values under the F-E decimal
grammar. A malformed or missing model value makes the task visibly ineligible;
complete-case row omission, imputation, and ambient formula type inference are
forbidden.

Positive targets must occur in the pinned Reactome 97 eligible-target roster,
which is generated only from the three parent-hashed mappings. It contains 1,686
targets for each RNA assay and 1,690 for proteomics after the registered
assay-specific 10-200-member rule. The compact roster, generator, and both
SHA-256 identities are committed before task selection.

Execution cardinality is assessed after all rule-derived post-parse
eligibility. Every retained held-out study needs at least two distinct eligible
analysis contrasts; every assay needs at least three eligible independence
groups, six positive tasks, and two targets. Any shortfall makes that assay
incomplete without replacement or hidden exclusion.

## Proteomics input boundary

The primary matrix may start from documented raw or source-normalized linear
intensity, abundance, or area quantities. maxLFQ is allowed only when
source-normalized linear. All quantities must be unlogged, non-negative,
unimputed, and auditable. Centered, signed, ratio-only, log-only, batch-adjusted,
or imputed-only matrices are excluded. Source protein inference is allowed only
when documented. Human species filtering is fixed before values; shared-species
and multi-accession groups are excluded; canonical UniProt mapping and
zero/missing codes are retained.

One assay-input row per source closes the adapter contract. Value adapters are
`delimited_wide`, `delimited_long`, `tenx_h5`,
`maxquant_protein_groups`, or `fragpipe_combined_protein`; metadata and mapping
adapters are `tsv_header`, `csv_header`, `json_records`, or `parquet_table`.
Direct objects use `whole_object` with container `none` or `gzip`; archives use
one exact safe member with `zip`, `tar`, or `tar_gzip`. Matrix orientation is
`features_by_observations`, `observations_by_features`, or
`long_observation_feature`. Generic wide selectors are respectively
`exact_observation_columns` or `exact_observation_rows`; long tables use
`exact_quantity_field`, and 10x uses `tenx_matrix`. MaxQuant uses
`maxquant_intensity_columns` or `maxquant_lfq_intensity_columns`; FragPipe uses
`fragpipe_intensity_columns` or `fragpipe_maxlfq_intensity_columns`.
Observation selection is `exact_sample_ids_from_sample_metadata`, or for
pseudobulk `exact_cell_ids_from_cell_metadata` against its distinct cell
metadata. Selector fields are exact, distinct where required, and never
wildcarded.

ID, join, unit, pair, block, filter, and categorical-covariate values from every
metadata adapter must decode as UTF-8 string scalars. Typed numeric, Boolean,
binary, list, or null values are invalid in those fields. Quantity parsing
selects the registered member, axes, and field first; resolves exact missing
tokens or typed nulls; parses remaining C-locale decimal values; rejects
malformed, non-finite, negative, non-integer RNA, and nonzero-significand
underflow values; then applies the registered zero, repeat, duplicate, mapping,
and preprocessing rules. Missing tokens are a unique, non-empty, exact
pipe-delimited list; `<blank>` alone marks an empty input field, and no decimal
lexeme that parses to zero or non-finite may be a missing token.

A technical-repeat ID and aggregation operator occur together. RNA permits
only `sum_raw_counts`; proteomics permits `arithmetic_mean_linear` or documented
`source_nonoverlapping_linear_sum`; aggregation precedes biological-unit
scoring, and the repeat ID cannot alias sample, unit, pair, block, filter, or
covariate fields. RNA's duplicate-accession operator is `not_applicable`.
Proteomics uses `require_unique` or documented
`source_nonoverlapping_linear_sum`; the latter sums unlogged, nonoverlapping
linear quantities in first-occurrence order. Optional mapping object fields
form one closed member/container/adapter/source/target group, source and target
names are distinct, and the rule is `strip_ensembl_version` for RNA or
`canonical_uniprot_strip_isoform` for proteomics. Sample IDs are unique and
total over selected observations; cell IDs are unique and every cell resolves
to exactly one registered sample. Joins are one-to-one, or many cells to one
sample, never many-to-many.

Sample, cell, and mapping metadata are separately planned and allowlisted,
never aliases of sealed molecular bytes. For archive access, F-S validates
member and container syntax, not decoded contents. F-E archive adapters reject
duplicate member names, nested decoding, and every non-regular entry, including
symbolic links, hard links, devices, FIFOs, and sockets. They deliver only the
registered member without filesystem extraction or link following. Both
declared aggregate size and observed cumulative decoder output—including
listing, skipped members, and the target—must fit one per-object 256 GiB and
1000:1-to-F-B-byte budget. Further limits are 100,000 members, 64 MiB name
bytes, an 8 MiB stream buffer, 2 GiB address space, 64 descriptors, six CPU and
wall hours, and no child process. Direct gzip uses the same output, ratio, and
resource bounds. Counters are checked before allocation or output; truncation,
checksum, timeout, decoder, or resource failure rejects. An access class is an
evidence-backed prospective content contract, not proof that a file contains no
molecular values; a post-open mismatch remains ineligible and cannot be
reclassified, replaced, or used for held-out analysis.

The parent-pinned derived native matrix remains the only RNA native secondary.
Bulk-proteomics native secondary is unavailable in F-2.1.0 without a separately
prospective pinned object/member/adapter/quantity/join contract; it is never
discovered or reconstructed from the primary linear input and never replaces or
rescues primary same-input failure.

Retrieval-and-analysis permission is required. Redistribution permission is a
separate recorded fact and is never inferred from public access.

## Fixed-panel primary decision

The estimand is the rule-derived eligible subset of the finite, prospectively
selected candidate source/contrast/target/pipeline roster. Eligibility may
depend on registered QC, coverage, pathway support, and control-matching rules;
the decision is conditional on every realized choice. An independence group is
the study unit, and at most one source per assay may use it. Targets receive
equal weight within analysis contrast, contrasts equal weight within study, and
studies equal weight within assay. This prevents nested targets or duplicate
source records from silently increasing study weight.

Studies and deliberately selected tasks are fixed, not resampled. In each of
one million shared stability draws, every globally sorted dependence-group and
biological-cluster key receives one raw `Exp(1)` weight. Every observation pool
reuses its clusters' raw weights, normalizes across its unique member-cluster
keys to mean one, and propagates one weight to every row of that cluster across
dependent tasks, targets, controls, assays, and methods.
Paired observations use joint pools; unpaired arms use separate pools. This
synchronizes exact and partial donor overlap across interventions or data types.
Weighted-model rank or non-finite failures fail the assay without redraw. Score
matrices and realized matched-control sets stay fixed, and method differences
use the same draw.

The type-8 `0.05/15` lower quantile and original effect thresholds define a
stringent nominal fixed-panel stability decision. It is not a confidence bound
and does not provide familywise-error or population-coverage control. The
every-study rule is a deterministic consistency filter, not statistical or
target-level replication. A pass supports only the registered panel. Broader
study-, target-, or assay-population inference needs a later protocol with a
defined sampling frame and validated cluster-level inference.

Every method must provide the identical rule-derived eligible task Cartesian
product. A missing, non-finite, or failed GeneFunnel or comparator cell fails
the entire assay comparison gate without changing any denominator.

## Negative-control diagnostic

A null contrast uses either a value-blind artificial two-level label on
independent biological units receiving the same source condition and handling,
or source-declared paired exchangeable control labels. Artificial assignment is
balanced within declared exchangeability blocks; paired labels swap within
complete pairs. Every arm has at least four biological units. Cells and
technical repeats are never assignment units. The effective design requires at
least 20 execution-eligible contrasts per assay across at least three held-out
independence groups and at least four contrasts from every contributing group.
Rank, estimability, unit, distinct-partition, or permutation-resolution
failures remain visible and do not count or trigger replacement. Reused units
keep one dependence group; multiple labelings are contrasts, not independent
biological evidence. For artificial assignments, unit and declared block IDs
sort under the `C` locale; absence of a block field means one global block. The
feasible set is all binary label vectors with per-block and global arm-count
differences at most one and at least four units per arm. Only the lexically
smaller member of each label-complement pair remains; canonical vectors sort
lexically, and the registered label seed selects one uniform integer rank
through the exact `F-E-1.0.0` rejection sampler. Duplicate realized canonical
partitions reject.

The full null model is intercept, artificial condition, then declared varying
pair/block effects and covariates; the reduced model omits artificial condition
only. One Freedman-Lane permutation schedule is shared by all pathways and
methods in a contrast. Its unique transformation encodings sort under `C` with
identity first; sampled schedules contain distinct nonidentity transformations
and use the separate permutation seed. BH at `q=0.05` is applied separately to
each null contrast-method family; the two-sided statistic is the absolute
studentized condition coefficient using the `F-E-1.0.0`-pinned `lm` QR,
residual-df, and standard-error computation.

For `K` tested pathways, sampled nonidentity permutations number at least

```text
max(9999, ceiling(K / 0.05) - 1).
```

`M` counts unique allowed transformations including identity. If `M-1 < B`, all
`M` transformations are enumerated and the full-orbit fraction has minimum
`1/M`. Otherwise `B` distinct nonidentity indices are sampled without
replacement and the plus-one p-value has minimum `1/(B+1)`. Resolution
eligibility is value-independent and requires the applicable minimum to be at
most `0.05/K`; it does not assert that statistic ties permit rejection. Fewer
than 20 resolution-, rank-, unit-, partition-, and execution-eligible contrasts
or fewer than three independence groups makes the assay validation incomplete
without altering its retained primary number.

For each contrast and method, report rejections/tests and whether any pathway
was discovered. Average contrasts within study, then studies equally; retain
all values and the leave-one-study-out range. These dependent conditional
stress-test summaries are not FDP/FDR estimates, receive no binomial or
biological confidence interval, never validate a score's false-positive rate,
and cannot promote or rescue the primary gate.

## Prospective artifact order

The selection, byte, execution-closure, and access contracts have separate
identities:

1. Commit metadata-only `F-S-1.0.0`: exact sources, tasks, design, resampling,
   filters, assay inputs, planned objects, access roles, licenses, and candidate
   order. Every row is selected; order is global execution order, never fallback.
2. At clean selection commit S, record acquisition start and stream every
   planned object in order without archive listing, decompression, parsing,
   preview, content-derived type detection, or other content-aware access. Hash
   and stat the stored opaque bytes.
3. At byte commit B, commit `F-B-1.0.0`: S and its F-S hash plus one terminal
   row per planned object, with exact locator, retrieval UTC, byte count,
   SHA-256, or enumerated failure. `verified` attests that the acquisition client
   observed a successful full-body terminal HTTP 200 and that stat-hash-stat
   verification matched the exact closed object root; the validator does not
   independently prove server completeness. Partial 206 and empty 204 responses
   are ineligible. F-B contains no sample mapping or preprocessing claim.
4. At B, record the first access to verified allowlisted metadata and development
   molecular objects in object order, then freeze complete implementation commit
   I with its committed `F-E-1.0.0` closure before held-out molecular access.
   Content-sealed means a declared guarded process state between opaque
   acquisition and I, not cryptographic sealing, OS isolation, or proof of
   non-access.
5. At I, re-hash every verified held-out molecular object as a preopen check,
   record exactly one execution start before enabling any held-out decoder, then
   record each held-out first content access in object order and execute the
   panel once. Sample and preprocessing mappings belong to a later post-parse
   manifest.

The future F-B transfer runner uses GET over HTTPS port 443 with TLS peer and
hostname verification on every hop, `Accept-Encoding: identity`, and no
transparent decoding. It strips proxy, netrc, cookie, credential, and auth
state; permits at most ten redirects, 8 MiB headers, 30 seconds to connect, 300
seconds read-idle, 48 hours wall time, 1 TiB per object, 4 TiB per attempt, one
transfer at a time, and a 4,100 GiB disk quota. Each vetted public address is
used for connection and the actual peer is rechecked; credentialed, rebound,
private, loopback, link-local, multicast, and reserved targets reject. The
stream counter, not `Content-Length`, is authoritative. Exclusive no-follow
temporary files, cleanup, fsync, stat-hash-stat, atomic rename, and isolated
network/filesystem/process policy complete acquisition. Locator syntax alone is
not SSRF protection.

`F-A-1.0.0` is one contiguous acquisition-through-held-out-first-access
transcript: acquisition start; every planned object's stream outcome; F-B
freeze; every allowlisted-metadata or development-content first access;
implementation freeze; every held-out preopen verification; one execution
start; and every held-out first content access. It has one attempt ID, no other
event or retry, and exact object byte/hash facts. Its
`execution_closure_sha256` is `none` before implementation freeze and thereafter
pins the exact committed `benchmark/validation-execution-closure.R` from I.

Production F-A validation requires `verify_git = TRUE`, the exact object root,
and live clean `HEAD == I`; S is a strict ancestor of B and B of I. The exact
closed F-S tree object and all registered semantic-runtime paths are byte- and
mode-identical regular blobs at S/B/I. F-B is the same regular blob at B/I. The
final F-E hash is absent at its canonical path in S/B and is a regular blob at
I. Shallow history, replace refs, grafts, alternates, indirect/worktree unsafe
Git config, or a changed runtime rejects. Evidence commands use absolute Git in
an empty fixed environment, so inherited environment and Git state are ignored.
Clean additionally means zero C-locale porcelain output including untracked
files and submodules.

Production loads only the logged F-E file from I into an environment whose
parent is `baseenv()`, verifies exact three required formals without defaults,
and accepts only the identical plain list `schema_version = F-E-1.0.0`,
`implementation_commit = I`, `validator_sha256 = <logged hash>`, and
`executable = TRUE`. It derives an object's only runtime path as
`object_root/access_role/object_id.bin`, opens it read-only and no-follow,
fstat-hash-fstats the same descriptor against F-B, and passes that descriptor
unchanged to its registered adapter. Reopen, alternate paths, fallback, and
post-guard substitution are forbidden; a staging copy must be byte-identical.
F-A rechecks the complete object and artifact closures after F-E validation.
Missing, modified, early-final, or nonconforming closure fails closed; until
F-E and the controlled runner exist, the design remains non-executable.

The on-disk hash and fresh base environment do not sandbox trusted project code
or protect against an already-modified calling process. Structural F-S/F-B/F-A
guarantees assume a trusted process that freshly sources the exact semantic
runtime. The future controlled runner uses a sanitized vanilla subprocess and
integrated guarded descriptors; a fresh environment supplies name-lookup
isolation only, not filesystem, network, or process isolation.

Structural F-A validation rejects omissions, reordering, mixed attempts, wrong
object facts, and invalid state transitions, but closes only the supplied
transcript. Git ancestry, clean-state records, access events, content-seal
guards, and future controlled-runner append/open guards provide auditable
ordering evidence; none proves absence of external or unlogged access.
Timestamps alone are insufficient, and no structural check excludes an
alternate, retrospective, or retried transcript. The future runner therefore
needs a unique external append-only attempt anchor.

## Structural validation posture

The hermetic selection fixture exercises 12 sources, 90 tasks including 60
null contrasts, 29 planned objects, and 37 F-A events. Its current 99
adversaries cover schema, UTF-8 scalar and numeric parsing, dependence and pool
identity, archive/resource limits, object closure, guarded descriptors, Git
tree/runtime/mode evidence, and access-state transitions. These are structural
contract tests, not selected sources, molecular evidence, or a comparative
result. F-2.1.0 comprises 128 amendment rows and 266 effective rows over the
immutable 184-row parent. The lineage remains design-only until exact F-S,
F-B, and F-E artifacts and the controlled runner exist.
