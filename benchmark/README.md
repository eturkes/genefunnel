<!-- Assisted-by: OpenAI Codex. -->

# GeneFunnel benchmark protocol

The tracked runners plus [`protocol.tsv`](protocol.tsv) define benchmark
protocol 1.0.0. The audit manifest records scenario seeds, dimensions, matrix
and gene-set construction, preprocessing, methods, assertions, and environment
artifacts; each execution also writes its actual expanded settings. A change to
fixture meaning, method, or assertion requires a protocol version change in
both `protocol.tsv` and `protocol.R`.

All fixtures are synthetic, deterministic, generated without network access,
and scored from an installed package in fresh R processes. Generated evidence
belongs under `benchmark/results/` and stays untracked.

The separate [`components-protocol.md`](components-protocol.md) and
[`components-protocol.tsv`](components-protocol.tsv) freeze A2's numerical and
default-path performance gates without changing scientific protocol 1.0.0.

Workstream B's separate
[`aggregation protocol B-1.0.3`](aggregation-protocol.md), machine
[`registry`](aggregation-protocol.tsv), and exact external
[`data manifest`](aggregation-data.tsv) freeze the internal audit schema,
controlled factorial, known-mixture validation, donor split, endpoints, and
go/fallback rules before prototype results. Its external data remain ignored;
the source manifest pins HTTPS bytes and adds no package dependency.

## Install

```sh
R_LIBS_USER="$PWD/.agent/R-library" R CMD INSTALL \
  --library="$PWD/.agent/R-library" .
```

## Controlled scientific protocol

```sh
R_LIBS_USER="$PWD/.agent/R-library" Rscript --vanilla \
  benchmark/run-controlled.R
```

The runner checks analytic or independently paired expectations for:

- equal member values;
- equal sums with low versus high within-set deviation;
- all-zero and one-nonzero sets;
- sample-specific missingness and effective set size;
- partial matrix coverage;
- a seeded 3%-stored single-cell-like sparse matrix versus its dense
  equivalent;
- sample and gene-set independence.

No preprocessing is performed: values are constructed directly on a
non-negative scale with a meaningful zero; `NA`/`NaN` cells remain missing.
The runner emits assertion-level `results.tsv`, aggregate `summary.tsv`, and a
plain Markdown `report.md`, then fails if any assertion fails.

## Computational performance protocol

Select a preset and isolated repeat count:

```sh
R_LIBS_USER="$PWD/.agent/R-library" Rscript --vanilla \
  benchmark/run.R --preset=full --repeats=3 --workers=2
```

Presets:

- `full` - 20,000 features; dense bulk-like matrices with 60 samples and
  sparse single-cell-like matrices with 600 samples; 1,000 sets of size 20;
  low/high overlap; serial/two-worker SOCK execution.
- `hotspot` - 20,000 x 200 sparse matrices across stored densities 0.1%, 1%,
  3%, and 10%, with both overlap patterns; isolates sparse traversal and
  catalogue preparation.
- `smoke` - the same execution paths at tiny dimensions. `ci-smoke.R` checks
  this preset plus the controlled protocol without time or memory thresholds.

`--output=DIR` selects a result directory. `--workers=N` changes the SOCK
worker count; serial cases use one worker.

Dense values are exponential, then exact seeded positions become zero and
sample-specific `NA`/`NaN`. Sparse matrices are constructed directly from
unique coordinates: unstored cells are implicit zeros and a seeded fraction
of stored cells becomes `NA`/`NaN`. Sparse fixtures are assembled from stored
coordinates without a logical dense intermediate. Low-overlap sets share no
members; high-overlap sets share 75% of their members. Serial and parallel
variants reuse fixture seeds and must produce identical serialized score
digests.

## Component default-path comparison

Prepare exact locked-baseline/current-candidate libraries, backed by the same
dependency library:

```sh
R_LIBS_USER="$PWD/.agent/R-library" Rscript --vanilla \
  benchmark/prepare-components.R \
  --output="$PWD/.agent/tmp/components-comparison"
```

The preparer archives exact Git trees, installs them separately, and records
source plus installed-tree fingerprints. Exercise orchestration with:

```sh
R_LIBS_USER="$PWD/.agent/R-library" Rscript --vanilla \
  benchmark/run-components.R \
  --baseline-library="$PWD/.agent/tmp/components-comparison/baseline-library" \
  --candidate-library="$PWD/.agent/tmp/components-comparison/candidate-library" \
  --mode=smoke
```

Gate mode additionally requires the locked full baseline SHA and the full
candidate SHA. Its 30 paired/interleaved repeats, workloads, one-sided timing
interval, allocation/RSS limits, and rejection policy are fixed in the
component protocol. Smoke mode uses reduced matrices and four repeats; it
checks orchestration and identity only and cannot support a performance claim.
Gate mode is Linux-only and aborts when recorded one-minute load per logical
CPU exceeds the protocol's quiescence limit.

## Aggregation synthetic validation

Run the complete frozen B-1.0.3 design only from a clean committed tree. Its
synthetic fields are identical to the completed B-1.0.2 design:

```sh
R_LIBS_USER="$PWD/.agent/R-library" Rscript --vanilla \
  benchmark/run-aggregation-synthetic.R --workers=4
```

The runner installs that exact tree into a temporary isolated library, checks
four deterministic generator/audit strata, then evaluates 62,208 shared latent
scenarios with independent A/B measurements. Scenario-local seeds make one- or
multi-worker results identical. `--chunk-size=N` controls atomic RDS
checkpoints; rerunning with the same `--output=DIR` safely resumes only when the
protocol, Git commit, chunk identities, and manifest match.

Outputs include every observation, co-primary endpoint, factor stratum, paired
fold prediction, QR coefficient/alias fact, 2,000 bootstrap estimates, compact
summary/report, and an artifact SHA-256 manifest. A failed scientific gate is a
valid completed result: the runner records `FAIL` and exits successfully. A
schema, provenance, installation, arithmetic, or orchestration failure aborts.
Generated evidence remains ignored under `benchmark/results/`.

The complete B-1.0.2 execution from clean commit `d1cbf15` passed every frozen
synthetic gate. Its compact tracked
[`result`](aggregation-synthetic-result.md) and
[`endpoint table`](aggregation-synthetic-result.tsv) also retain the severe-
dropout failure envelope: controlled success does not establish dropout
robustness or biological validity.

## CellBench known-mixture validation

After placing the four exact `RNAmix_celseq2.*` and `RNAmix_sortseq.*` files
from [`aggregation-data.tsv`](aggregation-data.tsv) in one external directory,
run only from a clean committed tree:

```sh
R_LIBS_USER="$PWD/.agent/R-library" Rscript --vanilla \
  benchmark/run-aggregation-cellbench.R --data-dir="$PWD/.agent/tmp"
```

The runner verifies all four payload hashes before parsing, installs the exact
commit into an isolated library, derives 12 sets only from CEL-seq2 pure
profiles, and retains every library/set error, all 384 fixed error groups, 384
condition/set medians, and five co-primary endpoint decisions. Missing fixed
groups or undefined metrics produce a retained scientific `FAIL`; malformed or
misaligned input aborts. A scientific failure is a completed result and exits
successfully. Generated output remains ignored under `benchmark/results/`.

The complete B-1.0.3 execution from clean commit `9b3cab6` failed both frozen
CellBench gates. The compact tracked
[`negative result`](aggregation-cellbench-result.md) and
[`endpoint table`](aggregation-cellbench-result.tsv) retain the zero-score,
pair-set error, and incomplete-condition failure envelope. This closes the
protocol's public-promotion route; complete-case correlations cannot rescue the
fixed-grid failure.

## Kang donor-perturbation characterization

After placing the exact `GSE96583_RAW.tar`, batch-2 gene/metadata files, and
Reactome v97 archive from [`aggregation-data.tsv`](aggregation-data.tsv) in one
external directory, run only from a clean committed tree:

```sh
R_LIBS_USER="$PWD/.agent/R-library" Rscript --vanilla \
  benchmark/run-aggregation-kang.R --data-dir="$PWD/.agent/tmp"
```

The runner verifies all four payload hashes, enforces the exact Matrix Market
and duplicated-barcode join contracts, installs the committed package in an
isolated library, and constructs full/odd/even raw-UMI cell-type profiles. It
retains the fixed 96-unit eligibility grid, every included Reactome pathway
audit, all 16 donor-condition correlations, 16 primary-pathway donor contrasts,
two held-out/sign-test decisions, and six endpoints. Undefined fixed-grid
metrics produce a retained scientific `FAIL`; malformed input aborts. This
characterization cannot rescue the already-failed CellBench promotion gate.

The complete B-1.0.3 execution from clean commit `62d42dc` passed technical
split stability but failed the combined held-out and biological-effect gates.
The compact tracked [`result`](aggregation-kang-result.md) and
[`endpoint table`](aggregation-kang-result.tsv) retain unit/pathway support,
the 16 correlations, donor directions, and exact sign-test failure. Workstream
B's ordered external validation is complete; no public audit or perturbation
claim follows.

## Observed-member sensitivity

Protocol [`E-1.0.0`](sensitivity-protocol.md) fixes the internal compact schema,
exact dyadic-rational brute oracle, profile-before-optimization rule, controlled
feature masks, measurement repeats, held-out folds, models, bootstrap, and
incremental prediction gates; it was frozen before a package sensitivity value
existed. The unexported exact brute prototype now implements the compact schema,
while the empirical reliability gates remain pending. The tracked
[`fixed profile`](sensitivity-profile-result.md) records median elapsed 213.585 s
and exact-stack share 0.999403, so both frozen optimization triggers pass. This
permits implementation research only. The tracked
[`sorted-prefix adoption check`](sensitivity-optimization-result.md) then
reproduces the fixed brute output digest from a clean isolated candidate. Check
the registry's exact grid without generating a diagnostic:

```sh
Rscript --vanilla -e \
  'source("benchmark/sensitivity-protocol.R"); \
   print(sensitivity_validate_protocol("."))'
```

The parent profile dimensions/seeds omitted exact fixture and measurement-pass
construction. Byte-pinned execution supplement
[`E-P-1.0.0`](sensitivity-profile-protocol.tsv) closes that gap after the brute
implementation but before any fixed profile call. Validate both identities:

```sh
Rscript --vanilla -e \
  'source("benchmark/sensitivity-protocol.R"); \
   source("benchmark/sensitivity-profile.R"); \
   print(sensitivity_validate_protocol(".")); \
   print(sensitivity_profile_validate("."))'
```

After committing the runner, execute the clean-SHA profile with:

```sh
candidate_id=$(git rev-parse HEAD)
Rscript --vanilla benchmark/run-sensitivity-profile.R \
  --mode=gate --candidate-id="$candidate_id" \
  --output=benchmark/results/sensitivity-profile
```

Gate mode installs the exact Git archive in an isolated library, records three
timed calls plus separate CPU/allocation profiles, requires one output digest,
and applies only the frozen optimization-eligibility rule. Its timing is
descriptive and cannot promote a reliability API.

The registry contains 5,760 latent scenarios, 345,600 feature-loss rows, 5,760
controlled-repeat rows, and exactly 576 scenarios per held-out fold. This first
stage is synthetic and internal: even a complete pass cannot support export or
a claim about real technical/biological replicates.

The parent controlled registry deliberately remains unchanged. Its terms did
not uniquely choose R's Bernoulli call, predictor dummy/scaling details, or
bootstrap draw order, so byte-pinned execution supplement
[`E-C-1.0.0`](sensitivity-controlled-protocol.tsv) closes those degrees before
any controlled scenario is constructed. Validate it with:

```sh
Rscript --vanilla -e \
  'source("benchmark/sensitivity-protocol.R"); \
   source("benchmark/sensitivity-controlled-protocol.R"); \
   print(sensitivity_controlled_validate_protocol("."))'
```

`sensitivity-controlled.R` implements the registered profile, Poisson/NB,
post-count Bernoulli dropout, weighted masks, paired absence encodings, package
scores, partial-input exact diagnostics, targets, and fixed-schema observation
rows. CI exercises four size/archetype/noise strata twice, explicitly checks
row-absence/`NA` identity, and rejects a mutated target.

## Generated artifacts

Every runner writes:

- `protocol.tsv` - a copy of the tracked audit manifest;
- `manifest.tsv` - actual scenarios, seeds, dimensions, and settings;
- `metadata.tsv` and `session-info.txt` - Git state, OS, hardware, R, package
  versions, and execution settings;
- `report.md` - human-readable methods, results, environment, limitations, and
  reproduction command.

The performance runner additionally writes `runs.tsv`, `summary.tsv`, and a
`runs/` directory containing each isolated fixture manifest, allocation trace,
stdout/stderr, and external resource measurement. The controlled runner writes
`results.tsv` and `summary.tsv`.

The component comparison additionally writes `packages.tsv`, paired
`pairs.tsv`, `load.tsv`, and `decision.tsv`. Every worker verifies that
`genefunnel` loaded from its assigned baseline/candidate library. Unlike the
general performance runner below, the component timer uses no concurrent RSS
sampler; passive GNU-time process RSS is its resource gate.

The sensitivity profile writes `manifest.tsv`, three-call `runs.tsv`,
`profile.tsv`, `allocation.tsv`, `decision.tsv`, raw `Rprof`/`Rprofmem` traces,
an isolated source/library record, and artifact hashes. Its decision concerns
optimization eligibility only.

The preparer writes `installations.tsv`, exact source archives, install logs,
per-library provenance markers, and installed-package tree fingerprints. Gate
mode recomputes the source and installed fingerprints before any worker runs.

The aggregation runner additionally writes `observations.tsv`, `endpoints.tsv`,
`fold-results.tsv`, `predictions.tsv`, `model-coefficients.tsv`,
`bootstrap.tsv`, `strata.tsv`, `summary.tsv`, `installation.log`, atomic
`checkpoints/`, and `artifacts.tsv`.

The CellBench runner writes verified `data-files.tsv`, `set-manifest.tsv`,
`set-membership.tsv`, `pure-profiles.tsv`, `references.tsv`,
`observations.tsv`, `curve-groups.tsv`, `condition-medians.tsv`,
`endpoints.tsv`, `summary.tsv`, `installation.log`, and `artifacts.tsv`.

The Kang runner writes verified `data-files.tsv`, sparse-matrix and
preprocessing facts, unit/cell assignments, pathway manifest/membership,
`audit-summary.tsv`, `stability.tsv`, `donor-contrasts.tsv`,
`pathway-decisions.tsv`, six endpoints, environment/report files, and artifact
hashes.

## Compiled catalogue comparison

Protocol [`C-1.0.2`](catalogue-protocol.md) tests the internal compiled-catalogue
spike without promoting an API. Prepare the exact named-list parent and current
candidate from Git snapshots:

```sh
R_LIBS_USER="$PWD/.agent/R-library" Rscript --vanilla \
  benchmark/prepare-catalogue.R \
  --output="$PWD/.agent/tmp/catalogue-comparison"
```

Exercise all isolated timing/resource, digest, parser/tamper, and fresh/reused
SOCK paths with downscaled fixtures:

```sh
R_LIBS_USER="$PWD/.agent/R-library" Rscript --vanilla \
  benchmark/run-catalogue.R \
  --baseline-library="$PWD/.agent/tmp/catalogue-comparison/list-library" \
  --candidate-library="$PWD/.agent/tmp/catalogue-comparison/compiled-library" \
  --mode=smoke
```

Smoke uses four balanced pairs, two batches, and smaller matrices. It must pass
identity/orchestration but records `performance_claim = FALSE` and applies no
timing/resource threshold. Gate mode requires a clean repository at the full
candidate SHA plus baseline
`a573c124909235e41bdbc3cfae950947465d8755`; it verifies both snapshot markers,
aborts above load/core `0.25`, and executes the fixed 20-pair/five-batch decision:

```sh
candidate_id="$(git rev-parse HEAD)"
R_LIBS_USER="$PWD/.agent/R-library" Rscript --vanilla \
  benchmark/run-catalogue.R \
  --baseline-library="$PWD/.agent/tmp/catalogue-comparison/list-library" \
  --candidate-library="$PWD/.agent/tmp/catalogue-comparison/compiled-library" \
  --baseline-id=a573c124909235e41bdbc3cfae950947465d8755 \
  --candidate-id="$candidate_id" \
  --mode=gate
```

The runner writes tracked-protocol and expanded manifests, exact package
records, correctness/load/session metadata, isolated run logs, raw/paired/
summary tables, a report, and a machine decision. Generated evidence remains
under ignored `benchmark/results/` or `.agent/tmp/` paths. Each catalogue
observation uses separate fresh timing, allocation, and RSS processes. GNU time
records the timing process's whole-process maximum; a ready-synchronized,
resource-only `/proc` sampler records compilation and scoring peak increments.
No concurrent sampler perturbs primary timing.

For the general computational runner, `Rprofmem()` reports cumulative
manager-side R allocations, not retained or native allocation. On Linux, a 10
ms `/proc` sampler records aggregate resident
memory for the scoring process tree while excluding itself. The baseline
includes the input fixture; the increment includes output, bounded work
buffers, backend startup, and workers. Compatible GNU `time` records maximum
RSS separately.

Elapsed time covers `genefunnel()` only: validation, set preparation, backend
startup, chunk transfer, and scoring. It excludes package loading and fixture
generation. Observations are hardware/software-specific evidence, never
universal thresholds.

## Provenance boundary

The repository has no thesis datasets, so this protocol neither reruns nor
claims reproduction of thesis results or historical runtimes. It includes no
competitor methods because the package makes no comparative claim and adding
their dependency/API surface would not strengthen a current acceptance gate.

If real thesis data become available, keep them external. Optional retrieval
and execution scripts must record source URLs, licensing, checksums, immutable
dataset identifiers, preprocessing, feature and gene-set versions, and the
full environment; results must remain clearly separated from this synthetic
protocol.
