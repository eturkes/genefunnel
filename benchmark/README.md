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

`Rprofmem()` reports cumulative manager-side R allocations, not retained or
native allocation. On Linux, a 10 ms `/proc` sampler records aggregate resident
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
