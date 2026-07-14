# GeneFunnel benchmarks

The harness runs installed-package scoring in a fresh R process for each
scenario and repeat. Fixtures are deterministic and generated without network
access. Generated results belong under `benchmark/results/` and stay untracked.

## Run

Install the current source first, then select a preset:

```sh
R_LIBS_USER="$PWD/.agent/R-library" R CMD INSTALL \
  --library="$PWD/.agent/R-library" .

R_LIBS_USER="$PWD/.agent/R-library" Rscript --vanilla \
  benchmark/run.R --preset=full --repeats=3 --workers=2
```

Presets:

- `full` - 20,000 features; dense bulk-like matrices with 60 samples and sparse
  single-cell-like matrices with 600 samples; 1,000 sets of size 20; low- and
  high-overlap catalogues; serial and two-worker SOCK execution.
- `hotspot` - 20,000 × 200 sparse matrices across stored densities 0.1%, 1%,
  3%, and 10%, with both overlap patterns. This isolates sparse traversal and
  catalogue preparation scaling.
- `smoke` - the same execution paths at tiny dimensions. `ci-smoke.R` verifies
  runner completion, metrics, and backend-identical output digests without
  enforcing time or memory thresholds.

`--output=DIR` selects a result directory. `--workers=N` changes the SOCK worker
count; serial cases always use one worker.

## Fixtures

Dense values are exponential, then exact seeded positions become zero and
sample-specific `NA`/`NaN`. Sparse matrices are constructed directly from
unique coordinate samples: unstored cells are implicit zeros and a seeded
fraction of stored cells become `NA`/`NaN`. Generation never creates a dense
copy of a sparse fixture.

Low-overlap sets partition sampled features and share no members. High-overlap
sets share 75% of their members. Serial and parallel variants reuse the same
matrix and catalogue seeds; the runner fails if their serialized score digests
differ.

## Outputs

- `manifest.tsv` - seeds, dimensions, density, missingness, overlap, and backend
  for every scenario.
- `runs.tsv` - elapsed/user/system time, score and member-cell throughput,
  manager-side R allocations, memory, object sizes, output digest, and package
  versions for every isolated run.
- `summary.tsv` - per-scenario minimum/median timing and maximum observed memory.
- `metadata.tsv` and `session-info.txt` - Git state, OS, CPU, installed memory,
  R, package versions, and session details.
- `runs/` - per-run fixture manifest, allocation trace, stdout/stderr, and
  external resource measurement.

[`Rprofmem()`](https://stat.ethz.ch/R-manual/R-devel/library/utils/html/Rprofmem.html)
reports cumulative manager-side R allocations, not retained or native
allocation. Allocation tracing runs as a separate scoring pass and its output
must be identical to the timed pass. On Linux, a 10 ms `/proc` sampler records aggregate resident
memory for the scoring process tree while excluding the sampler itself. The
baseline includes the input fixture; the increment includes scoring output,
bounded work buffers, backend startup, and worker processes. GNU `time` maximum
RSS is also recorded when a compatible
[`/usr/bin/time`](https://www.gnu.org/software/time/manual/time.html) exists.
These measures answer different questions and should remain separate.

Elapsed time covers `genefunnel()` only, including validation, set preparation,
backend startup, chunk transfer, and scoring. It excludes package loading and
fixture generation. Treat observations as hardware/software-specific evidence,
not universal performance promises.
