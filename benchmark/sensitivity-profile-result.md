<!-- Assisted-by: OpenAI Codex. -->

# Exact sensitivity profile result

**Decision:** optimization work is eligible. This is a descriptive engineering
decision, not a speed, reliability, biological, or public-API claim.

## Execution identity

- Parent protocol: `E-1.0.0`; profile supplement: `E-P-1.0.0`.
- Clean candidate: `37120c125efbdaf3e2c3b30e8049465c4874460b`.
- Source archive SHA-256:
  `714bac6613456dcc29940be2940821949248d6ee7bb54f10754510fdb3e41082`.
- Installed-package manifest MD5: `d85cf5751d014f656f2b980ccffbb371`.
- Fixture: 4,096 features x 12 samples; 48 sets x 128 members; matrix seed
  28072027; set seed 28072028; fixture MD5
  `32fa8d4843d024765a4adfc676793dbe`.
- Backend: `BiocParallel::SerialParam(progressbar = FALSE)`.

E-1.0.0 fixed the dimensions, seeds, metrics, and thresholds but omitted the
value/member constructors and exact pass boundaries. E-P-1.0.0 recorded that
omission and byte-pinned those mechanics after brute implementation but before
this fixed call. It changed no correctness or empirical threshold.

## Timed calls

| Repeat | Elapsed (s) | User (s) | System (s) |
|---:|---:|---:|---:|
| 1 | 239.364 | 206.192 | 0.189 |
| 2 | 197.017 | 186.536 | 0.036 |
| 3 | 213.585 | 201.542 | 0.138 |

Median elapsed = **213.585 s**, above the frozen strict 60 s trigger. Every
timed, CPU-profile, and allocation-profile output had MD5
`3d9635e779a9ed1eee453a2a04596369`.

## Attribution and allocation

The separate 10 ms `Rprof()` pass took 312.829 s. Of 23,446 sampled stacks,
23,432 contained a `.gf_` exact-arithmetic frame: **0.999402883221019**,
above the frozen strict 0.50 trigger.

The separate `Rprofmem(threshold = 0)` pass took 306.724 s and recorded
607,529 numeric allocation events totaling 378,561,864 manager-side R bytes.
Those cumulative allocations and profiled-pass runtimes are diagnostic context,
not primary timing or retained-memory measurements.

Both trigger branches pass, so `optimization_eligible = TRUE` and
`performance_claim = FALSE`. Any optimized implementation must remain internal
and match the brute result and status on every fixed and randomized cell before
adoption. The brute path remains the correctness oracle.

## Environment and evidence

Execution used R 4.6.1, Matrix 1.7.5, BiocParallel 1.46.0, an Intel Core Ultra
7 268V host with eight logical cores, and 32,329,932 KiB reported memory. Raw
ignored artifacts are authenticated by:

- `Rprof.out` SHA-256:
  `b44850aa8a5c84fbcde61461e0fcb6363ede5c831ab1f8c29b7ec066db536d8a`;
- `Rprofmem.out` SHA-256:
  `348560b236c7185bf7ff84ae8f3bc8e5394462d7baf5335595ee4e1d12a74010`;
- core artifact manifest SHA-256:
  `15a2aa727675df4b3f483f9db420d02a475a4666c162ea7b80de49318f513ebf`.

The machine-readable companion is
[`sensitivity-profile-result.tsv`](sensitivity-profile-result.tsv). Reproduce
from the candidate commit with a new output directory:

```sh
Rscript --vanilla benchmark/run-sensitivity-profile.R \
  --mode=gate \
  --candidate-id=37120c125efbdaf3e2c3b30e8049465c4874460b \
  --output=benchmark/results/sensitivity-profile-reproduction
```
