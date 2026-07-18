<!-- Assisted-by: OpenAI Codex. -->

# Compiled catalogue spike protocol C-1.0.2

**Frozen before implementation results.** Runtime control = named-list
`genefunnel()` at the exact spike parent; candidate = internal compiled path at
the eventual result commit. [`catalogue-protocol.tsv`](catalogue-protocol.tsv)
is the machine manifest. Any change to fixture meaning, endpoint, threshold,
repeat schedule, or multiplicity rule requires a new protocol version recorded
before inspecting replacement results.

`C-1.0.1` is a pre-result feasibility correction: `C-1.0.0` requested 50,000
disjoint low-overlap memberships from 20,000 features. The feature universe is
now 50,000, exactly matching that membership count; endpoints, thresholds,
seeds, orders, and every other fixture dimension remain fixed. No implementation
or performance result was inspected when this correction was committed.

`C-1.0.2` is a pre-run manifest-completeness correction: the already described
zero/missingness and duplicate/unmatched fixtures now have explicit TSV columns
for their fractions/strides. Values are dense zero `0.35`, dense cell missing
`0.01`, sparse stored missing `0.05`, unmatched stride `20`, and duplicate
stride `25`. All other protocol state remains unchanged; no runner or protocol
result existed when this correction was committed.

## Question and unit

Does one compilation plus five scoring calls safely reduce cumulative elapsed
time for repeated compatible batches while bounding object/wire/RSS growth?
Unit = one paired list/compiled execution of one deterministic workload in
fresh R processes. Four storage/overlap workloads are co-primary; all must pass.
Synthetic fixtures test computation only and support no biological claim.

## Fixed fixtures

- 50,000 exact ordered features; 2,000 sets; canonical set size 25 = 50,000
  memberships before overlap. Low overlap uses the existing disjoint
  constructor; high overlap shares 75% of each set. Every 20th set replaces its
  final member with a unique absent identifier; every 25th appends a duplicate
  of its first member. Stable deduplication therefore retains exactly 25
  canonical members per set for the memory denominator.
- Dense = 24 bulk-shaped samples, 35% exact zeros, and 1% sample-specific
  `NA`/`NaN`; sparse = 240 samples at 3% stored density with 5% of stored cells
  missing. Existing deterministic generators supply non-negative finite values.
  Inputs and five batches reuse identical
  feature order and catalogue; batch `k` uses matrix seed `matrix_seed + k - 1`
  while dimensions and names stay compatible.
- Primary execution uses explicit `SerialParam()` to isolate catalogue work.
  Separate correctness/resource cases use fresh and already-started two-worker
  SOCK backends. Exactly one scheduler is active.

## Endpoints and repeats

Twenty paired repeats per workload use the pre-generated order strings in the
TSV (10 list-first, 10 compiled-first). Workloads rotate within repeat. Each
method observation runs in a fresh `R --vanilla` process after fixture creation,
package load, and `gc()`.

Primary elapsed endpoint:

- list = five complete `genefunnel(mat_batch, gene_sets, SerialParam())` calls;
- compiled = one compilation + five complete internal compiled calls;
- estimator = geometric mean compiled/list cumulative ratio at call five;
- uncertainty = one-sided 95% Student interval over paired log ratios;
- multiplicity = intersection-union decision: every workload's upper bound must
  be `<= 0.85`, so no alpha split is required to claim all-workload success.

Secondary descriptive endpoints = compile, validation, canonical encoding,
serialization/deserialization, per-call scoring, and cumulative ratios at calls
1-4; median + paired ratios are reported without rescuing the primary decision.
Amortization additionally requires median cumulative compiled time `<=` list by
call five in every workload.

Every paired call must have identical serialized score digest, dimensions,
dimnames, missingness, warnings, and error class/message for fixed adversarial
cases. Serialization bytes/fingerprints must be deterministic across fresh
processes. Fresh/reused SOCK results must equal serial results and workers must
load the assigned installed candidate.

## Resource gates

Resource passes are separate from timing and therefore cannot contend with it.
`Rprofmem()` records manager allocation; compatible GNU time plus passive
`/proc` sampling records whole-process and compilation peak RSS. Unavailable
metrics reject the performance claim, never become zero.

- `object.size(compiled) / 50,000 <= 256` bytes/membership.
- serialized `GFCAT-1` bytes / 50,000 `<= 192` bytes/membership.
- compilation peak-RSS increment above the generated list/features fixture
  `<= 16 MiB + 512 * 50,000` bytes.
- candidate scoring allocation/RSS may not exceed the list path by more than
  `max(64 KiB, 1%)` and `max(4 MiB, 10%)`, respectively, after excluding the
  intentionally retained compiled object.

Malformed object/raw mutations, stale feature value/order, duplicate feature
IDs, truncation, trailing data, digest corruption, and incompatible schema/
formula versions must fail before native scoring.

## Environment and rejection

Linux timing uses the A2 quiescence rule: one-minute load/logical CPU `<= 0.25`
at start and before each pair. Record Git SHAs, dirty state, source/installed
fingerprints, R/dependency versions and paths, OS/CPU/RAM, external time tool,
load checks, seeds, orders, raw observations, paired table, and decision table.
Any load excursion, incomplete repeat, identity failure, unavailable required
metric, or threshold failure rejects promotion; repeats are fixed and are not
extended after inspection.

## Decision

**Go:** all correctness/portability/resource gates plus all four `0.85` upper
timing bounds pass. Then design the public compiler/scorer/provenance API without
changing `genefunnel()`.

**Fallback:** keep matching preparation internal. If compilation lacks the
minimum cumulative effect, omit compiled scoring; retain only a provenance
helper if its independent contract later passes.
