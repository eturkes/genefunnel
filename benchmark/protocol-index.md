<!-- Assisted-by: OpenAI Codex. -->

# Executable benchmark protocol index F-I-1.0.0

The machine [`index`](protocol-index.tsv), SHA-256
`3d34bf1bf40cadc0f64fdbc65fbe0e7ac6fff3e888a6cb74a6c069fc2d67577b`,
maps an explicit protocol version and suite to its complete tracked execution
closure. The harness accepts no implicit latest version. It validates the index,
every listed file, the embedded manifest version, and the unique runner before
starting a fresh R process.

Protocol `1.0.0` freezes ten LF-terminated files: its machine manifest/parser,
fixture/provenance/report libraries, performance runner/worker/measurement
library, and controlled runner/assertion library. Repository attributes force
LF checkout for benchmark R/TSV files, so the hashes do not depend on Windows
line-ending policy. Base R `tools::sha256sum()` supplies the byte fingerprints;
base `system2()` plus `shQuote()` supplies portable synchronous dispatch. This
adds no package or benchmark dependency.

The hashes detect content drift and identify registered bytes. They are not a
signature, origin proof, trust mechanism, or guarantee that the scientific
design is adequate. The Git commit and generated environment evidence remain
separate provenance facts.

Run an exact registered suite with:

```sh
Rscript --vanilla benchmark/run-protocol.R \
  --protocol=1.0.0 --suite=controlled

Rscript --vanilla benchmark/run-protocol.R \
  --protocol=1.0.0 --suite=performance --preset=smoke \
  --repeats=1 --workers=2
```

`--list` reports registered pairs. All arguments other than harness options pass
unchanged to the selected runner. A non-zero runner status is returned unchanged.

Protocol `1.0.0` contains analytic score assertions and descriptive synthetic
resource workloads only. It contains no competitor, external assay, biological
endpoint, inferential test, or thesis-data reproduction. An implementation-only
repair requires a new registered patch version with preserved scientific
meaning; any new method, scenario, or assertion starts protocol `2.0.0` and is
committed before its comparative result is inspected. Old version files remain
immutable and executable.

Prospective parent `F-2.0.0` and its immutable effective-design amendment
`F-2.1.0` are tracked separately in
[`validation-protocol.tsv`](validation-protocol.tsv) and
[`validation-amendment.tsv`](validation-amendment.tsv). They name future
methods, assays, assertions, fixed-panel scope, and selection/acquisition
contracts but have no comparative runner or method implementation, no
selected-source bundle or molecular data, and no results. F-2.1.0 is a 128-row
delta/266-row effective contract with a fixed 5,062-row Reactome 97 roster. Its
hermetic structural fixture currently covers 12 sources, 90 tasks, 60 null
contrasts, 29 objects, 37 F-A events, and 99 adversaries. Contracts include
explicit `field_type`, byte-exact UTF-8 scalar and underflow-aware numeric
parsing, dependence bijection, total archive decoder/resource and finite
acquisition bounds, guarded object descriptors, and exact S/B/I tree/runtime/
mode Git evidence. They assume a trusted fresh process and do not sandbox
project code. No candidate source, archive, or molecular table was selected or
accessed by this work; that process record is not proof against external
access. The committed readiness sequence is hashed F-S selection → F-B
opaque-byte closure → F-E implementation closure → controlled runner before
each F-A attempt. F-A evidence spans acquisition through held-out first access;
production validation requires `verify_git = TRUE` and the exact Git/object,
F-E, and runner closures. F-2.1.0 remains non-executable until those exact
artifacts exist. Dispatch continues to expose only `1.0.0`.
