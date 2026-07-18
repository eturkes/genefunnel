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
