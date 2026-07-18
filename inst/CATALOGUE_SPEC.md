<!-- Assisted-by: OpenAI Codex. -->

# Compiled catalogue spike contract

**Status:** internal research contract, schema `GFCAT-1`, formula `GF-1`.
No function or class described here is a supported public API. Promotion waits
for protocol `C-1.0.2`; failure retains the named-list path and may retain only
provenance machinery.

## Purpose and authority

A compiled catalogue caches exact identifier matching and a portable
row-to-set adjacency for repeated matrices with one ordered feature universe.
It changes no score. `genefunnel(mat, gene_sets, BPPARAM)` remains authoritative
and never accepts the compiled object. [`SCIENTIFIC_SPEC.md`](SCIENTIFIC_SPEC.md)
owns scoring; this file owns only the spike's state, integrity, reuse, and wire
contracts.

## Canonical state

Construction validates through the named-list contract, then copies value-only
state. Member duplicates are removed stably; element `names()` on feature/member
vectors are non-semantic and discarded. The object contains:

1. schema/formula versions;
2. exact ordered feature identifiers;
3. ordered set names + canonical member identifiers;
4. duplicate-member counts and the ordinary coverage table;
5. per-member feature indices (`0` = globally unmatched), matched per-set
   indices, and row-to-set compressed adjacency;
6. SHA-256 fingerprints for the catalogue, feature universe, and full compiled
   state.

Adjacency is portable R integer state: zero-based `offsets` of length
`length(features) + 1`, then one-based set indices ordered by feature row and
canonical set order. Native pointers and environments are absent. Construction
exposes no mutator; R copy-on-modify plus validation on every use provides
logical, not physical, immutability.

Every score use validates the complete schema, field types/order, identifier
contract, index ranges, matched identities, derived indices/coverage/adjacency,
and all fingerprints. It then compares matrix row names with the stored feature
vector in exact order. Any mismatch is stale reuse and fails before scoring.
Fingerprints accelerate integrity checks and identify content; exact validation
carries safety. They are neither origin proof nor a security signature. A caller
with package internals can deliberately construct different self-consistent
content, so provenance must not claim otherwise.

## Canonical bytes and fingerprints

R 4.5 and later - including the current
[Bioconductor 3.23 / R 4.6 runtime](https://www.bioconductor.org/news/bioc_3_23_release/) - provide
[`tools::sha256sum(bytes=)`](https://stat.ethz.ch/R-manual/R-patched/library/tools/html/sha256sum.html)
using a public-domain SHA-256 implementation. It is preferred over a new
third-party hash dependency, an OpenSSL system surface, or duplicate package
cryptography. SHA-256 collision resistance supports durable content labels;
exact state checks remain mandatory.

The package does not fingerprint `serialize()` output. R's
[`serialize()` contract](https://stat.ethz.ch/R-manual/R-devel/library/base/help/serialize.html)
uses portable XDR by default but explicitly reserves future format changes and
discourages long-term object storage. `GFCAT-1` instead frames these bytes:

```text
magic "GFCAT" | schema:u32be | formula:id | features:id[] |
sets:{name:id, duplicate_count:u32be, members:id[], member_index:u32be[]}[] |
adjacency_offsets:u32be[] | adjacency_set_indices:u32be[] | sha256(body)
```

All counts/indices are unsigned big-endian 32-bit values constrained to R's
non-negative integer range before allocation. Each identifier is
`encoding-tag:u8 | byte-length:u32be | bytes`; tags distinguish R's `unknown`,
`UTF-8`, `latin1`, and `bytes` representations. Values, names, set/member order,
and encoding marks are preserved; Unicode normalization and identifier mapping
are absent. Invalid tags, impossible lengths, truncation, trailing bytes,
fingerprint mismatch, or a derived-state mismatch fail closed.

The schema contains no floating-point field: empty-set coverage missingness is
derived from exact integer counts, and signed zero cannot be conflated. Any
future floating field must define missing/NaN bit canonicalization and signed
zero explicitly under a new schema. Result/provenance digests must frame set,
sample, and streamed viewport order explicitly; catalogue fingerprints alone do
not identify a result.

Domain-separated encodings produce:

- `catalogue`: ordered names, canonical members, duplicate counts, formula;
- `features`: ordered identifier universe; and
- `content`: the complete `GFCAT-1` body, including mappings and adjacency.

Custom deserialization verifies the appended content digest before parsing,
applies strict remaining-byte/allocation bounds, reconstructs matching and
adjacency from primary identifiers, and accepts the payload only when both
representations agree. Supported-platform tests fix byte/digest vectors and
exercise fresh R/SOCK processes.

## Scoring spike

The internal compiled entry point accepts `(mat, catalogue, BPPARAM)` only. Its
error order is matrix structure/value validation, catalogue validity + exact
feature reuse, then backend validation. It applies the same retained-set
universe, aggregate omission warning, chunks, native kernels, assembly, names,
and plain numeric-matrix return as the named-list path. Workers receive portable
integer memberships; fresh/reused SOCK workers never receive a native pointer.

No automatic catalogue compilation, implicit coercion, provenance attachment,
normalization, mapping, or cache exists. Caller labels and a public provenance
sidecar remain later work.

## Promotion gate

The repository's
[`catalogue-protocol.md`](https://github.com/eturkes/genefunnel/blob/main/benchmark/catalogue-protocol.md)
freezes the spike workloads and decision. Required outcomes: exact list-path
results/errors; deterministic portable bytes; malformed/stale rejection;
fresh/reused SOCK identity; compilation amortized by call five; a one-sided 95%
upper cumulative-time ratio at most `0.85` in every workload; compiled-object
size at most 256 bytes per canonical membership; serialized size at most 192
bytes per membership; and bounded compilation RSS. One failure keeps compiled
scoring internal and triggers the stated fallback.
