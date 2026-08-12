# Formal GeneFunnel scorer

The exact Lean rational scorer satisfies the reviewed GeneFunnel scoring
contract. `GeneFunnel.Proof.implementation_correct` has the exact type
`GeneFunnel.Spec.Contract GeneFunnel.Impl.run`.
`GeneFunnel.Bounds.bounds_correct` has the exact type
`GeneFunnel.Spec.BoundsContract GeneFunnel.Impl.run`.

## Verified boundary

- `Spec.lean` = reviewed semantics: optional entries model sample-local
  missingness; observed zeroes remain; negative rationals are rejected; fewer
  than two observations yield no score; scoreable cells use the normative
  absolute-deviation equation.
- `Impl.lean` = executable evaluator using the stable below-mean identity from
  the native scorer.
- `Proof.lean` = equality proof for every finite rational cell plus complete
  valid/rejected behavior; `Bounds.lean` = nonnegative and at-most-total score
  proof for the same executable root.
- `Audit.lean` = exact symbol-origin/type checks, reviewed-contract and safe
  executable closure, plus an allowance limited to the standard logical
  dependencies `Classical.choice`, `Quot.sound`, and `propext`.
- `verify.sh` = source-policy and review-budget checks, SHA-256 drift detection,
  isolated cache-free build, audit execution, and two fresh `leanchecker`
  replays under pinned Lean 4.32.2.

`inst/SCIENTIFIC_SPEC.md` and `src/calculateScores.cpp` are hashed review
evidence. Their drift invalidates the committed digest record. Hash binding
does not prove that either file is a faithful translation of the Lean model.

## Workflow

Prerequisites: Bash 4+, GNU coreutils, ripgrep, `elan`, and the pinned Lean
toolchain.

```bash
./formal/verify.sh
```

After an intentional verification-input change:

```bash
./formal/verify.sh --update-hashes
git diff -- formal/SHA256SUMS
./formal/verify.sh
```

`formal/verification-inputs.txt` defines the closed snapshot.
`formal/SHA256SUMS` is an unsigned local drift record; CI checks it and reruns
the proof from an isolated copy of those inputs.

## Non-claims

The theorems cover the pure exact cell model only. They do not establish
binary64 rounding/error behavior, global R validation, identifier matching or
deduplication, dense/sparse traversal, chunking, parallel scheduling, Rcpp/FFI,
compiler/code-generation correctness, resource bounds, biological validity,
or specification correspondence. Package tests provide separate conformance
evidence for the implemented paths.
