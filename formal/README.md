# Formal GeneFunnel scorer

The exact Lean rational scorer satisfies the reviewed GeneFunnel scoring
contract. `GeneFunnel.Proof.implementation_correct` has the exact type
`GeneFunnel.Spec.Contract GeneFunnel.Impl.run`.
`GeneFunnel.Bounds.bounds_correct` has the exact type
`GeneFunnel.Spec.BoundsContract GeneFunnel.Impl.run`.

## Verified boundary

- `Spec.lean` defines the reviewed semantics. Optional entries model
  sample-local missingness. Observed zeros remain, and negative rationals are rejected.
  Fewer than two observations yield no score. Scoreable cells use the normative
  absolute-deviation equation.
- `Impl.lean` contains the executable evaluator. It uses the stable below-mean
  identity from the native scorer.
- `Proof.lean` proves equality for every finite rational cell and covers all
  valid and rejected behavior. `Bounds.lean` proves non-negative and
  at-most-total scores for the same executable root.
- `Audit.lean` checks exact symbol origins, types, the reviewed contract, and
  the safe executable closure. It permits only `Classical.choice`, `Quot.sound`,
  and `propext` as standard logical dependencies.
- `verify.sh` checks the source policy, review budget, and SHA-256 drift. It
  runs an isolated cache-free build and the audit. It also performs two fresh
  `leanchecker` replays with Lean 4.32.2.

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

The theorems cover only the pure exact cell model. They do not establish
binary64 rounding or error behavior, global R validation, identifier matching,
or deduplication. They also exclude storage traversal, chunking, parallel
scheduling, Rcpp/FFI, compiler and code-generation correctness, resource
bounds, biological validity, and specification correspondence. Package tests
provide separate conformance evidence for the implemented paths.
