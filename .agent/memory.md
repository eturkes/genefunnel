# Project memory

- Formal scoring proof = `formal/GeneFunnel/*` + `formal/verify.sh`, Lean
  4.32.2. The theorems cover exact rational cell equivalence + bounds;
  native/R/storage/parallel paths remain conformance-tested boundaries. After
  a verification-input change: `./formal/verify.sh --update-hashes`, inspect
  the digest diff, then rerun `./formal/verify.sh`.
- Development R packages live in `.agent/R-library`; prefix package commands
  with `R_LIBS_USER="$PWD/.agent/R-library"`.
- Install into that library serially. Concurrent `install.packages()` processes
  contend on `00LOCK-*` and can leave an incomplete toolchain.
- Validate `Matrix` objects through represented columns, not raw `@x`: structured
  dense classes can hold ignored values and triplet sparse duplicates aggregate.
- ASan-instrumented package install/load requires
  `LD_PRELOAD="$(g++ -print-file-name=libasan.so)"`; disable leak detection to
  avoid R-runtime lifetime noise. UBSan needs no preload.
- Workstream F's effective prospective design is immutable parent `F-2.0.0`
  plus amendment `F-2.1.0`. No candidate source is selected and no candidate
  archive or molecular-table access is recorded or known; this process record
  is not proof of non-access.
- F-2.1.0 pins the generated 5,062-row Reactome 97 eligible-target roster and
  128-row amendment/266-row effective contract plus machine-validated assay-
  input/access-role schemas. Its current hermetic fixture = 12 sources, 90
  tasks, 60 nulls, 29 objects, 37 F-A events, 99 adversaries. It closes
  `field_type`, byte-exact UTF-8 scalar/underflow-aware numeric parsing, total
  decoder/resource and finite acquisition bounds, dependence bijection, guarded
  descriptors, and exact S/B/I tree/runtime/mode evidence. Structural guarantees
  assume a trusted fresh process; they are not a sandbox. Future execution
  order is `F-S-1.0.0` selection → `F-B-1.0.0` opaque bytes → `F-E-1.0.0`
  implementation closure → controlled runner before each `F-A-1.0.0` attempt;
  F-A covers acquisition through held-out first access. Content-sealed is a
  declared guarded process state, not cryptographic or OS isolation. Production
  F-A requires `verify_git = TRUE`, exact Git/object closure, F-E, and the runner;
  the design remains non-executable until those exact artifacts exist.
