# Project memory

- Development R packages live in `.agent/R-library`; prefix package commands
  with `R_LIBS_USER="$PWD/.agent/R-library"`.
- Install into that library serially. Concurrent `install.packages()` processes
  contend on `00LOCK-*` and can leave an incomplete toolchain.
- Validate `Matrix` objects through represented columns, not raw `@x`: structured
  dense classes can hold ignored values and triplet sparse duplicates aggregate.
- ASan-instrumented package install/load requires
  `LD_PRELOAD="$(g++ -print-file-name=libasan.so)"`; disable leak detection to
  avoid R-runtime lifetime noise. UBSan needs no preload.
