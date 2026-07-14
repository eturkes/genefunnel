# Project memory

- Development R packages live in `.agent/R-library`; prefix package commands
  with `R_LIBS_USER="$PWD/.agent/R-library"`.
- Install into that library serially. Concurrent `install.packages()` processes
  contend on `00LOCK-*` and can leave an incomplete toolchain.
