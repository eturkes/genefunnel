# GeneFunnel

**GeneFunnel** is an R gene-set scoring algorithm for feature-by-sample
expression matrices. Dense inputs stay dense; sparse inputs stay sparse through
bounded scoring chunks. Reproducible dense, sparse, serial, and parallel
performance measurements are defined in the [benchmark harness](benchmark/README.md).

```
This file is part of GeneFunnel.
Copyright (C) 2025  Emir Turkes, UK DRI at UCL

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

Emir Turkes can be contacted at emir.turkes@eturkes.com
```

## Installation

GeneFunnel is currently in development. You can install the development version from GitHub using:

```
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("eturkes/genefunnel")
```
