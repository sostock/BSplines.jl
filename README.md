# BSplines

[![PkgEval](https://juliaci.github.io/NanosoldierReports/pkgeval_badges/B/BSplines.svg)](https://juliaci.github.io/NanosoldierReports/pkgeval_badges/report.html)
[![CI](https://github.com/sostock/BSplines.jl/workflows/CI/badge.svg)](https://github.com/sostock/BSplines.jl/actions?query=workflow%3ACI)
[![codecov](https://codecov.io/gh/sostock/BSplines.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/sostock/BSplines.jl)

[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://sostock.github.io/BSplines.jl/stable)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://sostock.github.io/BSplines.jl/dev)

This package provides data types and functions for working with [B-splines](https://en.wikipedia.org/wiki/B-spline) as a means to approximate real functions.
For its usage, see the [documentation](https://sostock.github.io/BSplines.jl/stable).

**This package is not related to the package at [https://github.com/gusl/BSplines.jl](https://github.com/gusl/BSplines.jl).
While [gusl/BSplines.jl](https://github.com/gusl/BSplines.jl) is only compatible with Julia ≤ 0.4, this package works with Julia ≥ 1.0.**

## Installation

This package is compatible with Julia ≥ 1.0. It can be installed by typing
```
] add BSplines
```
in the Julia REPL.

## Acknowledgments

The algorithms used for evaluating B-splines and their derivatives are adapted from
the Fortran code found in Carl de Boor’s book *A practical Guide to Splines* [[1]](#deBoor1978),
in particular from the subroutines `BSPLVB`, `BSPLVD` and `BVALUE`.
This package is published with his permission.

<a name="deBoor1978">[1]</a> Carl de Boor, *A Practical Guide to Splines*, New York, N.Y.: Springer-Verlag, 1978.
