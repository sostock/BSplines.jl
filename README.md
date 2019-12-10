# BSplines

[![Build Status](https://travis-ci.com/sostock/BSplines.jl.svg?branch=master)](https://travis-ci.com/sostock/BSplines.jl)
[![Build status](https://ci.appveyor.com/api/projects/status/ruh7o1yalohqawbd/branch/master?svg=true)](https://ci.appveyor.com/project/sostock/bsplines-jl/branch/master)
[![codecov](https://codecov.io/gh/sostock/BSplines.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/sostock/BSplines.jl)
[![Coverage Status](https://coveralls.io/repos/github/sostock/BSplines.jl/badge.svg?branch=master)](https://coveralls.io/github/sostock/BSplines.jl?branch=master)

[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://sostock.github.io/BSplines.jl/dev)

This package provides data types and functions for working with [B-splines](https://en.wikipedia.org/wiki/B-spline) as a means to approximate real functions.
For its usage, see the [documentation](https://sostock.github.io/BSplines.jl/dev).

**This package is not related to the package at [https://github.com/gusl/BSplines.jl](https://github.com/gusl/BSplines.jl).
While [gusl/BSplines.jl](https://github.com/gusl/BSplines.jl) is only compatible with Julia ≤ 0.4, this package works with Julia ≥ 1.0.**

## Installation

This package is compatible with Julia ≥ 1.0. It is currently not registered and can be installed by typing
```
] add https://github.com/sostock/BSplines.jl
```
in the Julia REPL.

## Acknowledgments

The algorithms used for evaluating B-splines (and their derivatives) are adapted from 
the Fortran code found in Carl de Boor’s book *A practical Guide to Splines* [[1]](#deBoor1978),
in particular from the subroutines `BSPLVB`, `BSPLVD` and `BVALUE`.

<a name="deBoor1978">[1]</a> Carl de Boor, *A Practical Guide to Splines*, New York, N.Y.: Springer-Verlag, 1978.
