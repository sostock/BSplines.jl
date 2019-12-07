# BSplines.jl

This package provides data types and functions for working with [B-splines](https://en.wikipedia.org/wiki/B-spline) as a means to approximate real functions.

## Installation

This package is compatible with Julia ≥ 0.7. It is currently not registered and can be installed by typing
```
] add https://github.com/sostock/BSplines.jl
```
in the Julia REPL.

## Acknowledgments

The algorithms used for evaluating B-splines (and their derivatives) are adapted from 
the Fortran code found in Carl de Boor’s book *A practical Guide to Splines* [^deBoor1978],
in particular from the subroutines `BSPLVB`, `BSPLVD` and `BVALUE`.

[^deBoor1978]:
    Carl de Boor, *A Practical Guide to Splines*, New York, N.Y.: Springer-Verlag, 1978.

## Manual outline

```@contents
Pages = ["basis.md", "spline.md", "functions.md", "plotting.md"]
```
