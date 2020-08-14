# BSplines.jl changelog

## master

* ![BREAKING](https://img.shields.io/badge/-BREAKING-red) `bsplines!(dest, args...)` now returns an `OffsetArray` that wraps `dest`, making its output equal to that of `bsplines(args...)`. ([#8](https://github.com/sostock/BSplines.jl/pull/8))

* ![Feature](https://img.shields.io/badge/-feature-green) `splinevalue` now accepts a keyword argument `workspace` for providing a vector to store intermediate values in order to avoid unnecessary allocations. ([#10](https://github.com/sostock/BSplines.jl/pull/10))

* ![Enhancement](https://img.shields.io/badge/-enhancement-blue) The default printing of `BSplineBasis` and `Spline` now uses the compact style for printing the breakpoint and coefficient vectors. ([#9](https://github.com/sostock/BSplines.jl/pull/9))


## v0.2.5

* ![Enhancement](https://img.shields.io/badge/-enhancement-blue) `length(::BSplineBasis)` now always returns a value of type `Int`. ([#7](https://github.com/sostock/BSplines.jl/pull/7))

## v0.2.4

* ![Maintenance](https://img.shields.io/badge/-maintenance-grey) Fixed the tests on Julia ≥ 1.6.

## v0.2.3

* ![Maintenance](https://img.shields.io/badge/-maintenance-grey) Compatibility with RecipesBase v1.

## v0.2.2

* ![Maintenance](https://img.shields.io/badge/-maintenance-grey) Compatibility with RecipesBase v0.8.

## v0.2.1

* ![Maintenance](https://img.shields.io/badge/-maintenance-grey) Compatibility with OffsetArrays v1.

## v0.2.0

* ![Feature](https://img.shields.io/badge/-feature-green) New functions `basismatrix` and `basismatrix!`.

## v0.1.0

Initial release.
