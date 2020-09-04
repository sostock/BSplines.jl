# BSplines.jl changelog

## master

* ![Maintenance](https://img.shields.io/badge/-maintenance-grey) Fix tests on Julia ≥ 1.6. ([#18](https://github.com/sostock/BSplines.jl/pull/18))

## v0.3.0

* ![BREAKING](https://img.shields.io/badge/-BREAKING-red) `bsplines!(dest, args...)` now returns an `OffsetArray` that wraps `dest`, making its output equal to that of `bsplines(args...)`. ([#8](https://github.com/sostock/BSplines.jl/pull/8))
* ![Feature](https://img.shields.io/badge/-feature-green) `splinevalue` now accepts a keyword argument `workspace` for providing a vector to store intermediate values in order to avoid unnecessary allocations. ([#10](https://github.com/sostock/BSplines.jl/pull/10))
* ![Feature](https://img.shields.io/badge/-feature-green) When calculating derivatives, `bsplines` and `bsplines!` now accept a keyword argument `derivspace` for providing a matrix to store intermediate values in order to avoid unnecessary allocations. ([#16](https://github.com/sostock/BSplines.jl/pull/16))
* ![Feature](https://img.shields.io/badge/-feature-green) `BSplineBasis` and `IntervalIndices` now support reverse iteration via `Iterators.reverse`. ([#15](https://github.com/sostock/BSplines.jl/pull/15))
* ![Enhancement](https://img.shields.io/badge/-enhancement-blue) The default printing of `BSplineBasis` and `Spline` now uses the compact style for printing the breakpoint and coefficient vectors. ([#9](https://github.com/sostock/BSplines.jl/pull/9))
* ![Bugfix](https://img.shields.io/badge/-bugfix-purple) `BSplineBasis` and `IntervalIndices` now implement `IteratorSize` and `eltype` correctly. ([#14](https://github.com/sostock/HalfIntegers.jl/pull/14))

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
