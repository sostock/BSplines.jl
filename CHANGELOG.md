# BSplines.jl changelog

## master

* ![BREAKING](https://img.shields.io/badge/-BREAKING-red) The `BSplineBasis` type now stores the knot vector instead of the breakpoint vector. Consequentially, the type parameter `T` in `BSplineBasis{T}` now refers to the type of the knot vector.
* ![BREAKING](https://img.shields.io/badge/-BREAKING-red) The `breakpoints` function now always returns a vector of unique values and is generally not identical (`===`) to the vector that was used to create the basis.
* ![BREAKING](https://img.shields.io/badge/-BREAKING-red) The `bsplines` and `bsplines!` functions no longer return `nothing` when there are no non-zero B-splines at the specified point. Instead, both functions always return an `OffsetArray` wrapping a `view` (a `view` of the destination array in the case of `bsplines!`). When there are no non-zero B-splines at the specified point, the `view` is empty.
* ![Feature](https://img.shields.io/badge/-feature-green) ![Deprecation](https://img.shields.io/badge/-deprecation-orange) A new `BSplineBasis` constructor enables specifying either a knot vector or a breakpoint vector as keyword argument, i.e., `BSplineBasis(k, knots=vec)` specifies `vec` as the knot vector and `BSplineBasis(k, breakpoints=vec)` specifies `vec` as the breakpoint vector whose first and last elements are repeated to create the knot vector. When using the `breakpoints` keyword, `vec` can still contain repeated interior points which appear with the same multiplicity in the knot sequence. The old `BSplineBasis(k, vec)` constructor has been deprecated in favor of `BSplineBasis(k, breakpoints=vec)`.
* ![Feature](https://img.shields.io/badge/-feature-green) A `BSplineBasis` can now be sliced with an `AbstractUnitRange` to create a new basis, i.e., `basis[indices]` will create a `BSplineBasis` that contains the B-splines from `basis` with the specified indices. `view(basis, indices)` can be used as well, in which case the knot vector of the returned basis is a `view` of the original knot vector.

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
