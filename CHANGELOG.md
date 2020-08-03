# BSplines.jl changelog

## master

* ![BREAKING](https://img.shields.io/badge/-BREAKING-red) The `BSplineBasis` type now stores the knot vector instead of the breakpoint vector. Consequentially, the type parameter `T` in `BSplineBasis{T}` now refers to the type of the knot vector.
* ![BREAKING](https://img.shields.io/badge/-BREAKING-red) The `breakpoints` function now always returns a vector of unique values and is generally not identical (`===`) to the vector that was used to create the basis.
* ![Feature](https://img.shields.io/badge/-feature-green) ![Deprecation](https://img.shields.io/badge/-deprecation-orange) A new `BSplineBasis` constructor enables specifying either a knot vector or a breakpoint vector as keyword argument, i.e., `BSplineBasis(k, knots=vec)` specifies `vec` as the knot vector and `BSplineBasis(k, breakpoints=vec)` specifies `vec` as the breakpoint vector whose first and last elements are repeated to create the knot vector. When using the `breakpoints` keyword, `vec` can still contain repeated interior points which appear with the same multiplicity in the knot sequence. The old `BSplineBasis(k, vec)` constructor has been deprecated in favor of `BSplineBasis(k, breakpoints=vec)`.
* ![Feature](https://img.shields.io/badge/-feature-green) A `BSplineBasis` can now be sliced with an `AbstractUnitRange` to create a new basis, i.e., `basis[indices]` will create a `BSplineBasis` that contains the B-splines from `basis` with the specified indices. `view(basis, indices)` can be used as well, in which case the knot vector of the returned basis is a `view` of the original knot vector.
* ![Enhancement](https://img.shields.io/badge/-enhancement-blue) `length(::BSplineBasis)` now always returns a value of type `Int`.

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
