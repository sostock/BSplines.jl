# The `BSplineBasis` type and related functions

```@setup basis
using BSplines
```

To work with B-spline bases, this package defines the parametric [`BSplineBasis{T}`](@ref) type which represents B-spline bases with breakpoint vector of type `T`.

A B-spline basis is completely characterized by its order ``k`` and knot vector ``t``.
In the case of the `BSplineBasis` type, **the knot vector of a basis is generated from its breakpoint vector by repeating the first and last breakpoints so that they appear ``k`` times.**

!!! note
    Knot sequences where the first and last knot do not appear ``k`` times are not supported by the `BSplineBasis` type.
    The reason for this is that it simplifies the evaluation of B-splines, since it means that exactly ``k`` B-splines are non-zero on each interval.
    If the first or last knot would appear less than ``k`` times, this would not be the case.

## Properties of B-spline bases -- `order`, `breakpoints`, and `knots`

A `BSplineBasis` is constructed from its order and breakpoint vector.
The [`order`](@ref) and [`breakpoints`](@ref) functions return these properties.
The `length` function returns the number of B-splines in the basis:
```@repl basis
basis = BSplineBasis(4, 0:5)
order(basis)
breakpoints(basis)
length(basis)
```
The [`knots`](@ref) function returns the knot vector that is generated from the breakpoints.
In order to not allocate memory for a new array, a wrapper type around the original breakpoint vector is used:
```@repl basis
knots(basis)
```

The interval over which a `BSplineBasis` is defined (i.e., the interval between the first and the last breakpoint) can be obtained (as a `Tuple`) with the [`support`](@ref) function:
```@repl basis
support(basis)
```


## Evaluating B-Splines and their derivatives

The [`bsplines`](@ref) and [`bsplines!`](@ref) functions can be used to obtain the values (or derivatives) of all B-splines that are non-zero at a given point.

### Evaluating B-splines 

If `x` is within the support of the B-spline basis, `bsplines(basis, x)` returns an
[`OffsetArray`](https://github.com/JuliaArrays/OffsetArrays.jl)
which contains the value of the `i`-th B-spline at the index `i`.
The array always has the length `order(basis)`:

```@repl basis
basis = BSplineBasis(4, 0:5);
bsplines(basis, 3.2)
bsplines(basis, 7//3)
```

The type of the elements of the array depends on the type of the breakpoints and the type of `x`.
If the point `x` is outside of the support of `basis`, `nothing` is returned instead:

```@repl basis
bsplines(basis, 6) # returns `nothing`, which is not printed in the REPL
```

### Evaluating derivatives

The `bsplines` function can also calculate derivatives of the B-splines.
The `N`-th derivative is specified via an optional third argument of the singleton type [`Derivative{N}`](@ref).
Instead of `Derivative{N}()`, the shorter constructor `Derivative(N)` can be used:

```@repl basis
basis = BSplineBasis(4, 0:5);
bsplines(basis, 7//3, Derivative(1)) # calculate first derivatives of non-zero B-splines
```

To calculate the value of the non-zero B-splines as well as all of their derivatives up to the `N-1`-th, the [`AllDerivatives{N}`](@ref) singleton type can be used instead.
`AllDerivatives(N)` is equivalent to `AllDerivatives{N}()`.
If `x` is within the support of the B-spline basis, `bsplines(basis, x, AllDerivatives(N))` returns an `OffsetArray` which contains the `j`-th derivative of the `i`-th B-spline at the index `i, j` (the zeroth derivative is the value itself):

```@repl basis
bsplines(basis, 4, AllDerivatives(3)) # calculate zeroth, first and second derivatives
```

### Pre-allocating output arrays

The `bsplines` function allocates a new array to return the values (except when it returns `nothing`).
In order to write the values to a pre-allocated array instead, the [`bsplines!`](@ref) function can be used.

The `bsplines!` function returns an integer `offset` such that the `i`-th element of the destination array contains the value (or derivative) of the `i+offset`-th B-spline.
In the case of `AllDerivatives{N}`, the destination array contains the `j-1`-th derivative of the `i+offset`-th B-spline at the index `i, j`.
If the point `x` is outside of the support of the basis, `nothing` is returned instead and the destination array is not mutated.

```@repl basis
basis = BSplineBasis(4, 0:5);
vec = zeros(4);
bsplines!(vec, basis, 3.2)
vec
bsplines!(vec, basis, 7//3, Derivative(2))
vec
mat = zeros(Rational{Int}, 4, 2);
bsplines!(mat, basis, 4, AllDerivatives(2))
mat
```

When only calculating B-splines or their derivatives via `Derivative{N}`, the destination array must be of length `order(basis)`.
In the case of the `AllDerivatives{N}` argument, the destination array must be an `order(basis) × N` matrix.
In any case, the destination array must not have offset axes.

### Specifying the relevant interval

Evaluating B-splines at a point `x` requires finding the largest index `i` such that `t[i] ≤ x` and `t[i] < t[end]` where `t` is the knot vector.
The `bsplines` and `bsplines!` functions use the [`intervalindex`](@ref) function to find this index.
If the index is already known, it can be specified with the `leftknot` keyword argument to `bsplines`/`bsplines!` in order to speed up the computation.

## Constructing knot vectors

As stated above, the `knots` of a `BSplineBasis` are generated from its `breakpoints` by repeating the first and last breakpoints so that they appear `order(basis)` times.
However, in many situations it may be desirable to
  * repeat other knots than the first and last or
  * have the first and/or last knot repeated less than ``k`` times (e.g. to implement certain boundary conditions).

In order to have repeated interior knots, it is sufficient to include it in the breakpoint vector multiple times when constructing the basis:

```@repl basis
basis = BSplineBasis(3, [1, 2, 3, 4, 5, 5, 6, 7, 8]) # 5 appears twice
knots(basis)
```

There is currently no way to create a `BSplineBasis` where the first and last knots appear less than ``k`` times.
However, some functions provided by this package (see [Higher-level functions](@ref)) provide a `indices` keyword argument to select only a certain range of B-splines from a basis, thus achieving the same result as if the knot vector had been shortened.

Because of repeated knots, not every pair of knots `(t[i], t[i+1])` describes an actual interval.
The [`intervalindices`](@ref) helps with finding all indices which describe relevant intervals.
See its docstrings for more information.
