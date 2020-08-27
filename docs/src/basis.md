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
basis = BSplineBasis(4, breakpoints=0:5)
order(basis)
breakpoints(basis)
length(basis)
```

!!! warning
    To obtain a valid B-spline basis, the breakpoint vector must be sorted in ascending order.
    This is not checked by the `BSplineBasis` constructor.

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
basis = BSplineBasis(4, breakpoints=0:5);
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
basis = BSplineBasis(4, breakpoints=0:5);
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
In order to use a pre-allocated array instead, the [`bsplines!`](@ref) function can be used: `bsplines!(dest, args...)` behaves like `bsplines(args...)`, but the calculations are done in `dest` and the returned `OffsetArray` is a wrapper around `dest`.

```@repl basis
basis = BSplineBasis(4, breakpoints=0:5);
destvec = zeros(4);
bsplines!(destvec, basis, 3.2)
parent(ans) === destvec
bsplines!(destvec, basis, 7//3, Derivative(2))
parent(ans) === destvec
destmat = zeros(Rational{Int}, 4, 2);
bsplines!(destmat, basis, 4, AllDerivatives(2))
parent(ans) === destmat
```

When calculating values of B-splines or their derivatives via the `Derivative{N}` argument, `dest` must be a vector of length `order(basis)`.
In the case of the `AllDerivatives{N}` argument, it must be a matrix of size `(order(basis), N)`.
In any case, `dest` must not have offset axes itself.

### Improving performance

The `bsplines` and `bsplines!` functions accept two optional keyword arguments (one of them only when evaluating derivatives) can be used to speed up the evaluation of splines:

* `derivspace`:
  When evaluating derivatives, coefficients are stored in a matrix of size `(order(basis), order(basis))`.
  In order to avoid allocating a new matrix every time, a pre-allocated matrix can be supplied with the `derivspace` argument. It can only be used when calculating derivatives, i.e., with `Derivative(N)` where `N ≥ 1` or `AllDerivatives(N)` where `N ≥ 2`.
* `leftknot`:
  Evaluating B-splines at a point `x` requires finding the largest index `i` such that `t[i] ≤ x` and `t[i] < t[end]` where `t` is the knot vector.
  The `bsplines` and `bsplines!` functions use the [`intervalindex`](@ref) function to find this index.
  If the index is already known, it can be specified with the `leftknot` keyword argument.

```@repl basis
using BenchmarkTools
basis = BSplineBasis(4, breakpoints=0:5);
dest = zeros(order(basis));
space = zeros(order(basis), order(basis));
left = intervalindex(basis, 2.5);
@btime bsplines!($dest, $basis, 2.5, Derivative(1));
@btime bsplines!($dest, $basis, 2.5, Derivative(1), derivspace=$space);
@btime bsplines!($dest, $basis, 2.5, Derivative(1), derivspace=$space, leftknot=$left);
```

## Constructing knot vectors

As stated above, the `knots` of a `BSplineBasis` are generated from its `breakpoints` by repeating the first and last breakpoints so that they appear `order(basis)` times.
However, in many situations it may be desirable to
  * repeat other knots than the first and last or
  * have the first and/or last knot repeated less than ``k`` times (e.g. to implement certain boundary conditions).

In order to have repeated interior knots, it is sufficient to include it in the breakpoint vector multiple times when constructing the basis:

```@repl basis
basis = BSplineBasis(3, breakpoints=[1, 2, 3, 4, 5, 5, 6, 7, 8]) # 5 appears twice
knots(basis)
```

There is currently no way to create a `BSplineBasis` where the first and last knots appear less than ``k`` times.
However, some functions provided by this package (see [Higher-level functions](@ref)) provide a `indices` keyword argument to select only a certain range of B-splines from a basis, thus achieving the same result as if the knot vector had been shortened.

Because of repeated knots, not every pair of knots `(t[i], t[i+1])` describes an actual interval.
The [`intervalindices`](@ref) function helps with finding all indices which describe relevant intervals.
See its docstrings for more information.
