# The `Spline` type and related functions

```@setup spline
using BSplines
```

A spline consists of a B-spline basis and a coefficient vector of the same length.
The [`Spline{B, C}`](@ref) type represents a spline with a basis of type `B` and a coefficient vector of type `C`.
The constructor `Spline(basis, coeffs)` returns a spline with the given B-spline basis and coefficients.
The basis on which the spline is defined can be obtained with the [`basis`](@ref) function.
The [`coeffs`](@ref) function returns the coefficient vector of a spline.

```@repl spline
b = BSplineBasis(4, 0:5);
c = rand(length(b));
spl = Spline(b, c)
basis(spl) === b
coeffs(spl) === c
```

## Evaluating `Spline`s and their derivatives

To evaluate a `Spline`, it can be called as a function.
Alternatively, the [`splinevalue`](@ref) function can be used:

```@repl spline
spl = Spline(BSplineBasis(4, 0:5), rand(8));
spl(2.5)
splinevalue(spl, 2.5)
```

To evaluate derivatives of a spline, the optional [`Derivative(N)`](@ref) argument is used:

```@repl spline
spl(2.5, Derivative(1))
splinevalue(spl, 2.5, Derivative(2))
```

!!! info
    The [`AllDerivatives{N}`](@ref) type is not supported as an argument.

Evaluating a spline at a point `x` requires finding the largest index `i` such that `t[i] ≤ x` and `t[i] < t[end]` where `t` is the knot vector of the basis.
If the index is already known, it can be specified with the `leftknot` keyword argument to speed up the computation.

## Arithmetic with `Spline`s

A `Spline` can be multiplied or divided by a real number using `*`, `/`, or `\`.
These operations return a new spline:

```@repl spline
spl = Spline(BSplineBasis(3, 0:5), [1:7;])
3.0 * spl
spl / 4
4 \ spl == ans
```

Using the `lmul!`/`rmul!` and `ldiv!`/`rdiv!` functions from the `LinearAlgebra` stdlib, a `Spline` can be multiplied/divided by a real number *in-place*, which modifies the spline (if its coefficient vector is mutable):

```@repl spline
using LinearAlgebra: rmul!
spl = Spline(BSplineBasis(3, 0:5), [1:7;])
rmul!(spl, 3);
spl
```

If two splines are defined on the same basis, they can be added and subtracted from each other using `+` and `-`:

```@repl spline
spl1 = Spline(BSplineBasis(3, 0:5), [1:7;]);
spl2 = Spline(BSplineBasis(3, 0:5), [2:8;]);
spl1 + spl2
spl1 - spl2
```

## The `BSpline` type

Indexing into or iterating over a B-spline basis yields `BSpline`s.
The parametric type [`BSpline{T}`](@ref) represents a single B-spline from a B-spline basis of type `T`:
```@repl spline
b = BSplineBasis(4, 0:5);
b[3]
```

!!! note
    Actually, a `BSpline` is just a `Spline` with a coefficient vector of a special type, i.e., `BSpline{T}` is an alias for `Spline{T,BSplines.StandardBasisVector{Bool}}`.

`BSpline`s can be used like any other `Spline`.
However, they behave differently in two ways:
  * `BSpline`s are printed differently: instead of the coefficient vector, its index within the basis is shown.
  * Calling [`support`](@ref) on a `BSpline` yields its actual support (the interval on which it is non-zero) instead of the support of the basis.
The second point is illustrated by the following example: even though `spl` and `bspl` compare equal, their `support` differs:
```@repl spline
spl = Spline(b, [0,0,1,0,0,0,0,0])
bspl = b[3]
spl == bspl
support(spl) # support of the basis
support(bspl) # actual support of the B-spline
```
