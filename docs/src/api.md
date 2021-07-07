# API documentation

## Types and constructors

```@docs
AllDerivatives
AllDerivatives(::Integer)
BSpline
BSpline{B}(::B, ::Integer) where B<:BSplineBasis
BSplineBasis
BSplineBasis(::Integer; breakpoints, knots)
Derivative
Derivative(::Integer)
Spline
Spline(::Any, ::Any)
Function
```
## Functions

```@docs
approximate
averagebasis
basis
basismatrix
basismatrix!
breakpoints
bsplines!
bsplines
coeffs
interpolate
intervalindex
intervalindices
knotaverages!
knotaverages
knots
order
splinevalue
support
```

## Extended `Base` functions

```@docs
Base.view(::BSplineBasis, ::AbstractUnitRange)
```
