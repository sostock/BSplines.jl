# Higher-level functions

The functions listed below help with common use cases of B-splines.
See their docstrings for more information.

### `knotaverages`/`knotaverages!`

The [`knotaverages`](@ref) function returns a vector ``\tau`` of knot averages
```math
\tau_i = \frac{1}{k-1} \sum_{j=i+1}^{i+k-1} t_j
```
where ``t`` is the knot vector of the B-spline basis and ``k`` is its order.
The length of ``\tau`` equals the number of B-splines in the basis.
The knot averages are recommended in [^deBoor1978] (p. 214) as data points for interpolation.
Instead of creating the knot vector for all indices ``i``, a range of indices can be supplied with the keyword argument `indices`.

The [`knotaverages!`](@ref) function can be used to write the knot averages to a pre-allocated array.

[^deBoor1978]: Carl de Boor, *A Practical Guide to Splines*, New York, N.Y.: Springer, 1978.

### `averagebasis`

The [`averagebasis`](@ref) function returns a B-spline basis of a specified order that is well-suited for interpolating a function at a given set of data points.

### `interpolate`

The [`interpolate`](@ref) function interpolates data (vectors of ``x`` and ``y`` values) in a given B-spline basis.
It returns a `Spline`.
The `indices` keyword can be used to restrict the interpolation to a range of B-splines from the basis.

### `approximate`

The [`approximate`](@ref) function approximates a function ``f:\mathbb{R}\to\mathbb{R}`` in a given B-spline basis by sampling the function at the knot averages of the basis and interpolating the samples.
The `indices` keyword can be used to restrict the B-splines used for the interpolation.
