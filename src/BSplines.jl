module BSplines

import LinearAlgebra

using Base: @propagate_inbounds, has_offset_axes
using LinearAlgebra: I, ldiv!, lmul!, rdiv!, rmul!
using OffsetArrays: OffsetArray

export
    # Types
    AllDerivatives,
    BSpline,
    BSplineBasis,
    Derivative,
    Spline,

    # Functions
    approximate,
    averagebasis,
    basis,
    basismatrix,
    basismatrix!,
    breakpoints,
    bsplines,
    bsplines!,
    coeffs,
    interpolate,
    intervalindex,
    intervalindices,
    knotaverages,
    knotaverages!,
    knots,
    order,
    splinevalue,
    support

include("vectortypes.jl")
include("derivatives.jl")
include("bsplinebasis.jl")
include("spline.jl")
include("plotting.jl")

bspline_returntype(basis::BSplineBasis, x::Number) = bspline_returntype(basis, typeof(x))
bspline_returntype(basis::BSplineBasis, types::Type...) =
    bspline_returntype(eltype(breakpoints(basis)), types...)

bspline_returntype(spline::Spline, x::Number) = bspline_returntype(spline, typeof(x))
bspline_returntype(spline::Spline, xtype::Type) =
    bspline_returntype(basis(spline), eltype(coeffs(spline)), xtype)

function bspline_returntype(knottype::Type, types::Type...)
    T = promote_type(knottype, types...)
    typeof(one(T)/one(knottype))
end


"""
    approximate(f, basis::BSplineBasis; indices=eachindex(basis)) -> Spline

Approximate the function `f` in the B-spline basis `basis`. If `indices` is supplied, only
the basis functions at the given indices of `basis` are used.

The approximation is calculated by interpolating samples of `f` at the
[`knotaverages`](@ref) of the basis.

See also: [`interpolate`](@ref)

# Examples

```jldoctest
julia> basis = BSplineBasis(5, 0:5);

julia> spl = approximate(sin, basis, indices=2:length(basis))
Spline{BSplineBasis{UnitRange{Int64}},Array{Float64,1}}:
 basis: 9-element BSplineBasis{UnitRange{Int64}}:
  order: 5
  breakpoints: 0:5
 coeffs: [0.0, 0.24963094468700395, 0.7525104872191076, 1.229735980709192, 0.7385208497317045, -0.4328168377896504, -1.0125409246826416, -1.029692234224304, -0.9589242746631385]

julia> spl(π/4)
0.7071028397621081
```
"""
function approximate(f, basis::BSplineBasis; indices::Union{AbstractUnitRange,Colon}=Colon())
    xvalues = knotaverages(basis, indices=indices)
    interpolate(basis, xvalues, f.(xvalues); indices=indices)
end

"""
    interpolate(basis::BSplineBasis, xvalues, yvalues; indices=eachindex(basis)) -> Spline

Interpolate the data given by `xvalues` and `yvalues` in the B-spline basis `basis`. If
`indices` is supplied, only the basis functions at the given indices of `basis` are used.

The spline interpolation is calculated by creating the matrix
`B = [basis[i](x) for x=xvalues, i=indices]` and then calculating `B\\yvalues`.

See also: [`approximate`](@ref)

# Examples

```jldoctest
julia> basis = BSplineBasis(5, 1:10);

julia> xs = range(1, stop=10, length=length(basis)); ys = log.(xs);

julia> spl = interpolate(basis, xs, ys)
Spline{BSplineBasis{UnitRange{Int64}},Array{Float64,1}}:
 basis: 13-element BSplineBasis{UnitRange{Int64}}:
  order: 5
  breakpoints: 1:10
 coeffs: [0.0, 0.2480193113006778, 0.5968722382485209, 0.9466707785821219, 1.2689430722820556, 1.5140484163988175, 1.7114875128056504, 1.8766572742486973, 2.0185639127626325, 2.1429230073972407, 2.2259183994808813, 2.2775846911059, 2.302585092994046]

julia> spl(float(ℯ))
0.9999766059171411
```
"""
function interpolate(basis::BSplineBasis, xvalues::AbstractVector, yvalues::AbstractVector;
                     indices::Union{AbstractUnitRange,Colon}=Colon())
    if !checkbounds(Bool, basis, indices)
        throw(ArgumentError("indices must be valid indices for basis."))
    end
    if length(xvalues) != length(yvalues)
        throw(DimensionMismatch("lengths of x and y values do not match."))
    end
    bmatrix = basismatrix(basis, xvalues, indices=indices)
    _interpolate(basis, bmatrix, yvalues, indices)
end

function _interpolate(basis, basismatrix, yvalues, indices::AbstractUnitRange)
    indcoeffs = basismatrix \ yvalues
    coeffs = zeros(eltype(indcoeffs), length(basis))
    coeffs[indices] = indcoeffs
    Spline(basis, coeffs)
end

_interpolate(basis, basismatrix, yvalues, ::Colon) = Spline(basis, basismatrix\yvalues)

"""
    knotaverages(basis::BSplineBasis; indices=eachindex(basis))

Return the knot averages `τ[i] = mean(knots[i+1:i+order-1])` for `i ∈ indices` and the
knots and order of `basis`. The knot averages are recommended in [^deBoor1978] (p. 214) as
data points for interpolation.

See also: [`knotaverages!`](@ref)

[^deBoor1978]:
    Carl de Boor, *A Practical Guide to Splines*, New York, N.Y.: Springer, 1978.

# Examples

```jldoctest
julia> knotaverages(BSplineBasis(3, 0:5))
7-element Array{Float64,1}:
 0.0
 0.5
 1.5
 2.5
 3.5
 4.5
 5.0

julia> knotaverages(BSplineBasis(4, [1, 3//2, 5//2, 4]), indices=2:6)
5-element Array{Rational{Int64},1}:
 7//6
 5//3
 8//3
 7//2
 4//1
```
"""
function knotaverages(basis::BSplineBasis; indices=Colon())
    if !checkbounds(Bool, basis, indices)
        throw(ArgumentError("invalid indices for basis: $indices"))
    end
    T = bspline_returntype(basis, Bool)
    dest = Vector{T}(undef, length_indices(basis, indices))
    _knotaverages!(dest, basis, range_indices(basis, indices))
end

"""
    knotaverages!(dest, basis::BSplineBasis; indices=eachindex(basis))

Calculate the knot averages `τ[i] = mean(knots[i+1:i+order-1])` for `i ∈ indices` and the
knots and order of `basis` and store the result in `dest`. The knot averages are recommended
in [^deBoor1978] (p. 214) as data points for interpolation.

See also: [`knotaverages`](@ref)

[^deBoor1978]:
    Carl de Boor, *A Practical Guide to Splines*, New York, N.Y.: Springer, 1978.

# Examples

```jldoctest
julia> dest = Vector{Float64}(undef, 7);

julia> knotaverages!(dest, BSplineBasis(3, 0:5))
7-element Array{Float64,1}:
 0.0
 0.5
 1.5
 2.5
 3.5
 4.5
 5.0

julia> dest = Vector{Rational{Int}}(undef, 5);

julia> knotaverages!(dest, BSplineBasis(3, 0:5), indices=2:6)
5-element Array{Rational{Int64},1}:
 1//2
 3//2
 5//2
 7//2
 9//2
```
"""
function knotaverages!(dest, basis::BSplineBasis; indices=Colon())
    has_offset_axes(dest) && throw(ArgumentError("destination vector may not have offset axes."))
    if !checkbounds(Bool, basis, indices)
        throw(ArgumentError("invalid indices for basis: $indices"))
    end
    len = length_indices(basis, indices)
    if length(dest) != len
        throw(DimensionMismatch("length of destination, $(length(dest)), " * 
                                "does not match number of indices, $len."))
    end
    _knotaverages!(dest, basis, range_indices(basis, indices))
end

@propagate_inbounds function _knotaverages!(dest, basis, indices::AbstractUnitRange)
    km1 = order(basis) - 1
    moving_average!(dest, view(knots(basis), first(indices)+1:last(indices)+km1), km1)
end

@propagate_inbounds function moving_average!(dest, xs, n)
    n⁻¹ = inv(convert(eltype(dest), n))
    moving_sum = sum(xs[j] for j=1:n)
    dest[1] = n⁻¹ * moving_sum
    for i = firstindex(dest)+1:lastindex(dest)
        moving_sum -= xs[i-1]
        moving_sum += xs[i-1+n]
        dest[i] = n⁻¹ * moving_sum
    end
    dest
end

"""
    averagebasis(order, datapoints) -> BSplineBasis

Returns a B-spline basis with the specified `order` that is well-suited for interpolation on
the given `datapoints`. The `datapoints` vector is assumed to be sorted.

The calculated breakpoints are described in [^deBoor1978], p. 219, as a “reasonable
alternative” to the optimal breakpoint sequence since they are “often very close to the
optimum” and are computationally inexpensive.

[^deBoor1978]:
    Carl de Boor, *A Practical Guide to Splines*, New York, N.Y.: Springer, 1978.

# Examples

```jldoctest
julia> averagebasis(5, 0:10)
11-element BSplineBasis{Array{Float64,1}}:
 order: 5
 breakpoints: [0.0, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 10.0]
```
"""
function averagebasis(order::Integer, datapoints::AbstractVector{<:Real})
    if order > length(datapoints)
        throw(ArgumentError("order cannot not be greater than number of datapoints."))
    end
    BSplineBasis(order, averagebasis_breakpoints(order, datapoints))
end

function averagebasis_breakpoints(k, datapoints)
    T = bspline_returntype(Int, eltype(datapoints))
    bps = Vector{T}(undef, length(datapoints)-k+2)
    bps[1] = first(datapoints)
    bps[end] = last(datapoints)
    if length(bps) > 2
        moving_average!(@view(bps[2:end-1]), @view(datapoints[2:end-1]), k-1)
    end
    bps
end

end # module BSplines