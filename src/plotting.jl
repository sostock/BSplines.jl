import RecipesBase

using RecipesBase: @recipe

"""
    SplineFunctionWrapper{S<:Spline, D<:Derivative, OnlySupport} <: Function

A function type that wraps a spline of type `S`. (not exported)

Calling an object of this type evaluates the derivative `D` of the spline. The parameter
`OnlySupport` determines the value returned when evaluating the function at a point outside
of the [`support`](@ref) of the spline: If `OnlySupport === true`, a value equal to `NaN`
will be returned; if `OnlySupport === false`, a value equal to zero is returned.
"""
struct SplineFunctionWrapper{S<:Spline,D<:Derivative,OnlySupport} <: Function
    spline::S

    function SplineFunctionWrapper{S,D,OnlySupport}(spline) where {S<:Spline,D<:Derivative,OnlySupport}
        OnlySupport isa Bool || throw(ArgumentError("type parameter OnlySupport must be a Bool."))
        new(spline)
    end
end

SplineFunctionWrapper(s::Spline, d::Derivative, onlysupport::Bool=true) =
    SplineFunctionWrapper{typeof(s),typeof(d),onlysupport}(s)

SplineFunctionWrapper(s::Spline, onlysupport::Bool=true) =
    SplineFunctionWrapper(s, NoDerivative(), onlysupport)

function (f::SplineFunctionWrapper{S,D,true})(x) where {S,D}
    a,b = support(f.spline)
    val = f.spline(x, D())
    a ≤ x ≤ b ? val : convert(float(typeof(val)), NaN)
end
(f::SplineFunctionWrapper{S,D,false})(x) where {S,D} = f.spline(x, D())

support(s::SplineFunctionWrapper) = support(s.spline)

"""
    Function(spline::Spline, [deriv::Derivative], [onlysupport=true])

Wrap `spline` in an object that is a subtype of `Function`. Calling the returned function
with a single argument `x` will evaluate the spline (or one of its derivatives as specified
by `deriv`) at `x`.

If the optional argument `onlysupport` is set to `true` (the default), the returned function
will return `NaN` if evaluated outside of the [`support`](@ref) of `spline`. If
`onlysupport` is set to `false`, it will return zero there (as does calling the `spline`
directly).

Note that a `Spline` can be called without wrapping them as described here, although they
are not a subtype of `Function`. Wrapping a spline in a `Function` object is mainly intended
to aid in plotting, which is the rationale behind the `onlysupport=true` default: when using
[`Plots.jl`](https://github.com/JuliaPlots/Plots.jl), this will cause the spline to not be
drawn outside of its support.

# Examples

```jldoctest
julia> spline = approximate(sin, BSplineBasis(5, 0:5)); # create a Spline

julia> f = Function(spline, false); # f is zero outside of the interval [0,5]

julia> g = Function(spline, Derivative(1)); # g is NaN outside of the interval [0,5]

julia> f(1.5) == spline(1.5)
true

julia> g(1.5) == spline(1.5, Derivative(1))
true

julia> f(-1), g(-1)
(0.0, NaN)
```
"""
Base.Function(x::Spline, args...) = SplineFunctionWrapper(x, args...)

# Type recipe for wrapping a `Spline` in a `Function`
@recipe f(::Type{S}, s::S) where {S<:Spline} = Function(s, true)

# One argument
@recipe f(s::SplineFunctionWrapper) = (s, support(s)...)
@recipe f(s::Vector{S}) where {S<:Spline} = (Function.(s, true), joinsupport(s)...)
@recipe f(s::Vector{S}) where {S<:SplineFunctionWrapper} = (s, joinsupport(s)...)

joinsupport(splines) = mapreduce(support, _joinsupport, splines)
_joinsupport((xmin, xmax), (ymin, ymax)) = (min(xmin, ymin), max(xmax, ymax))

# Two arguments
@recipe f(x, s::Vector{S}) where {S<:Spline} = (x, Function.(s, true))
@recipe f(s::Vector{S}, x) where {S<:Spline} = (x, Function.(s, true))

# Three arguments
@recipe f(s::Vector{S}, xmin::Number, xmax::Number) where {S<:Spline} =
    (Function.(s, true), xmin, xmax)

# BSplineBasis
@recipe function f(basis::BSplineBasis, pointsperspline::Integer=100)
    label --> permutedims(eachindex(basis))
    xs = Matrix{Float64}(undef, pointsperspline, length(basis))
    ys = similar(xs)
    for (j, bspline) = enumerate(basis)
        xmin, xmax = support(bspline)
        for (i, x) = enumerate(range(xmin, stop=xmax, length=pointsperspline))
            xs[i,j] = x
            ys[i,j] = bspline(x)
        end
    end
    xs, ys
end
