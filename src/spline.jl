"""
    Spline{B<:AbstractBSplineBasis, C<:AbstractVector{<:Real}}

Type for a spline based on a B-spline basis of type `B` and coefficient vector of type `C`.
"""
struct Spline{B<:AbstractBSplineBasis,C<:AbstractVector{<:Real}}
    basis::B
    coeffs::C

    function Spline{B,C}(basis, coeffs) where {B<:AbstractBSplineBasis,C<:AbstractVector{<:Real}}
        @boundscheck if eachindex(basis) != eachindex(coeffs)
            throw(DimensionMismatch("basis has indices $(eachindex(basis)), coefficient vector " *
                                    "has indices $(eachindex(coeffs))."))
        end
        new(basis, coeffs)
    end
end

"""
    Spline(basis::AbstractBSplineBasis, coeffs)

Create a spline from a B-spline basis and a vector of coefficients.
"""
Spline(basis, coeffs) = Spline{typeof(basis),typeof(coeffs)}(basis, coeffs)

Base.:(==)(x::Spline, y::Spline) = basis(x) == basis(y) && coeffs(x) == coeffs(y)

Base.hash(x::Spline, h::UInt) =
    hash(coeffs(x), hash(breakpoints(basis(x)), hash(order(x), hash(:Spline, h))))

function Base.show(io::IO, ::MIME"text/plain", s::Spline)
    summary(io, s); println(io, ':')
    print(io, " basis: "); summary(io, basis(s)); println(io, ':')
    println(io, "  order: ", order(basis(s)))
    println(io, "  breakpoints: ", breakpoints(basis(s)))
    print(io, " coeffs: ", coeffs(s))
end

(s::Spline)(args...; kwargs...) = splinevalue(s, args...; kwargs...)

check_intervalindex(::Type{Bool}, s::Spline, x, index) =
    check_intervalindex(Bool, basis(s), x, index)

"""
    basis(spline::Spline)

Return the B-spline basis on which the spline is defined.
"""
basis(s::Spline) = s.basis

"""
    coeffs(spline::Spline)

Return the coefficient vector of `spline`.
"""
coeffs(s::Spline) = s.coeffs

order(s::Spline) = order(basis(s))

"""
    BSpline{B} <: Spline{B}

Type for a B-spline from a B-spline basis of type `B`. `BSpline`s can be obtained by
indexing into a B-spline basis.
"""
const BSpline = Spline{T, StandardBasisVector{Bool}} where T<:AbstractBSplineBasis

@propagate_inbounds BSpline{B}(basis, index::Integer) where B<:AbstractBSplineBasis =
    BSpline{B}(basis, StandardBasisVector(length(basis), index))

"""
    BSpline(basis::AbstractBSplineBasis, index)

Return the `index`-th B-spline of `basis`.
"""
@propagate_inbounds BSpline(basis::B, index::Integer) where B<:AbstractBSplineBasis =
    BSpline{B}(basis, index)

Base.show(io::IO, x::BSpline) =
    (summary(io, x); print(io, '(', basis(x), ", ", coeffs(x).index, ')'))

function Base.show(io::IO, ::MIME"text/plain", x::BSpline)
    summary(io, x); println(io, ':')
    print(io, " basis: "); summary(io, basis(x)); println(io, ':')
    println(io, "  order: ", order(basis(x)))
    println(io, "  breakpoints: ", breakpoints(basis(x)))
    print(io, " index: ", coeffs(x).index, " (knots: ", bsplineknots(x), ')')
end

Base.summary(io::IO, x::BSpline{B}) where B = print(io, "BSpline{", B, '}')

bsplineknots(x::BSpline) = bsplineknots(basis(x), coeffs(x).index)
bsplineknots(basis, index) = knots(basis)[index:index+order(basis)]

"""
    support(spline::Spline) -> a, b

If `spline` is a `BSpline`, return the interval ``[a,b]`` on which the B-spline is non-zero.
Otherwise, return `support(basis(spline))`.

# Examples

```jldoctest
julia> basis = BSplineBasis(3, 0:5);

julia> spline = Spline(basis, ones(7));

julia> zerospline = Spline(basis, zeros(7));

julia> support(spline)
(0, 5)

julia> support(zerospline) # even though the spline is zero everywhere
(0, 5)

julia> support(basis[4]) # for BSplines, return their actual support
(1, 4)
```
"""
support(s::Spline) = support(basis(s))
support(b::BSpline) = (t = knots(basis(b)); i = coeffs(b).index; (t[i], t[i+order(b)]))

"""
    splinevalue(spline::Spline, x; leftknot=intervalindex(basis(spline), x))

Calculate the values of `spline` at `x`.

If the index of the relevant interval is already known, it can be supplied with the optional
`leftknot` keyword to speed up the calculation.

Instead of calling `splinevalue`, a spline object can be called directly:
`spline(x; [leftknot])` is equivalent to `splinevalue(spline, x; [leftknot])`.

# Examples

```jldoctest
julia> spl = Spline(BSplineBasis(4, 0:5), 1:8);

julia> splinevalue(spl, 1.7)
3.69775

julia> splinevalue(spl, 3.6, leftknot=7)
5.618

julia> spl(18//5)
2809//500

julia> spl(5, leftknot=8)
8.0
```
"""
splinevalue(::Spline, x; kwargs...)

"""
    splinevalue(spline::Spline, x, ::Derivative{N}; leftknot=intervalindex(basis(spline), x))

Calculate the value of the `N`-th derivative of `spline` at `x`.

If the index of the relevant interval is already known, it can be supplied with the optional
`leftknot` keyword to speed up the calculation.

Instead of calling `splinevalue`, a spline object can be called directly:
`spline(x, Derivative(N); [leftknot])` is equivalent to
`splinevalue(spline, x, Derivative(N); [leftknot])`.

# Examples

```jldoctest
julia> spl = Spline(BSplineBasis(4, 0:5), 1:8);

julia> splinevalue(spl, 1.7, Derivative(1))
1.0225

julia> splinevalue(spl, 18//5, Derivative(2), leftknot=7)
3//10

julia> spl(3.6, Derivative(3))
0.5

julia> spl(5, Derivative(1))
3.0
```
"""
splinevalue(::Spline, x, ::Derivative; kwargs...)
