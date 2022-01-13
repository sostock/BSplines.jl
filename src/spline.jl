"""
    Spline{B<:BSplineBasis, C<:AbstractVector{<:Real}}

Type for a spline based on a B-spline basis of type `B` and coefficient vector of type `C`.
"""
struct Spline{B<:BSplineBasis,C<:AbstractVector{<:Real}}
    basis::B
    coeffs::C

    function Spline{B,C}(basis, coeffs) where {B<:BSplineBasis,C<:AbstractVector{<:Real}}
        @boundscheck if eachindex(basis) != eachindex(coeffs)
            throw(DimensionMismatch("basis has indices $(eachindex(basis)), coefficient vector " *
                                    "has indices $(eachindex(coeffs))."))
        end
        new(basis, coeffs)
    end
end

"""
    Spline(basis, coeffs)

Create a spline from a B-spline basis and a vector of coefficients.
"""
Spline(basis, coeffs) = Spline{typeof(basis),typeof(coeffs)}(basis, coeffs)

Base.:(==)(x::Spline, y::Spline) = basis(x) == basis(y) && coeffs(x) == coeffs(y)

Base.hash(x::Spline, h::UInt) =
    hash(coeffs(x), hash(knots(basis(x)), hash(order(x), hash(:Spline, h))))

function Base.show(io::IO, ::MIME"text/plain", s::Spline)
    cio = IOContext(io, :compact=>true)
    summary(io, s); println(io, ':')
    print(io, " basis: "); summary(io, basis(s)); println(io, ':')
    println(io, "  order: ", order(basis(s)))
    print(io, "  knots: ")
    show(cio, knots(basis(s)))
    if s isa BSpline
        print(io, "\n index: ", coeffs(s).index, " (knots: ")
        show(cio, bsplineknots(s))
        print(io, ')')
    else
        print(io, "\n coeffs: ")
        show(cio, coeffs(s))
    end
end

# Splines act as scalars for broadcasting
Base.Broadcast.broadcastable(s::Spline) = Ref(s)

for op in (:+, :-)
    @eval Base.$op(x::Spline) = Spline(basis(x), $op(coeffs(x)))
    @eval function Base.$op(x::Spline, y::Spline)
        basis(x) == basis(y) || throw(ArgumentError("splines must be defined on the same basis."))
        Spline(basis(x), $op(coeffs(x), coeffs(y)))
    end
end

Base.:*(x::Spline, y::Real) = Spline(basis(x), coeffs(x)*y)
Base.:/(x::Spline, y::Real) = Spline(basis(x), coeffs(x)/y)

Base.:*(x::Real, y::Spline) = Spline(basis(y), x*coeffs(y))
Base.:\(x::Real, y::Spline) = Spline(basis(y), x\coeffs(y))

LinearAlgebra.lmul!(x::Real, y::Spline) = (lmul!(x, coeffs(y)); y)
LinearAlgebra.rmul!(x::Spline, y::Real) = (rmul!(coeffs(x), y); x)

@static if VERSION ≥ v"1.2"
    LinearAlgebra.ldiv!(x::Real, y::Spline) = (ldiv!(x, coeffs(y)); y)
    LinearAlgebra.rdiv!(x::Spline, y::Real) = (rdiv!(coeffs(x), y); x)
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

Type for a B-spline from a B-spline basis of type `B`.
"""
const BSpline = Spline{B, StandardBasisVector{Bool}} where B<:BSplineBasis

@propagate_inbounds BSpline{B}(basis, index::Integer) where B<:BSplineBasis =
    BSpline{B}(basis, StandardBasisVector(length(basis), index))

"""
    BSpline(basis::BSplineBasis, index)

Return the `index`-th B-spline of `basis`. Instead of `BSpline(basis, index)`, the shorthand
`basis[index]` can be used.

# Example

```jldoctest
julia> basis = BSplineBasis(5, breakpoints=1:10);

julia> BSpline(basis, 3) # or: basis[3]
$(BSpline{BSplineBasis{BSplines.KnotVector{Int64, UnitRange{Int64}}}}):
 basis: 13-element $(BSplineBasis{BSplines.KnotVector{Int64, UnitRange{Int64}}}):
  order: 5
  knots: [1, 1, 1, 1, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 10, 10, 10, 10]
 index: 3 (knots: [1, 1, 1, 2, 3, 4])
```
"""
@propagate_inbounds BSpline(basis::B, index::Integer) where B<:BSplineBasis =
    BSpline{B}(basis, index)

Base.show(io::IO, x::BSpline) =
    (summary(io, x); print(io, '(', basis(x), ", ", coeffs(x).index, ')'))

@static if VERSION < v"1.6.0-DEV.347"
    Base.summary(io::IO, x::BSpline{B}) where B = print(io, "BSpline{", B, '}')
end

bsplineknots(x::BSpline) = bsplineknots(basis(x), coeffs(x).index)
bsplineknots(basis, index) = knots(basis)[index:index+order(basis)]

"""
    support(spline::Spline) -> a, b

If `spline` is a `BSpline`, return the interval ``[a,b]`` on which the B-spline is non-zero.
Otherwise, return `support(basis(spline))`.

# Examples

```jldoctest
julia> basis = BSplineBasis(3, breakpoints=0:5);

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
    splinevalue(spline::Spline, x[, ::Derivative{N}]; kwargs...)

Calculate the value of `spline`, or its `N`-th derivative, at `x`.

Two optional keyword arguments can be used to increase performance:
* `workspace`: By default, the function allocates a vector of length `order(spline)` in
  which the calculation is performed. To avoid this, a pre-allocated vector can be supplied
  with the `workspace` keyword. In this case, the returned value is always of type
  `eltype(workspace)`.
* `leftknot`: If the index of the relevant interval (i.e.,
  `intervalindex(basis(spline), x)`) is already known, it can be supplied with the
  `leftknot` keyword.

Instead of calling `splinevalue`, a spline object can be called directly:
`spline(x[, Derivative(N)]; kwargs...)` is equivalent to
`splinevalue(spline, x[, Derivative(N)]; kwargs...)`.

# Examples

```jldoctest
julia> spl = Spline(BSplineBasis(4, breakpoints=0:5), 1:8);

julia> splinevalue(spl, 1.7)
3.69775

julia> splinevalue(spl, 1.7, Derivative(1))
1.0225

julia> splinevalue(spl, 3.6, leftknot=7)
5.618

julia> spl(18//5)
2809//500

julia> spl(3.6, Derivative(3), leftknot=7)
0.5
```
"""
function splinevalue(spline::Spline, x, drv::Derivative{N}=NoDerivative();
                     leftknot=intervalindex(basis(spline), x),
                     workspace=nothing) where N
    check_intervalindex(spline, x, leftknot)
    T = workspace === nothing ? bspline_returntype(spline, x) : eltype(workspace)
    leftknot === nothing && return isnan(x) ? T(NaN) : zero(T)
    if spline isa BSpline
        iszero_bspline(spline, leftknot) && return zero(T)
    end
    N ≥ order(spline) && return zero(T)
    if workspace === nothing
        workspace = Vector{T}(undef, order(spline))
    else
        exp = (Base.OneTo(order(spline)),)
        got = axes(workspace)
        if exp != got
            throw(DimensionMismatch("workspace has wrong axes: expected $exp, got $got"))
        end
    end
    @inbounds _splinevalue(workspace, basis(spline), coeffs(spline), x, leftknot, drv)
end

# Return true iff `spline` is zero on the interval given by `leftknot`
function iszero_bspline(spline::BSpline, leftknot::Integer)
    i = coeffs(spline).index
    (i < leftknot-order(spline)+1) | (i > leftknot)
end

# The implementation of this method (and the methods iterate_splinevalue_derivatives! and
# iterate_splinevalue_bsplines! it calls) is adapted from the Fortran subroutine BVALUE from
# Carl de Boor’s book *A practical Guide to Splines* [^deBoor1978].
#
# [^deBoor1978]:
#     Carl de Boor, *A Practical Guide to Splines*, New York, N.Y.: Springer-Verlag, 1978.
@propagate_inbounds function _splinevalue(A, basis::BSplineBasis, coeffs, x, leftknot::Integer, ::Derivative{N}) where N
    t = knots(basis)
    k = order(basis)
    xtyped = convert(eltype(A), x)
    A .= @view coeffs[leftknot-k+1:leftknot]
    # Difference the coefficients (stored in A) N times
    for j = 1:N
        iterate_splinevalue_derivatives!(A, t, k, j, leftknot)
    end
    # Compute value of N-th derivative at x from its B-spline coefficients in A[1:k-N]
    for j = N+1:k-1
        iterate_splinevalue_bsplines!(A, t, k, j, xtyped, leftknot)
    end
    A[1]
end

@propagate_inbounds function iterate_splinevalue_bsplines!(coeffs::AbstractArray{T}, knots, k, j, x::T, leftknot) where T
    kmj = k-j
    for i = 1:kmj
        rindex = leftknot + i
        lindex = rindex - kmj
        tl = convert(T, knots[lindex])
        tr = convert(T, knots[rindex])
        coeffs[i] = (coeffs[i+1]*(x-tl) + coeffs[i]*(tr-x)) / (tr - tl)
    end
end

@propagate_inbounds function iterate_splinevalue_derivatives!(coeffs::AbstractArray{T}, knots, k, j, leftknot) where T
    kmj = k-j
    for i = 1:kmj
        rindex = leftknot + i
        lindex = rindex - kmj
        coeffs[i] = (coeffs[i+1] - coeffs[i]) / (convert(T, knots[rindex]) - convert(T, knots[lindex])) * kmj
    end
end
