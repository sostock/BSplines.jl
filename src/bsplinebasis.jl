"""
    BSplineBasis{T<:AbstractVector{<:Real}}

Type for a B-spline basis with breakpoint vector of type `T`.

Here, a B-spline basis is completely specified by its order ``k`` and breakpoint sequence.
The knot sequence is derived from the breakpoint sequence by duplicating the first and last
breakpoints so they each appear ``k`` times. Knot sequences where the first and last
breakpoints do not appear ``k`` times are not supported by this data type.
"""
struct BSplineBasis{T<:AbstractVector{<:Real}}
    order::Int
    breakpoints::T

    function BSplineBasis{T}(order, breakpoints) where T<:AbstractVector{<:Real}
        order ≥ 1 || throw(DomainError(order, "order of B-splines must be positive."))
        length(breakpoints) ≥ 2 || throw(ArgumentError("length of breakpoint vector must be at least 2."))
        has_offset_axes(breakpoints) && throw(ArgumentError("breakpoint vector must not have offset axes."))
        issorted(breakpoints) || throw(ArgumentError("breakpoint sequence must be non-decreasing."))
        new(order, breakpoints)
    end
end

"""
    BSplineBasis(order, breakpoints)

Create a B-spline basis with order `order` and breakpoint vector `breakpoints`.
"""
BSplineBasis(order, breakpoints) = BSplineBasis{typeof(breakpoints)}(order, breakpoints)

Base.:(==)(x::BSplineBasis, y::BSplineBasis) =
    order(x) == order(y) && breakpoints(x) == breakpoints(y)

Base.hash(x::BSplineBasis, h::UInt) =
    hash(breakpoints(x), hash(order(x), hash(:BSplineBasis, h)))

function Base.checkbounds(b::BSplineBasis, i)
    checkbounds(Bool, b, i) || Base.throw_boundserror(b, i)
    nothing
end
Base.checkbounds(::Type{Bool}, b::BSplineBasis, i) = checkindex(Bool, eachindex(b), i)

Base.eachindex(b::BSplineBasis) = Base.OneTo(lastindex(b))

Base.eltype(b::BSplineBasis) = BSpline{typeof(b)}

Base.length(b::BSplineBasis) = length(breakpoints(b)) + order(b) - 2

Base.iterate(b::BSplineBasis, i=1) = i-1 < length(b) ? (@inbounds b[i], i+1) : nothing

Base.firstindex(b::BSplineBasis) = 1
Base.lastindex(b::BSplineBasis)  = length(b)

Base.keys(b::BSplineBasis) = eachindex(b)

@propagate_inbounds function Base.getindex(b::BSplineBasis, i::Integer)
    @boundscheck checkbounds(b, i)
    @inbounds BSpline(b, i)
end

function Base.show(io::IO, ::MIME"text/plain", basis::BSplineBasis)
    summary(io, basis); println(io, ':')
    println(io, " order: ", order(basis))
    print(io, " breakpoints: ", breakpoints(basis))
end

Base.summary(io::IO, basis::BSplineBasis) =
    print(io, length(basis), "-element ", typeof(basis))

"""
    breakpoints(basis::BSplineBasis)

Return the breakpoint sequence of the B-spline basis.

# Examples

```jldoctest
julia> breakpoints(BSplineBasis(3, 0:5))
0:5

julia> breakpoints(BSplineBasis(4, [1.0, 1.5, 2.5, 4.0]))
4-element Array{Float64,1}:
 1.0
 1.5
 2.5
 4.0
```
"""
breakpoints(b::BSplineBasis) = b.breakpoints

"""
    knots(basis::BSplineBasis)

Return the knot sequence of the B-spline basis.

The knot sequence is the breakpoint sequence except that the first and last values are
duplicated so they appear `order(basis)` times.

# Examples

```jldoctest
julia> knots(BSplineBasis(3, 0:5))
10-element BSplines.KnotVector{Int64,UnitRange{Int64}}:
 0
 0
 0
 1
 2
 3
 4
 5
 5
 5
```
"""
knots(b::BSplineBasis) = @inbounds KnotVector(breakpoints(b), order(b)-1)

"""
    order(spline::Spline)
    order(basis::BSplineBasis)

Return the order of a spline or a B-spline basis.

# Examples

```jldoctest
julia> order(BSplineBasis(3, 0:5))
3

julia> order(BSplineBasis(4, [1.0, 1.5, 2.5, 4.0]))
4
```
"""
function order end

order(b::BSplineBasis) = b.order

"""
    intervalindex(vec, x[, start])

If `v` is an `AbstractVector`, return the largest index `i` so that `vec[i] ≤ x` and
`vec[i] < vec[end]`. Return `nothing` if `x < first(vec)` or `x > last(vec)` or `isnan(x)`.
The vector `vec` is assumed to be sorted in ascending order.

If `vec` is a [`BSplineBasis`](@ref), return `intervalindex(knots(vec), x[, start])`.

If `start` is given, a linear search is performed, starting from the index `start` going
forward or backward. If `start` is not given, a binary search is performed.

# Examples

```jldoctest
julia> intervalindex([1,1,2,3,4,4,4,5,6,6], 2.5)
3

julia> intervalindex([1,1,2,3,4,4,4,5,6,6], 7) # returns nothing

julia> intervalindex([1,1,2,3,4,4,4,5,6,6], 1)
2

julia> intervalindex([1,1,2,3,4,4,4,5,6,6], 4)
7

julia> intervalindex([1,1,2,3,4,4,4,5,6,6], 6.0)
8
```
"""
function intervalindex end

function intervalindex(vec::AbstractVector, x)
    isnan(x) && return nothing
    lastelement = last(vec)
    x < first(vec) && return nothing
    x > lastelement && return nothing
    if x == lastelement
        for i = lastindex(vec)-1:-1:firstindex(vec)
            vec[i] < lastelement && return Int(i)
        end
        return nothing
    end
    return Int(searchsortedlast(vec, x))
end

function intervalindex(vec::AbstractRange, x)
    isnan(x) && return nothing
    firstelement = first(vec)
    lastelement = last(vec)
    x < firstelement && return nothing
    x > lastelement && return nothing
    firstelement == lastelement && return nothing
    x == lastelement && return Int(lastindex(vec)-1)
    return Int(searchsortedlast(vec, x))
end

function intervalindex(vec::AbstractVector, x, start::Integer)
    if vec[start] ≤ x
        ind = findnext(v -> v > x, vec, start+1)
        ind !== nothing && return Int(ind-1)
        lastelement = last(vec)
        if x == lastelement
            for i = lastindex(vec)-1:-1:firstindex(vec)
                vec[i] < lastelement && return Int(i)
            end
        end
        return nothing
    else
        ind = findprev(v -> v ≤ x, vec, start-1)
        ind === nothing && return nothing
        return Int(ind)
    end
end

function intervalindex(vec::AbstractRange, x, start::Integer)
    if vec[start] ≤ x
        ind = findnext(v -> v > x, vec, start+1)
        ind !== nothing && return Int(ind-1)
        lastelement = last(vec)
        first(vec) == lastelement && return nothing
        x == lastelement && return Int(lastindex(vec)-1)
        return nothing
    else
        ind = findprev(v -> v ≤ x, vec, start-1)
        ind === nothing && return nothing
        return Int(ind)
    end
end

intervalindex(b::BSplineBasis, args...) = intervalindex(knots(b), args...)

"""
    check_intervalindex(knots, x, index)

Throw an `ArgumentError` if `index` is not equal to `intervalindex(knots, x)`. (not
exported)
"""
function check_intervalindex(knots, x, index)
    if !check_intervalindex(Bool, knots, x, index)
        throw(ArgumentError("wrong intervalindex: $(repr(index))"))
    end
end

"""
    check_intervalindex(Bool, knots, x, index)

Return `true` if `index` is equal to `intervalindex(knots, x)`. (not exported)
"""
check_intervalindex(::Type{Bool}, knots, x, index)

check_intervalindex(::Type{Bool}, knots::AbstractVector, x, ::Nothing) =
    (x < first(knots)) | (x > last(knots)) | isnan(x)

function check_intervalindex(::Type{Bool}, knots::AbstractVector, x, index::Integer)
    checkbounds(Bool, knots, index)   || return false
    checkbounds(Bool, knots, index+1) || return false
    tl = knots[index]
    tr = knots[index+1]
    (tl ≤ x) & (tl < tr) & ((x < tr) | (x == tr == last(knots)))
end

check_intervalindex(::Type{Bool}, b::BSplineBasis, x, index) =
    check_intervalindex(Bool, knots(b), x, index)

struct IntervalIndices{T<:AbstractVector{<:Real}}
    vec::T
    indices::UnitRange{Int}
    offset::Int
end

"""
    IntervalIndices(vec, indices, offset)

Return an iterator that yields the numbers
`(i + offset for i = indices[1:end-1] if vec[i] < vec[i+1])`. (not exported)
"""
IntervalIndices(vec::T, indices, offset) where T<:AbstractVector{<:Real} =
    IntervalIndices{T}(vec, indices, offset)

Base.IteratorSize(::IntervalIndices) = Base.SizeUnknown()

Base.eltype(i::IntervalIndices) = Int

function Base.iterate(i::IntervalIndices, (index,value)=(first(i.indices),i.vec[first(i.indices)]))
    while index < last(i.indices)
        nextindex = index + 1
        nextvalue = i.vec[nextindex]
        value < nextvalue && return (index+i.offset, (nextindex, nextvalue))
        index = nextindex
    end
    return nothing
end

"""
    intervalindices(basis::BSplineBasis, indices=eachindex(basis))

Return an iterator that yields the indices of all intervals on which `basis` is defined,
i.e., it produces all indices `ind` (in ascending order) for which
`(knots(basis)[ind], knots(basis)[ind+1])` is such an interval. 

If a range of `indices` is supplied, the iterator yields only those intervals on which *at
least one* of the B-splines `basis[j] for j=indices` is non-zero.

# Examples

```jldoctest
julia> intervalindices(BSplineBasis(3, 0:5))
3:7

julia> intervalindices(BSplineBasis(3, 0:5), 1:4)
3:6

julia> intervalindices(BSplineBasis(4, [1,2,3,4,4,4,5,6]))
BSplines.IntervalIndices{Array{Int64,1}}([1, 2, 3, 4, 4, 4, 5, 6], 1:8, 3)

julia> collect(ans)
5-element Array{Int64,1}:
  4
  5
  6
  9
 10
```
"""
intervalindices(::BSplineBasis, ::Union{AbstractUnitRange,Colon}=Colon())

function intervalindices(basis::BSplineBasis, indices::Colon=Colon())
    bps = breakpoints(basis)
    _intervalindices(bps, eachindex(bps), order(basis)-1)
end

function intervalindices(basis::BSplineBasis, indices::AbstractUnitRange)
    checkbounds(basis, indices)
    bps = breakpoints(basis)
    km1 = order(basis)-1
    if isempty(indices)
        _intervalindices(bps, 1:0, km1)
    else
        firstknotindex = max(first(indices)-km1, firstindex(bps))
        lastknotindex  = min(last(indices)+1, lastindex(bps))
        _intervalindices(bps, firstknotindex:lastknotindex, km1)
    end
end

"""
    intervalindices(basis::BSplineBasis, i, j, ...)

For integers `i`, `j`, …, return an iterator that yields the indices of all intervals on
which *all* of the B-splines `basis[i]`, `basis[j]`, … are non-zero, i.e., it produces all
indices `ind` (in ascending order) for which `(knots(basis)[ind], knots(basis)[ind+1])` is
such an interval.

# Examples

```jldoctest
julia> intervalindices(BSplineBasis(3, 0:5), 3)
3:5

julia> intervalindices(BSplineBasis(3, 0:5), 4, 5)
5:6

julia> intervalindices(BSplineBasis(3, 0:5), 2, 6) # B-splines do not overlap
6:5

julia> intervalindices(BSplineBasis(3, 0:5), 3, 5, 4)
5:5

julia> intervalindices(BSplineBasis(4, [1,2,3,4,4,4,5,6]), 3, 5)
BSplines.IntervalIndices{Array{Int64,1}}([1, 2, 3, 4, 4, 4, 5, 6], 2:4, 3)

julia> collect(ans)
2-element Array{Int64,1}:
 5
 6
```
"""
function intervalindices(basis::BSplineBasis, indices::Integer...)
    for i in indices
        checkbounds(basis, i)
    end
    bps = breakpoints(basis)
    km1 = order(basis)-1
    firstknotindex = max(max(indices...)-km1, firstindex(bps))
    lastknotindex  = min(min(indices...)+1, lastindex(bps))
    _intervalindices(bps, firstknotindex:lastknotindex, km1)
end

_intervalindices(vec, indices, offset) = IntervalIndices(vec, indices, offset)

_intervalindices(vec::AbstractRange, indices, offset) =
    (iszero(step(vec)) ? indices[1:0] : indices[1:end-1]) .+ offset

length_indices(basis::BSplineBasis, ::Colon) = length(basis)
length_indices(basis::BSplineBasis, indices) = length(indices)

range_indices(basis::BSplineBasis, ::Colon) = eachindex(basis)
range_indices(basis::BSplineBasis, indices) = indices

"""
    support(basis::BSplineBasis) -> a, b

Return the interval ``[a,b]`` on which the B-spline basis is defined, i.e., `a` is the first
and `b` the last breakpoint of the `basis`.

# Examples

```jldoctest
julia> support(BSplineBasis(3, 0:5))
(0, 5)

julia> support(BSplineBasis(4, [1.0, 1.5, 2.5, 4.0]))
(1.0, 4.0)
```
"""
support(b::BSplineBasis) = (bps = breakpoints(b); (first(bps), last(bps)))

"""
    bsplines(basis, x; leftknot=intervalindex(basis, x))

Calculate the values of all non-zero B-splines of `basis` at `x`.

If any B-splines are non-zero at `x`, an `OffsetVector` is returned that contains the value
of the `i`-th B-spline at the index `i`. If no B-splines are non-zero at `x`, `nothing` is
returned.

If the index of the relevant interval is already known, it can be supplied with the optional
`leftknot` keyword to speed up the calculation.

# Examples

```jldoctest
julia> bsplines(BSplineBasis(4, 0:5), 2.4)
4-element OffsetArray(::Array{Float64,1}, 3:6) with eltype Float64 with indices 3:6:
 0.03600000000000002
 0.5386666666666667 
 0.41466666666666663
 0.01066666666666666

julia> bsplines(BSplineBasis(4, 0:5), 6) # returns nothing

julia> bsplines(BSplineBasis(3, 0:5), 17//5, leftknot=6)
3-element OffsetArray(::Array{Rational{Int64},1}, 4:6) with eltype Rational{Int64} with indices 4:6:
  9//50
 37//50
  2//25
```
"""
bsplines(basis, x; kwargs...)

"""
    bsplines(basis, x, ::Derivative{N}; leftknot=intervalindex(basis, x))

Calculate the `N`-th derivatives of all B-splines of `basis` that are non-zero at `x`.

If any B-splines are non-zero at `x`, an `OffsetVector` is returned that contains the `N`-th
derivative of the `i`-th B-spline at the index `i`. If no B-splines are non-zero at `x`,
`nothing` is returned.

If the index of the relevant interval is already known, it can be supplied with the optional
`leftknot` keyword to speed up the calculation.

# Examples

```jldoctest
julia> bsplines(BSplineBasis(3, 0:5), 2.4, Derivative(1))
3-element OffsetArray(::Array{Float64,1}, 3:5) with eltype Float64 with indices 3:5:
 -0.6000000000000001 
  0.20000000000000018
  0.3999999999999999 

julia> bsplines(BSplineBasis(3, 0:5), 6, Derivative(1)) # returns nothing

julia> bsplines(BSplineBasis(4, 0:5), 17//5, Derivative(2), leftknot=7)
4-element OffsetArray(::Array{Rational{Int64},1}, 4:7) with eltype Rational{Int64} with indices 4:7:
  3//5
 -4//5
 -2//5
  3//5
```
"""
bsplines(basis, x, ::Derivative; kwargs...)

"""
    bsplines(basis, x, ::AllDerivatives{N}; leftknot=intervalindex(basis, x))

Calculate all `m`-th derivatives (`0 ≤ m < N`) of all B-splines of `basis` that are non-zero
at `x`.

If any B-splines are non-zero at `x`, an `OffsetMatrix` is returned that contains the `m`-th
derivative of the `i`-th B-spline at the index `i, m`. If no B-splines are non-zero at `x`,
`nothing` is returned.

If the index of the relevant interval is already known, it can be supplied with the optional
`leftknot` keyword to speed up the calculation.

# Examples

```jldoctest
julia> bsplines(BSplineBasis(3, 0:5), 2.4, AllDerivatives(3))
3×3 OffsetArray(::Array{Float64,2}, 3:5, 0:2) with eltype Float64 with indices 3:5×0:2:
 0.18  -0.6   1.0
 0.74   0.2  -2.0
 0.08   0.4   1.0

julia> bsplines(BSplineBasis(3, 0:5), 6.0, AllDerivatives(3)) # returns nothing

julia> bsplines(BSplineBasis(4, 0:5), 17//5, AllDerivatives(4), leftknot=7)
4×4 OffsetArray(::Array{Rational{Int64},2}, 4:7, 0:3) with eltype Rational{Int64} with indices 4:7×0:3:
   9//250   -9//50   3//5  -1//1
 202//375  -14//25  -4//5   3//1
 307//750   31//50  -2//5  -7//2
   2//125    3//25   3//5   3//2
```
"""
bsplines(basis, x, ::AllDerivatives; kwargs...)

"""
    bsplines!(dest, basis, x; leftknot=intervalindex(basis, x)) -> offset

Calculate the values of all B-splines of `basis` that are non-zero at `x` and store the
result in `dest`. The destination vector `dest` must have the length `order(basis)`.

If any B-splines are non-zero at `x`, an integer `offset` is returned and the value of the
`i`-th B-spline is written to `dest[i-offset]`. If no B-splines are non-zero at `x`,
`nothing` is returned and `dest` is not mutated.

If the index of the relevant interval is already known, it can be supplied with the optional
`leftknot` keyword to speed up the calculation.

# Examples

```jldoctest
julia> dest = zeros(4);

julia> bsplines!(dest, BSplineBasis(4, 0:5), 2.4)
2

julia> dest # dest[i] contains value of (i+2)-th B-spline
4-element Array{Float64,1}:
 0.03600000000000002
 0.5386666666666667 
 0.41466666666666663
 0.01066666666666666
```
"""
bsplines!(dest, basis, x; kwargs...)

"""
    bsplines!(dest, basis, x, ::Derivative{N}; leftknot=intervalindex(basis, x))

Calculate the values of all B-splines of `basis` that are non-zero at `x` and store the
result in `dest`. The destination vector `dest` must have the length `order(basis)`.

If any B-splines are non-zero at `x`, an integer `offset` is returned and the `N`-th
derivative of the `i`-th B-spline is written to `dest[i-offset]`. If no B-splines are
non-zero at `x`, `nothing` is returned and `dest` is not mutated.

If the index of the relevant interval is already known, it can be supplied with the optional
`leftknot` keyword to speed up the calculation.

# Examples

```jldoctest
julia> dest = zeros(4);

julia> bsplines!(dest, BSplineBasis(4, 0:5), 7.0, Derivative(2)) # returns nothing

julia> bsplines!(dest, BSplineBasis(4, 0:5), 4.2, Derivative(2))
4

julia> dest # dest[i] contains 2nd derivative of (i+4)-th B-spline
4-element Array{Float64,1}:
  0.7999999999999998
 -1.399999999999999 
 -0.6000000000000019
  1.200000000000001 
```
"""
bsplines!(dest, basis, x, ::Derivative; kwargs...)

"""
    bsplines!(dest, basis, x, ::AllDerivatives{N}; leftknot=intervalindex(basis, x))

Calculate the values of all B-splines of `basis` that are non-zero at `x` and store the
result in `dest`. The destination matrix `dest` must have the dimensions `order(basis)`×`N`.

If any B-splines are non-zero at `x`, an integer `offset` is returned and the `m`-th
derivative of the `i`-th B-spline is written to `dest[i-offset, m+1]`. If no B-splines are
non-zero at `x`, `nothing` is returned and `dest` is not mutated.

If the index of the relevant interval is already known, it can be supplied with the optional
`leftknot` keyword to speed up the calculation.

# Examples

```jldoctest
julia> dest = zeros(4, 3);

julia> bsplines!(dest, BSplineBasis(4, 0:5), -1.0, AllDerivatives(3)) # returns nothing

julia> bsplines!(dest, BSplineBasis(4, 0:5), 3.75, AllDerivatives(3))
3

julia> dest # dest[i,m] contains (m-1)-th derivative of (i+3)-th B-spline
4×3 Array{Float64,2}:
 0.00260417  -0.03125    0.25 
 0.315104    -0.65625    0.25 
 0.576823     0.265625  -1.625
 0.105469     0.421875   1.125
```
"""
bsplines!(dest, basis, x, ::AllDerivatives; kwargs...)

function bsplines(basis::BSplineBasis, x, drv=NoDerivative(); leftknot=intervalindex(basis, x))
    check_intervalindex(basis, x, leftknot)
    leftknot === nothing && return nothing
    dest = bsplines_destarray(basis, x, drv)
    offset = @inbounds _bsplines!(dest, basis, x, leftknot, drv)
    bsplines_offsetarray(dest, offset, drv)
end

bsplines_offsetarray(arr, offset, ::Derivative)     = OffsetArray(arr, offset)
bsplines_offsetarray(arr, offset, ::AllDerivatives) = OffsetArray(arr, offset, -1)

function bsplines!(dest, basis::BSplineBasis, x, drv=NoDerivative(); leftknot=intervalindex(basis, x))
    check_intervalindex(basis, x, leftknot)
    leftknot === nothing && return nothing
    check_destarray_axes(dest, basis, drv)
    @inbounds _bsplines!(dest, basis, x, leftknot, drv)
end

bsplines_destarray(basis, x, ::Derivative) =
    Vector{bspline_returntype(basis, x)}(undef, order(basis))
bsplines_destarray(basis, x, ::AllDerivatives{N}) where N =
    Matrix{bspline_returntype(basis, x)}(undef, order(basis), N)

function check_destarray_axes(dest, basis, drv)
    got = axes(dest)
    exp = destarray_axes(basis, drv)
    if got != exp
        throw(DimensionMismatch("destination array has wrong axes: expected $exp, got $got"))
    end
    nothing
end

destarray_axes(basis, ::Derivative)                = (Base.OneTo(order(basis)),)
destarray_axes(basis, ::AllDerivatives{N}) where N = (Base.OneTo(order(basis)), Base.OneTo(N))

# The implementation of this method (and the method iterate_bsplines! it calls) is adapted
# from the Fortran subroutine BSPLVB from Carl de Boor’s book *A practical Guide to
# Splines* [^deBoor1978].
#
# [^deBoor1978]:
#     Carl de Boor, *A Practical Guide to Splines*, New York, N.Y.: Springer-Verlag, 1978.
@propagate_inbounds function _bsplines!(dest, basis::BSplineBasis, x, leftknot::Integer, ::NoDerivUnion)
    t = knots(basis)
    k = order(basis)
    dest[1] = one(eltype(dest))
    for j = 1:k-1
        iterate_bsplines!(dest, dest, t, j, x, leftknot)
    end
    return leftknot-k
end

# The implementation of this method (and the methods iterate_bsplines! and
# iterate_derivatives! it calls) is adapted from the Fortran subroutine BSPLVD from Carl de
# Boor’s book *A practical Guide to Splines* [^deBoor1978].
#
# [^deBoor1978]:
#     Carl de Boor, *A Practical Guide to Splines*, New York, N.Y.: Springer-Verlag, 1978.
@propagate_inbounds function _bsplines!(dest, basis::BSplineBasis, x, leftknot::Integer, ::AllDerivatives{N}) where N
    t = knots(basis)
    k = order(basis)
    N′ = min(k, N)
    # Set k-th and higher derivatives to zero
    dest[:, N′+1:N] .= zero(eltype(dest))
    # Successively calculate values of B-splines and store them in dest
    lastcol = @view dest[N′:k, N′]
    lastcol[1] = one(eltype(dest))
    for j = 1:k-N′
        iterate_bsplines!(lastcol, lastcol, t, j, x, leftknot)
    end
    newcol = lastcol
    for j = k-N′+1:k-1
        oldcol = newcol
        newcol = @view dest[k-j:k, k-j]
        iterate_bsplines!(newcol, oldcol, t, j, x, leftknot)
    end
    # Successively calculate coefficients for derivatives and combine them with B-splines
    drvcoeffs = Matrix{eltype(dest)}(I, k, k)
    for col = 2:N′
        # Column `col` contains `col-1`-th derivatives
        iterate_derivatives!(drvcoeffs, t, k, col-1, leftknot)
        for row = 1:k
            sum = zero(eltype(dest))
            for j = max(row,col):k
                sum += oftype(sum, drvcoeffs[j, row] * dest[j, col])
            end
            dest[row, col] = sum
        end
    end
    return leftknot-k
end

# The implementation of this method (and the methods iterate_bsplines! and
# iterate_derivatives! it calls) is adapted from the Fortran subroutine BSPLVD from Carl de
# Boor’s book *A practical Guide to Splines* [^deBoor1978].
#
# [^deBoor1978]:
#     Carl de Boor, *A Practical Guide to Splines*, New York, N.Y.: Springer-Verlag, 1978.
@propagate_inbounds function _bsplines!(dest, basis::BSplineBasis, x, leftknot::Integer, ::Derivative{N}) where N
    t = knots(basis)
    k = order(basis)
    if N ≥ k
        # k-th and higher derivatives are zero
        dest .= zero(eltype(dest))
    else
        # Calculate B-splines of order k-N and store them in dest[N+1:k]
        col = @view dest[N+1:k]
        # col = uview(dest, N+1:k)
        col[1] = one(eltype(dest))
        for j = 1:k-N-1
            iterate_bsplines!(col, col, t, j, x, leftknot)
        end
        # Calculate coefficients for derivatives and combine them with B-splines
        drvcoeffs = Matrix{eltype(dest)}(I, k, k)
        for drv = 1:N
            iterate_derivatives!(drvcoeffs, t, k, drv, leftknot)
        end
        # Combine coefficients with B-splines
        for i = 1:k
            sum = zero(eltype(dest))
            for j = max(i,N+1):k
                sum += oftype(sum, drvcoeffs[j,i] * dest[j])
            end
            dest[i] = sum
        end
    end
    return leftknot-k
end

# Take the B-splines of order `j` from `src[1:j]` and write the B-splines of order `j+1` to `dest[1:j+1]`
@propagate_inbounds function iterate_bsplines!(dest, src, knots, j, x, leftknot)
    saved = zero(eltype(dest))
    for r = 1:j
        tₗ₊ᵣ   = knots[leftknot+r]
        tₗ₊ᵣ₋ⱼ = knots[leftknot+r-j]
        term = src[r] / (tₗ₊ᵣ - tₗ₊ᵣ₋ⱼ)
        dest[r] = saved + (tₗ₊ᵣ - x) * term
        saved = oftype(saved, (x - tₗ₊ᵣ₋ⱼ) * term)
    end
    dest[j+1] = saved
end

# Iterate the coefficients in `coeffs` from the `deriv-1`-th to the `deriv`-th derivative
@propagate_inbounds function iterate_derivatives!(coeffs, knots, k, deriv, leftknot)
    kmd = k-deriv
    for j = 1:kmd
        lp1mj = leftknot+1-j
        factor = kmd / (knots[lp1mj+kmd] - knots[lp1mj])
        row = k+1-j
        for col = 1:row
            coeffs[row, col] = (coeffs[row, col] - coeffs[row-1, col]) * factor
        end
    end
end

"""
    basismatrix(basis::BSplineBasis, xvalues, indices)

Calculate the matrix `[basis[i](x) for x=xvalues, i=indices]`. (not exported)
"""
function basismatrix(basis::BSplineBasis, xvalues, indices)
    T = bspline_returntype(basis, eltype(xvalues))
    dest = zeros(T, length(xvalues), length_indices(basis, indices))
    workspace = similar(dest, order(basis))
    _basismatrix!(dest, workspace, basis, xvalues, indices)
    dest
end

function _basismatrix!(dest, workspace, basis::BSplineBasis, xvalues, indices::AbstractUnitRange)
    indicesoffset = first(indices) - 1
    start = 1
    for (xindex, x) = enumerate(xvalues)
        leftknot = intervalindex(basis, x, start)
        @assert leftknot !== nothing "xvalue outside of the support of basis: $x."
        start = leftknot
        boffset = bsplines!(workspace, basis, x, leftknot=leftknot)
        bindices = boffset .+ axes(workspace, 1)
        for bindex = bindices ∩ indices
            dest[xindex, bindex-indicesoffset] = workspace[bindex-boffset]
        end
    end
end

function _basismatrix!(dest, workspace, basis::BSplineBasis, xvalues, ::Colon)
    start = 1
    for (xindex, x) = enumerate(xvalues)
        leftknot = intervalindex(basis, x, start)
        @assert leftknot !== nothing "xvalue outside of the support of basis: $x."
        start = leftknot
        boffset = bsplines!(workspace, basis, x, leftknot=leftknot)
        for index = axes(workspace, 1)
            dest[xindex, index+boffset] = workspace[index]
        end
    end
end
