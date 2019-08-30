"""
    AbstractBSplineBasis

Abstract supertype for B-spline bases. Here, a B-spline basis is completely specified by its
order ``k`` and breakpoint sequence. The knot sequence is derived from the breakpoint
sequence by duplicating the first and last breakpoints so they each appear ``k`` times. Knot
sequences where the first and last breakpoints do not appear ``k`` times are not supported
by this data type.
"""
abstract type AbstractBSplineBasis end

Base.:(==)(x::AbstractBSplineBasis, y::AbstractBSplineBasis) =
    order(x) == order(y) && breakpoints(x) == breakpoints(y)

Base.hash(x::AbstractBSplineBasis, h::UInt) =
    hash(breakpoints(x), hash(order(x), hash(:AbstractBSplineBasis, h)))

function Base.checkbounds(b::AbstractBSplineBasis, i)
    checkbounds(Bool, b, i) || Base.throw_boundserror(b, i)
    nothing
end
Base.checkbounds(::Type{Bool}, b::AbstractBSplineBasis, i) = checkindex(Bool, eachindex(b), i)

Base.eachindex(b::AbstractBSplineBasis) = Base.OneTo(lastindex(b))

Base.eltype(b::AbstractBSplineBasis) = BSpline{typeof(b)}

Base.length(b::AbstractBSplineBasis) = length(breakpoints(b)) + order(b) - 2

Base.iterate(b::AbstractBSplineBasis, i=1) = i-1 < length(b) ? (@inbounds b[i], i+1) : nothing

Base.firstindex(b::AbstractBSplineBasis) = 1
Base.lastindex(b::AbstractBSplineBasis) = length(b)

Base.keys(b::AbstractBSplineBasis) = eachindex(b)

@inline @propagate_inbounds function Base.getindex(b::AbstractBSplineBasis, i::Integer)
    @boundscheck checkbounds(b, i)
    @inbounds BSpline(b, i)
end

function Base.show(io::IO, ::MIME"text/plain", basis::AbstractBSplineBasis)
    summary(io, basis); println(io, ':')
    println(io, " order: ", order(basis))
    print(io, " breakpoints: ", breakpoints(basis))
end

Base.summary(io::IO, basis::AbstractBSplineBasis) =
    print(io, length(basis), "-element ", typeof(basis))

"""
    intervalindex(vec, x)
    intervalindex(vec, x, start)

If `v` is an `AbstractVector`, return the largest index `i` so that `vec[i] ≤ x` and
`vec[i] < vec[end]`. Return `nothing` if `x < first(vec)` or `x > last(vec)` or `isnan(x)`.
The vector `vec` is assumed to be sorted in ascending order.

If `vec` is a [`AbstractBSplineBasis`](@ref), return
`intervalindex(knots(vec), x)`/`intervalindex(knots(vec), x, start)`.

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

intervalindex(b::AbstractBSplineBasis, args...) = intervalindex(knots(b), args...)

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

check_intervalindex(::Type{Bool}, b::AbstractBSplineBasis, x, index) =
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
    intervalindices(basis::AbstractBSplineBasis, indices=eachindex(basis))

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
intervalindices(::AbstractBSplineBasis, ::Union{AbstractUnitRange,Colon}=Colon())

function intervalindices(basis::AbstractBSplineBasis, indices::Colon=Colon())
    bps = breakpoints(basis)
    _intervalindices(bps, eachindex(bps), order(basis)-1)
end

function intervalindices(basis::AbstractBSplineBasis, indices::AbstractUnitRange)
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
    intervalindices(basis::AbstractBSplineBasis, i, j, ...)

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
function intervalindices(basis::AbstractBSplineBasis, indices::Integer...)
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

length_indices(basis::AbstractBSplineBasis, ::Colon) = length(basis)
length_indices(basis::AbstractBSplineBasis, indices) = length(indices)

range_indices(basis::AbstractBSplineBasis, ::Colon) = eachindex(basis)
range_indices(basis::AbstractBSplineBasis, indices) = indices

"""
    breakpoints(basis::AbstractBSplineBasis)

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
function breakpoints end

"""
    knots(basis::AbstractBSplineBasis)

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
knots(b::AbstractBSplineBasis) = @inbounds KnotVector(breakpoints(b), order(b)-1)

"""
    order(spline::Spline)
    order(basis::AbstractBSplineBasis)

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

"""
    support(basis::AbstractBSplineBasis) -> a, b

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
support(b::AbstractBSplineBasis) = (bps = breakpoints(b); (first(bps), last(bps)))

"""
    basismatrix(basis::AbstractBSplineBasis, xvalues, indices)

Calculate the matrix `[basis[i](x) for x=xvalues, i=indices]`. (not exported)
"""
basismatrix(::AbstractBSplineBasis, ::Any, ::Union{AbstractUnitRange, Colon})

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
OffsetArray(::Array{Float64,1}, 3:6) with eltype Float64 with indices 3:6:
 0.03600000000000002
 0.5386666666666667 
 0.41466666666666663
 0.01066666666666666

julia> bsplines(BSplineBasis(4, 0:5), 6) # returns nothing

julia> bsplines(BSplineBasis(3, 0:5), 17//5, leftknot=6)
OffsetArray(::Array{Rational{Int64},1}, 4:6) with eltype Rational{Int64} with indices 4:6:
  9//50
 37//50
  2//25
```
"""
bsplines(::AbstractBSplineBasis, x; kwargs...)

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
OffsetArray(::Array{Float64,1}, 3:5) with eltype Float64 with indices 3:5:
 -0.6000000000000001 
  0.20000000000000018
  0.3999999999999999 

julia> bsplines(BSplineBasis(3, 0:5), 6, Derivative(1)) # returns nothing

julia> bsplines(BSplineBasis(4, 0:5), 17//5, Derivative(2), leftknot=7)
OffsetArray(::Array{Rational{Int64},1}, 4:7) with eltype Rational{Int64} with indices 4:7:
  3//5
 -4//5
 -2//5
  3//5
```
"""
bsplines(::AbstractBSplineBasis, x, ::Derivative; kwargs...)

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
OffsetArray(::Array{Float64,2}, 3:5, 0:2) with eltype Float64 with indices 3:5×0:2:
 0.18  -0.6   1.0
 0.74   0.2  -2.0
 0.08   0.4   1.0

julia> bsplines(BSplineBasis(3, 0:5), 2.4, AllDerivatives(3)) # returns nothing

julia> bsplines(BSplineBasis(4, 0:5), 17//5, AllDerivatives(4), leftknot=7)
OffsetArray(::Array{Rational{Int64},2}, 4:7, 0:3) with eltype Rational{Int64} with indices 4:7×0:3:
   9//250   -9//50   3//5  -1//1
 202//375  -14//25  -4//5   3//1
 307//750   31//50  -2//5  -7//2
   2//125    3//25   3//5   3//2
```
"""
bsplines(::AbstractBSplineBasis, x, ::AllDerivatives; kwargs...)

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
bsplines!(dest, ::AbstractBSplineBasis, x; kwargs...)

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
bsplines!(dest, ::AbstractBSplineBasis, x, ::Derivative; kwargs...)

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
bsplines!(dest, ::AbstractBSplineBasis, x, ::AllDerivatives; kwargs...)
