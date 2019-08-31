struct KnotVector{T,P<:AbstractVector{T}} <: AbstractVector{T}
    parent::P
    front::Int
    back::Int

    @inline function KnotVector{T,P}(parent, front, back) where {T,P<:AbstractVector{T}}
        @boundscheck begin
            front ≥ 0 || throw(DomainError(front, "front padding must be non-negative."))
            back  ≥ 0 || throw(DomainError(back,  "back padding must be non-negative."))
            has_offset_axes(parent) && throw(ArgumentError("parent array must not have offset axes."))
            isempty(parent) && throw(ArgumentError("parent array must not be empty."))
        end
        new(parent, front, back)
    end
end

"""
    KnotVector(parent::AbstractVector, front[, back=front])

Create a view into `parent` where the first element of `parent` is appended `front` times at
the front and the last element `back` times at the back. (not exported)

# Examples 

```jldoctest; setup = :(using BSplines: KnotVector)
julia> KnotVector(0:4, 2)
9-element KnotVector{Int64,UnitRange{Int64}}:
 0
 0
 0
 1
 2
 3
 4
 4
 4

julia> KnotVector(["first", "second", "last"], 2, 3)
8-element KnotVector{String,Array{String,1}}:
 "first"
 "first"
 "first"
 "second"
 "last"
 "last"
 "last"
 "last"
```
"""
@propagate_inbounds KnotVector(parent::AbstractVector, front, back=front) =
    KnotVector{eltype(parent), typeof(parent)}(parent, front, back)
@propagate_inbounds KnotVector(p::KnotVector, front, back=front) =
    KnotVector(parent(p), front + p.front, back + p.back)

@inline @propagate_inbounds function Base.getindex(p::KnotVector, i::Int)
    @boundscheck checkbounds(p, i)
    @inbounds parent(p)[clamp(i-p.front, 1, length(parent(p)))]
end

Base.IndexStyle(p::KnotVector) = IndexLinear()
Base.parent(p::KnotVector) = p.parent
Base.size(p::KnotVector) = (length(parent(p)) + p.front + p.back,)

Base.first(p::KnotVector) = @inbounds first(parent(p))
Base.last(p::KnotVector) = @inbounds last(parent(p))
Base.issorted(p::KnotVector) = issorted(parent(p))

function intervalindex(p::KnotVector, x)
    index = intervalindex(parent(p), x)
    index === nothing && return nothing
    return index + p.front
end

function intervalindex(p::KnotVector, x, start::Integer)
    checkbounds(p, start)
    index = intervalindex(parent(p), x, clamp(start-p.front, 1, length(parent(p))))
    index === nothing && return nothing
    return index + p.front
end

struct StandardBasisVector{T<:Number} <: AbstractVector{T}
    length::Int
    index::Int

    @inline function StandardBasisVector{T}(length, index) where T
        @boundscheck if !(1 ≤ index ≤ length)
            throw(DomainError(index, "Index not in range 1:$length."))
        end
        new(length, index)
    end
end

"""
    StandardBasisVector([T=Bool,] length, i)

Create a vector of the specified `length` where the `i`-th element is equal to `1` and all
other elements are equal to `0`. The elements are of type `T`. (not exported)

# Examples

```jldoctest; setup = :(using BSplines: StandardBasisVector)
julia> StandardBasisVector(5, 3)
5-element StandardBasisVector{Bool}:
 0
 0
 1
 0
 0

julia> StandardBasisVector(Float64, 6, 2)
6-element StandardBasisVector{Float64}:
 0.0
 1.0
 0.0
 0.0
 0.0
 0.0
```
"""
@propagate_inbounds StandardBasisVector(length, index) = StandardBasisVector(Bool, length, index)
@propagate_inbounds StandardBasisVector(T::Type, length, index) =
    StandardBasisVector{T}(length, index)

Base.size(A::StandardBasisVector) = (A.length,)

@inline @propagate_inbounds function Base.getindex(A::StandardBasisVector, i::Int)
    @boundscheck checkbounds(A, i)
    ifelse(i == A.index, one(eltype(A)), zero(eltype(A)))
end

Base.:(==)(x::StandardBasisVector, y::StandardBasisVector) =
    (x.length == y.length) & (x.index == y.index)
