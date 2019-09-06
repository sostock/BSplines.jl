"""
    BSplineBasis{T<:AbstractVector{<:Real}} <: AbstractBSplineBasis

Type for a B-spline basis with breakpoint vector of type `T`.
"""
struct BSplineBasis{T<:AbstractVector{<:Real}} <: AbstractBSplineBasis
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

breakpoints(b::BSplineBasis) = b.breakpoints

order(b::BSplineBasis) = b.order

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

function splinevalue(spline::Spline, x, drv::Derivative{N}=NoDerivative();
                     leftknot=intervalindex(basis(spline), x)) where N
    check_intervalindex(spline, x, leftknot)
    T = bspline_returntype(spline, x)
    leftknot === nothing && return isnan(x) ? T(NaN) : zero(T)
    if spline isa BSpline
        iszero_bspline(spline, leftknot) && return zero(T)
    end
    N ≥ order(spline) && return zero(T)
    @inbounds _splinevalue(T, basis(spline), coeffs(spline), x, leftknot, drv)
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
@propagate_inbounds function _splinevalue(T::Type, basis::BSplineBasis, coeffs, x, leftknot::Integer, ::Derivative{N}) where N
    t = knots(basis)
    k = order(basis)
    A = Vector{T}(undef, k)
    A .= @view coeffs[leftknot-k+1:leftknot]
    # Difference the coefficients (stored in A) N times
    for j = 1:N
        iterate_splinevalue_derivatives!(A, t, k, j, leftknot)
    end
    # Compute value of N-th derivative at x from its B-spline coefficients in A[1:k-N]
    for j = N+1:k-1
        iterate_splinevalue_bsplines!(A, t, k, j, x, leftknot)
    end
    A[1]
end

@propagate_inbounds function iterate_splinevalue_bsplines!(coeffs, knots, k, j, x, leftknot)
    kmj = k-j
    for i = 1:kmj
        rindex = leftknot + i
        lindex = rindex - kmj
        tl = knots[lindex]
        tr = knots[rindex]
        coeffs[i] = (coeffs[i+1]*(x-tl) + coeffs[i]*(tr-x)) / (tr - tl)
    end
end

@propagate_inbounds function iterate_splinevalue_derivatives!(coeffs, knots, k, j, leftknot)
    kmj = k-j
    for i = 1:kmj
        rindex = leftknot + i
        lindex = rindex - kmj
        coeffs[i] = (coeffs[i+1] - coeffs[i]) / (knots[rindex] - knots[lindex]) * kmj
    end
end

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
