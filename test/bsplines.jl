bsplines_exact(funs, x::Irrational, indices) = bsplines_exact(funs, big(x), indices)
bsplines_exact(funs, x::Irrational, indices, drv::AllDerivatives{N}) where N =
    bsplines_exact(funs, big(x), indices, drv)

function bsplines_exact(funs, x, indices)
    bigx = maybebig(x)
    OffsetArray([funs[i,1](bigx) for i=indices], indices)
end
function bsplines_exact(funs, x, indices, ::AllDerivatives{N}) where N
    bigx = maybebig(x)
    k = length(indices)
    k == size(funs, 2) || error("got $k indices, expected $(size(funs,2))")
    arr = [funs[i,d](bigx) for i=indices, d=1:min(k,N)]
    arr = hcat(arr, zeros(eltype(arr), k, max(0,N-k)))
    OffsetArray(arr, indices, 0:N-1)
end

# Enlarge x for calculating exact B-splines:
# * for Integer/Rational: use at least 64-bit to avoid overflow in rational arithmetic
# * for floating-point: use BigFloat to minimize error
maybebig(x::Union{Integer,Rational}) = x*one(Int64)
maybebig(x::AbstractFloat) = big(x)

@time @testset "bsplines" begin
    @testset "Return types" begin
        # b1 = BSplineBasis(1, 0:5)
        b1_types = [(Int,Float64), (Float64,Float64), (Float32,Float32),
                    (Rational{Int},Rational{Int}), (BigInt,BigFloat),
                    (BigFloat,BigFloat), (Rational{BigInt},Rational{BigInt}),
                    (Rational{Int16},Rational{Int})]
        # b2 = BSplineBasis(2, big(0):1:big(5))
        b2_types = [(Int,BigFloat), (Float32,BigFloat), (Rational{Int},Rational{BigInt})]
        # b3 = BSplineBasis(3, 0//1:5//1)
        b3_types = [(Int,Rational{Int}), (Float64,Float64), (BigInt,Rational{BigInt})]
        # b4 = BSplineBasis(4, Int128[0,1,2,3,4,5])
        b4_types = [(Int,Float64), (Float64,Float64), (Rational{Int},Rational{Int128})]
        # b5 = BSplineBasis(3, [-5.0, -3.5, -2.2, 1.0, 2.0, 3.6])
        b5_types = [(Int,Float64), (Float32,Float64), (Rational{Int},Float64)]
        # b6 = BSplineBasis(4, Rational{Int64}[-5//1, -7//2, -11//5, 1//1, 2//1, 18//5])
        b6_types = [(Int,Rational{Int64}), (Float32,Float32),
                    (Rational{Int},Rational{Int64})]
        # b7 = BSplineBasis(3, Int8[-6,-4,-2,0,0,1,2])
        b7_types = [(Int,Float64), (Float64,Float64), (Rational{Int16},Rational{Int16})]
        for (basis, types) = [(b1, b1_types), (b2, b2_types), (b3, b3_types),
                              (b4, b4_types), (b5, b5_types), (b6, b6_types),
                              (b7, b7_types)]
            leftknot = intervalindex(basis, 1)
            for (xtype, bspline_type) in types
                @test eltype(bsplines(basis, one(xtype), leftknot=leftknot)) == bspline_type
                @test eltype(bsplines(basis, one(xtype), AllDerivatives(2), leftknot=leftknot)) == bspline_type
                @test eltype(bsplines(basis, one(xtype), Derivative(0), leftknot=leftknot)) == bspline_type
            end
        end

        function test_bsplines_nothing(basis, x)
            @test bsplines(basis, x)                                      === nothing
            @test bsplines(basis, x, leftknot=nothing)                    === nothing
            @test bsplines(basis, x, AllDerivatives(2))                   === nothing
            @test bsplines(basis, x, AllDerivatives(2), leftknot=nothing) === nothing
            @test bsplines(basis, x, Derivative(2))                       === nothing
            @test bsplines(basis, x, Derivative(2), leftknot=nothing)     === nothing
        end
        test_bsplines_nothing(b1, NaN16)
        test_bsplines_nothing(b2, -Inf32)
        test_bsplines_nothing(b3, Inf64)
        test_bsplines_nothing(b4, 1//0)
        test_bsplines_nothing(b5, -1//0)
        test_bsplines_nothing(b6, -100)
        test_bsplines_nothing(b7, 10.0)
    end

    function test_bsplines(basis, basis_exact, x, isexact)
        NDERIV = 4
        leftknot = intervalindex(basis, x)
        range = leftknot-order(basis)+1:leftknot
        # AllDerivatives(NDERIVS)
        bspl_nderiv = bsplines(basis, x, AllDerivatives(NDERIV), leftknot=leftknot)
        bspl_nderiv_exact = bsplines_exact(basis_exact, x, range, AllDerivatives(NDERIV))
        if isexact
            @test bspl_nderiv == bspl_nderiv_exact
        else
            @test bspl_nderiv ≈ₑₗ bspl_nderiv_exact
        end
        # No Derivative
        @test bsplines(basis, x, leftknot=leftknot) == bspl_nderiv[:,0]
        # AllDerivatives(N) and Derivative(N) for N = 1:5
        for N = 1:(NDERIV-1)
            # Derivative(N)
            @test bsplines(basis, x, Derivative(N), leftknot=leftknot) == bspl_nderiv[:, N]
            # AllDerivatives(N)
            @test bsplines(basis, x, AllDerivatives(N), leftknot=leftknot) == OffsetArray(bspl_nderiv[range, 0:N-1], range, 0:N-1)
        end
    end

    # b1 = BSplineBasis(1, 0:5)
    for x = Real[0, 1//2, 1.0f0, 4.1, Int8(5)]
        test_bsplines(b1, b1_exact, x, true)
    end
    # b2 = BSplineBasis(2, big(0):1:big(5))
    for x = Real[1.625, 2.6f0, big(5)]
        test_bsplines(b2, b2_exact, x, false)
    end
    for x = Real[22//7, big(9//2)]
        test_bsplines(b2, b2_exact, x, true)
    end
    # b3 = BSplineBasis(3, 0//1:5//1)
    for x = Real[0.123, 1.625f0, π, Float16(5.0)]
        test_bsplines(b3, b3_exact, x, false)
    end
    for x = Real[false, 1, 13//5]
        test_bsplines(b3, b3_exact, x, true)
    end
    # b4 = BSplineBasis(4, Int128[0,1,2,3,4,5])
    for x = Real[0, 2.6f0, 4.5, Int8(5)]
        test_bsplines(b4, b4_exact, x, false)
    end
    for x = Real[13//8]
        test_bsplines(b4, b4_exact, x, true)
    end
    # b5 = BSplineBasis(3, [-5.0, -3.5, -2.2, 1.0, 2.0, 3.6])
    for x = Real[-5, -4//1, -3.0f0, π, 3.6]
        test_bsplines(b5, b5_exact, x, false)
    end
    # b6 = BSplineBasis(4, Rational{Int64}[-5//1, -7//2, -11//5, 1//1, 2//1, 18//5])
    for x = Real[-3.0f0, 1.5]
        test_bsplines(b6, b6_exact, x, false)
    end
    for x = Real[-5, Int128(-2), Int8(-3), 18//5]
        test_bsplines(b6, b6_exact, x, true)
    end
    # b7 = BSplineBasis(3, Int8[-6,-4,-2,0,0,1,2])
    for x = Real[-6.0, -2.6f0, big(0.375), Int8(2)]
        test_bsplines(b7, b7_exact, x, false)
    end
    for x = Real[Rational{Int16}(0), 5//7]
        test_bsplines(b7, b7_exact, x, true)
    end
    leftknot = intervalindex(b5, 0)
    @test_throws ArgumentError bsplines(b5, 0, leftknot=leftknot-1)
    @test_throws ArgumentError bsplines(b5, 0, leftknot=leftknot+1)
    @test_throws ArgumentError bsplines(b5, 0, leftknot=nothing)
    @test_throws ArgumentError bsplines(b5, 0, Derivative(0), leftknot=leftknot-1)
    @test_throws ArgumentError bsplines(b5, 0, Derivative(1), leftknot=leftknot+1)
    @test_throws ArgumentError bsplines(b5, 0, Derivative(2), leftknot=nothing)
    @test_throws ArgumentError bsplines(b5, 0, AllDerivatives(1), leftknot=leftknot-1)
    @test_throws ArgumentError bsplines(b5, 0, AllDerivatives(2), leftknot=leftknot+1)
    @test_throws ArgumentError bsplines(b5, 0, AllDerivatives(3), leftknot=nothing)
    @test bsplines(b5, 0) == bsplines(b5, 0, leftknot=leftknot) == bsplines(b5, 0, Derivative(0), leftknot=leftknot)
    @test bsplines(b5, 0, Derivative(1)) == bsplines(b5, 0, Derivative(1), leftknot=leftknot)
    @test bsplines(b5, 0, AllDerivatives(3)) == bsplines(b5, 0, AllDerivatives(3), leftknot=leftknot)
end

@time @testset "bsplines!" begin
    function test_bsplines!_nothing(basis, x)
        dest     = Array{Float64}(undef, order(basis))
        dest_all = Array{Float32}(undef, order(basis), 2)
        @test bsplines!(dest, basis, x) === nothing
        @test bsplines!(dest, basis, x, leftknot=nothing) === nothing
        @test bsplines!(dest_all, basis, x, AllDerivatives(2)) === nothing
        @test bsplines!(dest_all, basis, x, AllDerivatives(2), leftknot=nothing) === nothing
        @test bsplines!(dest, basis, x, Derivative(2)) === nothing
        @test bsplines!(dest, basis, x, Derivative(2), leftknot=nothing) === nothing
    end
    test_bsplines!_nothing(b1, NaN16)
    test_bsplines!_nothing(b2, -Inf32)
    test_bsplines!_nothing(b3, Inf64)
    test_bsplines!_nothing(b4, -100)
    test_bsplines!_nothing(b5, 10.0)

    function test_bsplines!(basis, basis_exact, T, x, isexact)
        NDERIV = 6
        leftknot = intervalindex(basis, x)
        offset = leftknot-order(basis)
        range = offset+1:leftknot

        # AllDerivatives(NDERIV)
        bspl_nderiv = fill(T(Inf), order(basis), NDERIV)
        bspl_nderiv_exact = parent(bsplines_exact(basis_exact, x, range, AllDerivatives(NDERIV)))
        off = bsplines!(bspl_nderiv, basis, x, AllDerivatives(NDERIV), leftknot=leftknot)
        if isexact
            @test off == offset && bspl_nderiv == bspl_nderiv_exact
        else
            @test off == offset && bspl_nderiv ≈ₑₗ bspl_nderiv_exact
        end

        # No derivative
        dest = fill(T(Inf), order(basis))
        @test bsplines!(dest, basis, x, leftknot=leftknot) == offset && dest == bspl_nderiv[:,1]

        # AllDerivatives(N) and Derivative(N) for N = 1:NDERIV-1
        for N = 1:NDERIV-1
            # Derivative(N)
            dest = fill(T(Inf), order(basis))
            @test bsplines!(dest, basis, x, Derivative(N), leftknot=leftknot) == offset && dest == bspl_nderiv[:,N+1]
            # AllDerivatives(N)
            dest = fill(T(Inf), order(basis), N)
            @test bsplines!(dest, basis, x, AllDerivatives(N), leftknot=leftknot) == offset && dest == bspl_nderiv[:,1:N]
        end
    end

    # b1 = BSplineBasis(1, 0:5)
    for (T, x) = [(Float16, 0), (Float32, 0.5f0), (BigFloat, big(2)), (Rational{Int}, 5//2)]
        test_bsplines!(b1, b1_exact, T, x, true)
    end
    # b2 = BSplineBasis(2, big(0):1:big(5))
    for (T, x) = [(Rational{BigInt}, 8//9), (Rational{Int}, 22//7)]
        test_bsplines!(b2, b2_exact, T, x, true)
    end
    for (T, x) = [(Float32, 0.123), (Float64, 3//2), (BigFloat, 13//8)]
        test_bsplines!(b2, b2_exact, T, x, false)
    end
    # b3 = BSplineBasis(3, 0//1:5//1)
    for (T, x) = [(Rational{Int64}, 8//9), (Rational{Int32}, 1), (Rational{Int128}, 13//8)]
        test_bsplines!(b3, b3_exact, T, x, true)
    end
    for (T, x) = [(Float32, 0), (BigFloat, 2.6f0), (Float64, π)]
        test_bsplines!(b3, b3_exact, T, x, false)
    end
    # b4 = BSplineBasis(4, Int128[0,1,2,3,4,5])
    for (T, x) = [(Rational{Int}, 13//5), (Rational{BigInt}, 5)]
        test_bsplines!(b4, b4_exact, T, x, true)
    end
    for (T, x) = [(Float64, 1.625), (Float32, 9//2)]
        test_bsplines!(b4, b4_exact, T, x, false)
    end
    # b5 = BSplineBasis(3, [-5.0, -3.5, -2.2, 1.0, 2.0, 3.6])
    for (T, x) = [(Float32, -3.0f0), (Float64, Int8(-2)), (BigFloat, π)]
        test_bsplines!(b5, b5_exact, T, x, false)
    end
    # b6 = BSplineBasis(4, Rational{Int64}[-5//1, -7//2, -11//5, 1//1, 2//1, 18//5])
    for (T, x) = [(Rational{Int32}, -5//1), (Rational{Int64}, big(18//5))]
        test_bsplines!(b6, b6_exact, T, x, true)
    end
    for (T, x) = [(Float32, -4), (Float64, -1.0)]
        test_bsplines!(b6, b6_exact, T, x, false)
    end
    # b7 = BSplineBasis(3, Int8[-6,-4,-2,0,0,1,2])
    for (T,x) = [(Rational{Int16}, -6), (Rational{Int64}, big(1//8)),
                 (Rational{BigInt}, Int16(1))]
        test_bsplines!(b7, b7_exact, T, x, true)
    end
    for (T,x) = [(Float32, -9//2), (Float64, Int16(2))]
        test_bsplines!(b7, b7_exact, T, x, false)
    end

    leftknot = intervalindex(b5, 0)
    @test_throws DimensionMismatch bsplines!(zeros(order(b5)+1), b5, 0)
    @test_throws DimensionMismatch bsplines!(zeros(order(b5)-1), b5, 0, leftknot=leftknot)
    @test_throws DimensionMismatch bsplines!(zeros(order(b5), 0), b5, 0, Derivative(1))
    @test_throws DimensionMismatch bsplines!(zeros(order(b5), 1), b5, 0, Derivative(1), leftknot=leftknot)
    @test_throws DimensionMismatch bsplines!(zeros(order(b5), 2), b5, 0, AllDerivatives(3))
    @test_throws DimensionMismatch bsplines!(zeros(order(b5), 3, 1), b5, 0, AllDerivatives(3), leftknot=leftknot)
    dest = zeros(order(b5))
    @test_throws ArgumentError bsplines!(dest, b5, 0, leftknot=leftknot-1)
    @test_throws ArgumentError bsplines!(dest, b5, 0, leftknot=leftknot+1)
    @test_throws ArgumentError bsplines!(dest, b5, 0, leftknot=nothing)
    @test_throws ArgumentError bsplines!(dest, b5, 0, Derivative(0), leftknot=leftknot-1)
    @test_throws ArgumentError bsplines!(dest, b5, 0, Derivative(1), leftknot=leftknot+1)
    @test_throws ArgumentError bsplines!(dest, b5, 0, Derivative(2), leftknot=nothing)
    dest = zeros(order(b5), 3)
    @test_throws ArgumentError bsplines!(dest, b5, 0, AllDerivatives(3), leftknot=leftknot-1)
    @test_throws ArgumentError bsplines!(dest, b5, 0, AllDerivatives(3), leftknot=leftknot+1)
    @test_throws ArgumentError bsplines!(dest, b5, 0, AllDerivatives(3), leftknot=nothing)
    dest1 = fill(NaN, order(b5))
    dest2 = fill(NaN, order(b5))
    @test bsplines!(dest1, b5, 0) == bsplines!(dest2, b5, 0, leftknot=leftknot) && dest1 == dest2
    dest1 = fill(NaN, order(b5))
    dest2 = fill(NaN, order(b5))
    @test bsplines!(dest1, b5, 0, Derivative(0), leftknot=leftknot) == bsplines!(dest2, b5, 0, leftknot=leftknot) && dest1 == dest2
    dest1 = fill(NaN, order(b5))
    dest2 = fill(NaN, order(b5))
    @test bsplines!(dest1, b5, 0, Derivative(1)) == bsplines!(dest2, b5, 0, Derivative(1), leftknot=leftknot) && dest1 == dest2
    dest1 = fill(NaN, order(b5), 3)
    dest2 = fill(NaN, order(b5), 3)
    @test bsplines!(dest1, b5, 0, AllDerivatives(3)) == bsplines!(dest2, b5, 0, AllDerivatives(3), leftknot=leftknot) && dest1 == dest2
end
