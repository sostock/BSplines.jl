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

maybebig(x::AbstractFloat) = big(x)
maybebig(x) = x

@time @testset "bsplines" begin
    @testset "Return types" begin
        function test_bsplines_nothing(basis, x)
            @test bsplines(basis, x)                                      === nothing
            @test bsplines(basis, x, leftknot=nothing)                    === nothing
            @test bsplines(basis, x, AllDerivatives(2))                   === nothing
            @test bsplines(basis, x, AllDerivatives(2), leftknot=nothing) === nothing
            @test bsplines(basis, x, Derivative(2))                       === nothing
            @test bsplines(basis, x, Derivative(2), leftknot=nothing)     === nothing
        end

        # b1 = BSplineBasis(1, 0:5)
        b1_types = [(Int,Float64), (Float64,Float64), (Float32,Float32),
                    (Rational{Int},Rational{Int}), (BigInt,BigFloat),
                    (BigFloat,BigFloat), (Rational{BigInt},Rational{BigInt}),
                    (Rational{Int16}, Rational{Int})]
        # b2 = BSplineBasis(2, big(0):1:big(5))
        b2_types = [(Int,BigFloat), (Float64,BigFloat), (Float32,BigFloat),
                    (Rational{Int},Rational{BigInt}), (BigInt,BigFloat),
                    (BigFloat,BigFloat), (Rational{BigInt},Rational{BigInt})]
        # b3 = BSplineBasis(3, 0//1:5//1)
        b3_types = [(Int,Rational{Int}), (Float64,Float64), (Float32,Float32),
                    (Rational{Int},Rational{Int}), (BigInt,Rational{BigInt}),
                    (BigFloat,BigFloat), (Rational{BigInt},Rational{BigInt})]
        # b4 = BSplineBasis(4, Int128[0,1,2,3,4,5])
        b4_types = [(Int,Float64), (Float64,Float64), (Float32,Float32),
                    (Rational{Int},Rational{Int128}), (BigInt,BigFloat),
                    (BigFloat,BigFloat), (Rational{BigInt},Rational{BigInt})]
        # b5 = BSplineBasis(3, [-5.0, -3.5, -2.2, 1.0, 2.0, 3.6])
        b5_types = [(Int,Float64), (Float64,Float64), (Float32,Float64),
                    (Rational{Int},Float64), (BigInt,BigFloat),
                    (BigFloat,BigFloat), (Rational{BigInt},BigFloat)]
        # b6 = BSplineBasis(4, Rational{BigInt}[-5//1, -7//2, -11//5, 1//1, 2//1, 18//5])
        b6_types = [(Int,Rational{BigInt}), (Float64,BigFloat), (Float32,BigFloat),
                    (Rational{Int},Rational{BigInt}), (BigInt,Rational{BigInt}),
                    (BigFloat,BigFloat), (Rational{BigInt},Rational{BigInt})]
        # b7 = BSplineBasis(3, Int8[-6,-4,-2,0,0,1,2])
        b7_types = [(Int,Float64), (Float64,Float64), (Float32,Float32),
                    (Rational{Int},Rational{Int}), (BigInt,BigFloat),
                    (BigFloat,BigFloat), (Rational{BigInt},Rational{BigInt}),
                    (Rational{Int16},Rational{Int16})]
        for (basis, types) = [(b1, b1_types), (b2, b2_types), (b3, b3_types),
                              (b4, b4_types), (b5, b5_types), (b6, b6_types),
                              (b7, b7_types)]
            leftknot = intervalindex(basis, 1)
            for (xtype, bspline_type) in types
                @test eltype(bsplines(basis, one(xtype)))                    == bspline_type
                @test eltype(bsplines(basis, one(xtype), leftknot=leftknot)) == bspline_type
                @test eltype(bsplines(basis, one(xtype), AllDerivatives(2)))                    == bspline_type
                @test eltype(bsplines(basis, one(xtype), AllDerivatives(2), leftknot=leftknot)) == bspline_type
                @test eltype(bsplines(basis, one(xtype), Derivative(0)))                    == bspline_type
                @test eltype(bsplines(basis, one(xtype), Derivative(0), leftknot=leftknot)) == bspline_type
                test_bsplines_nothing(basis, xtype(-10))
                test_bsplines_nothing(basis, xtype(10))
            end
            min, max = extrema(breakpoints(basis))
            for T = (Float16, Float32, Float64, BigFloat)
                for x = (T(-Inf), T(Inf), T(NaN), prevfloat(T(min)), nextfloat(T(max)))
                    test_bsplines_nothing(basis, x)
                end
            end
        end
    end

    function test_bsplines(basis, basis_exact, x, isexact)
        leftknot = intervalindex(basis, x)
        range = leftknot-order(basis)+1:leftknot
        exact = bsplines_exact(basis_exact, x, range, AllDerivatives(6))
        # AllDerivatives(6)
        bspl6 = bsplines(basis, x, AllDerivatives(6))
        if isexact
            @test bspl6 == exact
        else
            @test bspl6 ≈ₑₗ exact
        end
        @test bsplines(basis, x, AllDerivatives(6), leftknot=leftknot) == bspl6
        @test_throws ArgumentError bsplines(basis, x, AllDerivatives(6), leftknot=leftknot-1)
        @test_throws ArgumentError bsplines(basis, x, AllDerivatives(6), leftknot=leftknot+1)
        @test_throws ArgumentError bsplines(basis, x, AllDerivatives(6), leftknot=nothing)
        # No Derivative
        ex = bspl6[:,0]
        @test bsplines(basis, x)                    == ex
        @test bsplines(basis, x, leftknot=leftknot) == ex
        @test_throws ArgumentError bsplines(basis, x, leftknot=leftknot-1)
        @test_throws ArgumentError bsplines(basis, x, leftknot=leftknot+1)
        @test_throws ArgumentError bsplines(basis, x, leftknot=nothing)
        # Derivative(0)
        @test bsplines(basis, x, Derivative(0))                    == ex
        @test bsplines(basis, x, Derivative(0), leftknot=leftknot) == ex
        @test_throws ArgumentError bsplines(basis, x, Derivative(0), leftknot=leftknot-1)
        @test_throws ArgumentError bsplines(basis, x, Derivative(0), leftknot=leftknot+1)
        @test_throws ArgumentError bsplines(basis, x, Derivative(0), leftknot=nothing)
        # AllDerivatives(N) and Derivative(N) for N = 1:5
        for N = 1:5
            # Derivative(N)
            ex = bspl6[:, N]
            @test bsplines(basis, x, Derivative(N))                    == ex
            @test bsplines(basis, x, Derivative(N), leftknot=leftknot) == ex
            @test_throws ArgumentError bsplines(basis, x, Derivative(N), leftknot=leftknot-1)
            @test_throws ArgumentError bsplines(basis, x, Derivative(N), leftknot=leftknot+1)
            @test_throws ArgumentError bsplines(basis, x, Derivative(N), leftknot=nothing)
            # AllDerivatives(N)
            ex_all = OffsetArray(bspl6[range, 0:N-1], range, 0:N-1)
            @test bsplines(basis, x, AllDerivatives(N))                    == ex_all
            @test bsplines(basis, x, AllDerivatives(N), leftknot=leftknot) == ex_all
            @test_throws ArgumentError bsplines(basis, x, AllDerivatives(N), leftknot=leftknot-1)
            @test_throws ArgumentError bsplines(basis, x, AllDerivatives(N), leftknot=leftknot+1)
            @test_throws ArgumentError bsplines(basis, x, AllDerivatives(N), leftknot=nothing)
        end
    end

    # b1 = BSplineBasis(1, 0:5)
    for x = Real[0, 1//2, 1.0f0, 4.1, Int8(5)]
        test_bsplines(b1, b1_exact, x, true)
    end
    # b2 = BSplineBasis(2, big(0):1:big(5))
    for x = Real[Float16(0.0), 1.625, 2.6f0, big(5)]
        test_bsplines(b2, b2_exact, x, false)
    end
    for x = Real[8//9, 22//7, big(9//2), 5//1]
        test_bsplines(b2, b2_exact, x, true)
    end
    # b3 = BSplineBasis(3, 0//1:5//1)
    for x = Real[0.123, Float16(1.625), π, 5.0]
        test_bsplines(b3, b3_exact, x, false)
    end
    for x = Real[0//1, big(1), 13//5, big(9//2)]
        test_bsplines(b3, b3_exact, x, true)
    end
    # b4 = BSplineBasis(4, Int128[0,1,2,3,4,5])
    for x = Real[0, 2.6f0, big(4.5), Int8(5)]
        test_bsplines(b4, b4_exact, x, false)
    end
    for x = Real[13//8, 13//5, big(9//2), 5//1]
        test_bsplines(b4, b4_exact, x, true)
    end
    # b5 = BSplineBasis(3, [-5.0, -3.5, -2.2, 1.0, 2.0, 3.6])
    for x = Real[-5, -4//1, -3.0f0, π, 3.6]
        test_bsplines(b5, b5_exact, x, false)
    end
    # b6 = BSplineBasis(4, Rational{BigInt}[-5//1, -7//2, -11//5, 1//1, 2//1, 18//5])
    for x = Real[-3.0f0, Float16(-2.0), big(1.5)]
        test_bsplines(b6, b6_exact, x, false)
    end
    for x = Real[-5, Int8(-3), big(3//2), 18//5]
        test_bsplines(b6, b6_exact, x, true)
    end
    # b7 = BSplineBasis(3, Int8[-6,-4,-2,0,0,1,2])
    for x = Real[-6.0, -2.6f0, big(0.375), Int8(2)]
        test_bsplines(b7, b7_exact, x, false)
    end
    for x = Real[-9//2, Int16(0)//Int16(1), big(5//7), 5//4]
        test_bsplines(b7, b7_exact, x, true)
    end
end

@time @testset "bsplines!" begin
    function test_bsplines!_nothing(basis, x, dest, dest_all)
        @test bsplines!(dest, basis, x) === nothing
        @test bsplines!(dest, basis, x, leftknot=nothing) === nothing
        @test bsplines!(dest_all, basis, x, AllDerivatives(2)) === nothing
        @test bsplines!(dest_all, basis, x, AllDerivatives(2), leftknot=nothing) === nothing
        @test bsplines!(dest, basis, x, Derivative(2)) === nothing
        @test bsplines!(dest, basis, x, Derivative(2), leftknot=nothing) === nothing
    end

    for basis = (b1, b2, b3, b4, b5, b6, b7)
        dest     = Array{Float32}(undef, order(basis))
        dest_all = Array{Float64}(undef, order(basis), 2)
        for T = (Int, BigInt, Rational{Int}, Rational{BigInt})
            test_bsplines!_nothing(basis, T(-10), dest, dest_all)
            test_bsplines!_nothing(basis, T(10),  dest, dest_all)
        end
        min, max = extrema(breakpoints(basis))
        for T = (Float16, Float32, Float64, BigFloat)
            for x = (T(-Inf), T(Inf), T(NaN), prevfloat(T(min)), nextfloat(T(max)))
                test_bsplines!_nothing(basis, x, dest, dest_all)
            end
        end
    end

    function test_bsplines!(basis, basis_exact, T, x, isexact)
        leftknot = intervalindex(basis, x)
        offset = leftknot-order(basis)
        range = offset+1:leftknot
        exact = parent(bsplines_exact(basis_exact, x, range, AllDerivatives(6)))

        # AllDerivatives(6)
        bspl6 = Array{T}(undef, order(basis), 6)
        dest  = Array{T}(undef, order(basis), 6)
        if isexact
            @test begin
                fill!(bspl6, T(Inf))
                off = bsplines!(bspl6, basis, x, AllDerivatives(6))
                off == offset && bspl6 == exact
            end
        else
            @test begin
                fill!(bspl6, T(Inf))
                off = bsplines!(bspl6, basis, x, AllDerivatives(6))
                off == offset && bspl6 ≈ₑₗ exact
            end
        end
        @test begin
            fill!(dest, T(Inf))
            off = bsplines!(dest, basis, x, AllDerivatives(6), leftknot=leftknot)
            off == offset && dest == bspl6
        end
        @test_throws ArgumentError bsplines!(dest, basis, x, AllDerivatives(6), leftknot=leftknot-1)
        @test_throws ArgumentError bsplines!(dest, basis, x, AllDerivatives(6), leftknot=leftknot+1)
        @test_throws ArgumentError bsplines!(dest, basis, x, AllDerivatives(6), leftknot=nothing)
        for dest = [Array{T}(undef, order(basis)-1, 6), Array{T}(undef, order(basis)+1, 6),
                    Array{T}(undef, order(basis), 5), Array{T}(undef, order(basis), 7),
                    Array{T}(undef, order(basis)), Array{T}(undef, order(basis), 6, 1)]
            @test_throws DimensionMismatch bsplines!(dest, basis, x, AllDerivatives(6))
            @test_throws DimensionMismatch bsplines!(dest, basis, x, AllDerivatives(6), leftknot=leftknot)
        end

        # No derivative/Derivative(0)
        dest = Array{T}(undef, order(basis))
        ex = bspl6[:,1]
        @test begin
            fill!(dest, T(Inf))
            off = bsplines!(dest, basis, x)
            off == offset && dest == ex
        end
        @test begin
            fill!(dest, T(Inf))
            off = bsplines!(dest, basis, x, leftknot=leftknot)
            off == offset && dest == ex
        end
        @test begin
            fill!(dest, T(Inf))
            off = bsplines!(dest, basis, x, Derivative(0))
            off == offset && dest == ex
        end
        @test begin
            fill!(dest, T(Inf))
            off = bsplines!(dest, basis, x, Derivative(0), leftknot=leftknot)
            off == offset && dest == ex
        end
        @test_throws ArgumentError bsplines!(dest, basis, x, leftknot=leftknot-1)
        @test_throws ArgumentError bsplines!(dest, basis, x, leftknot=leftknot+1)
        @test_throws ArgumentError bsplines!(dest, basis, x, leftknot=nothing)
        @test_throws ArgumentError bsplines!(dest, basis, x, Derivative(0), leftknot=leftknot-1)
        @test_throws ArgumentError bsplines!(dest, basis, x, Derivative(0), leftknot=leftknot+1)
        @test_throws ArgumentError bsplines!(dest, basis, x, Derivative(0), leftknot=nothing)
        for dest = [Array{T}(undef, order(basis)+1), Array{T}(undef, order(basis)-1),
                    Array{T}(undef, order(basis), 1), Array{T}(undef, order(basis), 2)]
            @test_throws DimensionMismatch bsplines!(dest, basis, x)
            @test_throws DimensionMismatch bsplines!(dest, basis, x, leftknot=leftknot)
            @test_throws DimensionMismatch bsplines!(dest, basis, x, Derivative(0))
            @test_throws DimensionMismatch bsplines!(dest, basis, x, Derivative(0), leftknot=leftknot)
        end

        # AllDerivatives(N) and Derivative(N) for N = 1:5
        for N = 1:5
            # Derivative(N)
            dest = Array{T}(undef, order(basis))
            ex = bspl6[:,N+1]
            @test begin
                fill!(dest, T(Inf))
                off = bsplines!(dest, basis, x, Derivative(N))
                off == offset && dest == ex
            end
            @test begin
                fill!(dest, T(Inf))
                off = bsplines!(dest, basis, x, Derivative(N), leftknot=leftknot)
                off == offset && dest == ex
            end
            @test_throws ArgumentError bsplines!(dest, basis, x, Derivative(N), leftknot=leftknot-1)
            @test_throws ArgumentError bsplines!(dest, basis, x, Derivative(N), leftknot=leftknot+1)
            @test_throws ArgumentError bsplines!(dest, basis, x, Derivative(N), leftknot=nothing)
            for dest = [Array{T}(undef, order(basis)+1), Array{T}(undef, order(basis)-1),
                        Array{T}(undef, order(basis), 1)]
                @test_throws DimensionMismatch bsplines!(dest, basis, x, Derivative(N))
                @test_throws DimensionMismatch bsplines!(dest, basis, x, Derivative(N), leftknot=leftknot)
            end
            # AllDerivatives(N)
            dest = Array{T}(undef, order(basis), N)
            ex = bspl6[:,1:N]
            @test begin
                fill!(dest, T(Inf))
                off = bsplines!(dest, basis, x, AllDerivatives(N))
                off == offset && dest == ex
            end
            @test begin
                fill!(dest, T(Inf))
                off = bsplines!(dest, basis, x, AllDerivatives(N), leftknot=leftknot)
                off == offset && dest == ex
            end
            @test_throws ArgumentError bsplines!(dest, basis, x, AllDerivatives(N), leftknot=leftknot-1)
            @test_throws ArgumentError bsplines!(dest, basis, x, AllDerivatives(N), leftknot=leftknot+1)
            @test_throws ArgumentError bsplines!(dest, basis, x, AllDerivatives(N), leftknot=nothing)
            for dest = [Array{T}(undef, order(basis)-1, N), Array{T}(undef, order(basis)+1, N),
                        Array{T}(undef, order(basis), N-1), Array{T}(undef, order(basis), N+1),
                        Array{T}(undef, order(basis)), Array{T}(undef, order(basis), N, 1)]
                @test_throws DimensionMismatch bsplines!(dest, basis, x, AllDerivatives(N))
                @test_throws DimensionMismatch bsplines!(dest, basis, x, AllDerivatives(N), leftknot=leftknot)
            end
        end
    end

    # b1 = BSplineBasis(1, 0:5)
    for (T, x) = [(Float16, 0), (Float32, 0.5f0), (BigFloat, big(2)), (Rational{Int}, 5//2)]
        test_bsplines!(b1, b1_exact, T, x, true)
    end
    # b2 = BSplineBasis(2, big(0):1:big(5))
    for (T, x) = [(Rational{Int}, 0), (Rational{BigInt}, 8//9), (Rational{Int}, Int8(1)), (Rational{Int}, big(22//7))]
        test_bsplines!(b2, b2_exact, T, x, true)
    end
    for (T, x) = [(Float32, 0.123), (Float64, 3//2), (BigFloat, 13//8), (Float16, big(5))]
        test_bsplines!(b2, b2_exact, T, x, false)
    end
    # b3 = BSplineBasis(3, 0//1:5//1)
    for (T, x) = [(Rational{Int}, 8//9), (Rational{BigInt}, 1), (Rational{BigInt}, 13//8), (Rational{Int}, big(5))]
        test_bsplines!(b3, b3_exact, T, x, true)
    end
    for (T, x) = [(Float32, 0), (BigFloat, big(0.123)), (BigFloat, 2.6f0), (Float64, π)]
        test_bsplines!(b3, b3_exact, T, x, false)
    end
    # b4 = BSplineBasis(4, Int128[0,1,2,3,4,5])
    for (T, x) = [(Rational{Int}, 8//9), (Rational{BigInt}, 13//8), (Rational{Int}, big(13//5)),
                  (Rational{BigInt}, 1), (Rational{Int}, big(5))]
        test_bsplines!(b4, b4_exact, T, x, true)
    end
    for (T, x) = [(Float64, 1.625), (Float64, big(9//2)), (Rational{Int}, big(5))]
        test_bsplines!(b4, b4_exact, T, x, false)
    end
    # b5 = BSplineBasis(3, [-5.0, -3.5, -2.2, 1.0, 2.0, 3.6])
    for (T, x) = [(Float32, -3.0f0), (Float64, Int8(-2)), (Float64, π), (Float32, 3.6)]
        test_bsplines!(b5, b5_exact, T, x, false)
    end
    # b6 = BSplineBasis(4, Rational{BigInt}[-5//1, -7//2, -11//5, 1//1, 2//1, 18//5])
    for (T, x) = [(Rational{Int}, -5//1), (Rational{Int}, Int8(-1)), (Rational{BigInt}, Int128(2)),
                  (Rational{BigInt}, big(18//5))]
        test_bsplines!(b6, b6_exact, T, x, true)
    end
    for (T, x) = [(Float32, -4), (BigFloat, -2.0), (Float64, 3//2), (BigFloat, π)]
        test_bsplines!(b6, b6_exact, T, x, false)
    end
    # b7 = BSplineBasis(3, Int8[-6,-4,-2,0,0,1,2])
    for (T,x) = [(Rational{Int16}, -6), (Rational{Int}, -3//1), (Rational{BigInt}, big(1//8)), (Rational{BigInt}, Int16(1))]
        test_bsplines!(b7, b7_exact, T, x, true)
    end
    for (T,x) = [(Float32, -9//2), (Float64, 0.0), (BigFloat, -big(π)), (Float64, Int16(2))]
        test_bsplines!(b7, b7_exact, T, x, false)
    end
end
