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
# * for Integer/Rational: use 64-bit to avoid overflow in rational arithmetic
# * for floating-point: use BigFloat to minimize error
maybebig(x::Integer) = Int64(x)
maybebig(x::Rational) = Rational{Int64}(x)
maybebig(x::AbstractFloat) = big(x)

bsplines!_offsetviewtype(parent, ::Derivative=Derivative(0)) = typeof(OffsetArray(view(parent, 1:0), 1))
bsplines!_offsetviewtype(parent, ::AllDerivatives) = typeof(OffsetArray(view(parent, 1:0, :), 1, -1))

bsplines_offsetviewtype(T, d::Derivative=Derivative(0)) = bsplines!_offsetviewtype(zeros(T,0), d)
bsplines_offsetviewtype(T, d::AllDerivatives) = bsplines!_offsetviewtype(zeros(T,0,0), d)

pparent(x) = parent(parent(x))

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
                @test typeof(@inferred(bsplines(basis, one(xtype), leftknot=leftknot))) === bsplines_offsetviewtype(bspline_type)
                @test typeof(@inferred(bsplines(basis, one(xtype), AllDerivatives(2), leftknot=leftknot))) === bsplines_offsetviewtype(bspline_type, AllDerivatives(2))
                @test typeof(@inferred(bsplines(basis, one(xtype), Derivative(1), leftknot=leftknot))) === bsplines_offsetviewtype(bspline_type, Derivative(1))
            end
            for T = [Float32, Float64, BigFloat, Rational{Int128}]
                derivspace = Matrix{T}(undef, order(basis), order(basis))
                @test typeof(@inferred(bsplines(basis, 1, AllDerivatives(2), leftknot=leftknot, derivspace=derivspace))) === bsplines_offsetviewtype(T, AllDerivatives(2))
                @test typeof(@inferred(bsplines(basis, 1, Derivative(2), leftknot=leftknot, derivspace=derivspace))) === bsplines_offsetviewtype(T, Derivative(2))
            end
        end

        function test_bsplines_empty(basis, x)
            derivspace = Matrix{typeof(x)}(undef, order(basis), order(basis))
            @test size(bsplines(basis, x))                                                             == (0,)
            @test size(bsplines(basis, x, leftknot=nothing))                                           == (0,)
            @test size(bsplines(basis, x, AllDerivatives(2)))                                          == (0, 2)
            @test size(bsplines(basis, x, AllDerivatives(2), leftknot=nothing))                        == (0, 2)
            @test size(bsplines(basis, x, AllDerivatives(2), derivspace=derivspace))                   == (0, 2)
            @test size(bsplines(basis, x, AllDerivatives(2), leftknot=nothing, derivspace=derivspace)) == (0, 2)
            @test size(bsplines(basis, x, Derivative(2)))                                              == (0,)
            @test size(bsplines(basis, x, Derivative(2), leftknot=nothing))                            == (0,)
            @test size(bsplines(basis, x, Derivative(2), derivspace=derivspace))                       == (0,)
            @test size(bsplines(basis, x, Derivative(2), leftknot=nothing, derivspace=derivspace))     == (0,)
        end
        test_bsplines_empty(b1, NaN16)
        test_bsplines_empty(b2, -Inf32)
        test_bsplines_empty(b3, Inf64)
        test_bsplines_empty(b4, 1//0)
        test_bsplines_empty(b5, -1//0)
        test_bsplines_empty(b6, -100)
        test_bsplines_empty(b7, 10.0)
    end

    function test_bsplines(basis, basis_exact, x, isexact)
        NDERIV = 4
        leftknot = intervalindex(basis, x)
        range = leftknot-order(basis)+1:leftknot
        bsp0 = bsplines(basis, x, leftknot=leftknot)
        derivspace = Matrix{eltype(bsp0)}(undef, order(basis), order(basis))
        # AllDerivatives(NDERIVS)
        bspl_nderiv = bsplines(basis, x, AllDerivatives(NDERIV), leftknot=leftknot, derivspace=derivspace)
        bspl_nderiv_exact = bsplines_exact(basis_exact, x, range, AllDerivatives(NDERIV))
        if isexact
            @test bspl_nderiv == bspl_nderiv_exact
        else
            @test bspl_nderiv ≈ₑₗ bspl_nderiv_exact
        end
        @test bsp0 == view(bspl_nderiv,:,0)
        for N = 1:(NDERIV-1)
            # Derivative(N)
            @test bsplines(basis, x, Derivative(N), leftknot=leftknot, derivspace=derivspace) == view(bspl_nderiv,:,N)
            # AllDerivatives(N)
            if N == 1
                @test bsplines(basis, x, AllDerivatives(N), leftknot=leftknot) == OffsetArray(view(pparent(bspl_nderiv),:,1:N), range, 0:N-1)
            else
                @test bsplines(basis, x, AllDerivatives(N), leftknot=leftknot, derivspace=derivspace) == OffsetArray(view(pparent(bspl_nderiv),:,1:N), range, 0:N-1)
            end
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
    derivspace = Matrix{Float64}(undef, order(b5), order(b5))
    @test_throws ArgumentError bsplines(b5, 0, Derivative(0), derivspace=derivspace)
    @test_throws ArgumentError bsplines(b5, 0, AllDerivatives(1), derivspace=derivspace)
    @test_throws DimensionMismatch bsplines(b5, 0, Derivative(1), derivspace=zeros(order(b5)+1, order(b5)))
    @test_throws DimensionMismatch bsplines(b5, 0, AllDerivatives(2), derivspace=zeros(order(b5)))
    @test_throws DimensionMismatch bsplines(b5, 0, AllDerivatives(2), derivspace=zeros(order(b5), order(b5), 1))
    @test bsplines(b5, 0) == bsplines(b5, 0, leftknot=leftknot) == bsplines(b5, 0, Derivative(0), leftknot=leftknot)
    @test bsplines(b5, 0, Derivative(1)) == bsplines(b5, 0, Derivative(1), leftknot=leftknot) ==
                                            bsplines(b5, 0, Derivative(1), derivspace=derivspace) ==
                                            bsplines(b5, 0, Derivative(1), leftknot=leftknot, derivspace=derivspace)
    @test bsplines(b5, 0, AllDerivatives(3)) == bsplines(b5, 0, AllDerivatives(3), leftknot=leftknot) ==
                                                bsplines(b5, 0, AllDerivatives(3), derivspace=derivspace) ==
                                                bsplines(b5, 0, AllDerivatives(3), leftknot=leftknot, derivspace=derivspace)
end

@time @testset "bsplines!" begin
    resultsize_empty(deriv) = (0,)
    resultsize_empty(deriv::AllDerivatives{N}) where N = (0,N)
    function _test_bsplines!_empty(result, dest, deriv=Derivative(0))
        @test typeof(result) === bsplines!_offsetviewtype(dest, deriv)
        @test size(result) === resultsize_empty(deriv)
        @test pparent(result) === dest
    end
    function test_bsplines!_empty(basis, x)
        dest     = Array{Float64}(undef, order(basis))
        dest_all = Array{Float32}(undef, order(basis), 2)
        derivspace     = Array{Float64}(undef, order(basis), order(basis))
        derivspace_all = Array{Float32}(undef, order(basis), order(basis))
        _test_bsplines!_empty(@inferred(bsplines!(dest, basis, x)),                   dest)
        _test_bsplines!_empty(@inferred(bsplines!(dest, basis, x, leftknot=nothing)), dest)
        for (dst, drv, drvspace) = [(dest, Derivative(2), derivspace), (dest_all, AllDerivatives(2), derivspace_all)]
            _test_bsplines!_empty(@inferred(bsplines!(dst, basis, x, drv)),                                        dst, drv)
            _test_bsplines!_empty(@inferred(bsplines!(dst, basis, x, drv, leftknot=nothing)),                      dst, drv)
            _test_bsplines!_empty(@inferred(bsplines!(dst, basis, x, drv, derivspace=drvspace)),                   dst, drv)
            _test_bsplines!_empty(@inferred(bsplines!(dst, basis, x, drv, leftknot=nothing, derivspace=drvspace)), dst, drv)
        end
    end
    test_bsplines!_empty(b1, NaN16)
    test_bsplines!_empty(b2, -Inf32)
    test_bsplines!_empty(b3, Inf64)
    test_bsplines!_empty(b4, -100)
    test_bsplines!_empty(b5, 10.0)

    function test_bsplines!(basis, basis_exact, T, x, isexact)
        NDERIV = 6
        leftknot = intervalindex(basis, x)
        offset = leftknot-order(basis)
        range = offset+1:leftknot
        derivspace = Matrix{T}(undef, order(basis), order(basis))

        # AllDerivatives(NDERIV)
        dest0 = fill(T(Inf), order(basis), NDERIV)
        bspl_nderiv_exact = bsplines_exact(basis_exact, x, range, AllDerivatives(NDERIV))
        bspl_nderiv = bsplines!(dest0, basis, x, AllDerivatives(NDERIV), leftknot=leftknot, derivspace=derivspace)
        if isexact
            @test bspl_nderiv == bspl_nderiv_exact
        else
            @test bspl_nderiv ≈ₑₗ bspl_nderiv_exact
        end
        @test pparent(bspl_nderiv) === dest0

        # No derivative
        dest1 = fill(T(Inf), order(basis))
        bspl1 = bsplines!(dest1, basis, x, leftknot=leftknot)
        @test bspl1 == view(bspl_nderiv,:,0)
        @test pparent(bspl1) === dest1

        # AllDerivatives(N) and Derivative(N) for N = 1:NDERIV-1
        for N = 1:NDERIV-1
            # Derivative(N)
            fill!(dest1, T(Inf))
            bspl1 = bsplines!(dest1, basis, x, Derivative(N), leftknot=leftknot, derivspace=derivspace)
            @test bspl1 == view(bspl_nderiv,:,N)
            @test pparent(bspl1) === dest1
            # AllDerivatives(N)
            dest2 = fill(T(Inf), order(basis), N)
            if N == 1
                bspl2 = bsplines!(dest2, basis, x, AllDerivatives(N), leftknot=leftknot)
            else
                bspl2 = bsplines!(dest2, basis, x, AllDerivatives(N), leftknot=leftknot, derivspace=derivspace)
            end
            @test bspl2 == OffsetArray(view(pparent(bspl_nderiv),:,1:N), offset, -1)
            @test pparent(bspl2) === dest2
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
    derivspace = zeros(order(b5), order(b5))
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
    @test_throws ArgumentError bsplines!(dest, b5, 0, Derivative(0), derivspace=derivspace)
    @test_throws DimensionMismatch bsplines!(dest, b5, 0, Derivative(1), derivspace=zeros(order(b5)+1, order(b5)))
    @test_throws DimensionMismatch bsplines!(dest, b5, 0, Derivative(2), derivspace=zeros(order(b5), order(b5), 1))
    dest = zeros(order(b5), 1)
    @test_throws ArgumentError bsplines!(dest, b5, 0, AllDerivatives(1), derivspace=derivspace)
    dest = zeros(order(b5), 3)
    @test_throws ArgumentError bsplines!(dest, b5, 0, AllDerivatives(3), leftknot=leftknot-1)
    @test_throws ArgumentError bsplines!(dest, b5, 0, AllDerivatives(3), leftknot=leftknot+1)
    @test_throws ArgumentError bsplines!(dest, b5, 0, AllDerivatives(3), leftknot=nothing)
    @test_throws DimensionMismatch bsplines!(dest, b5, 0, Derivative(1), derivspace=zeros(order(b5), order(b5)-1))
    @test_throws DimensionMismatch bsplines!(dest, b5, 0, Derivative(2), derivspace=zeros(order(b5)))
    dest1 = fill(NaN, order(b5))
    dest2 = fill(NaN, order(b5))
    bsp1 = bsplines!(dest1, b5, 0)
    bsp2 = bsplines!(dest2, b5, 0, leftknot=leftknot)
    @test bsp1 == bsp2
    @test pparent(bsp1) === dest1
    @test pparent(bsp2) === dest2
    fill!(dest1, NaN)
    fill!(dest2, NaN)
    dest3 = fill(NaN, order(b5))
    dest4 = fill(NaN, order(b5))
    bsp1 = bsplines!(dest1, b5, 0, Derivative(1))
    bsp2 = bsplines!(dest2, b5, 0, Derivative(1), leftknot=leftknot)
    bsp3 = bsplines!(dest3, b5, 0, Derivative(1), derivspace=derivspace)
    bsp4 = bsplines!(dest4, b5, 0, Derivative(1), leftknot=leftknot, derivspace=derivspace)
    @test bsp1 == bsp2 == bsp3 == bsp4
    @test pparent(bsp1) === dest1
    @test pparent(bsp2) === dest2
    @test pparent(bsp3) === dest3
    @test pparent(bsp4) === dest4
    dest1 = fill(NaN, order(b5), 3)
    dest2 = fill(NaN, order(b5), 3)
    dest3 = fill(NaN, order(b5), 3)
    dest4 = fill(NaN, order(b5), 3)
    bsp1 = bsplines!(dest1, b5, 0, AllDerivatives(3))
    bsp2 = bsplines!(dest2, b5, 0, AllDerivatives(3), leftknot=leftknot)
    bsp3 = bsplines!(dest3, b5, 0, AllDerivatives(3), derivspace=derivspace)
    bsp4 = bsplines!(dest4, b5, 0, AllDerivatives(3), leftknot=leftknot, derivspace=derivspace)
    @test bsp1 == bsp2 == bsp3 == bsp4
    @test pparent(bsp1) === dest1
    @test pparent(bsp2) === dest2
    @test pparent(bsp3) === dest3
    @test pparent(bsp4) === dest4
end
