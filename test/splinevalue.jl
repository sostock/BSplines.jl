@time @testset "splinevalue" begin
    @testset "Return types" begin
        function test_splinevalue_returntype(spline, xtype, bspline_type)
            leftknot = intervalindex(basis(spline), 1)
            x = one(xtype)
            @test @inferred(splinevalue(spline, x)) isa bspline_type
            @test @inferred(splinevalue(spline, x, leftknot=leftknot)) isa bspline_type
            @test @inferred(splinevalue(spline, x, Derivative(2))) isa bspline_type
            @test @inferred(splinevalue(spline, x, Derivative(2), leftknot=leftknot)) isa bspline_type
            @test @inferred(spline(x)) isa bspline_type
            @test @inferred(spline(x, leftknot=leftknot)) isa bspline_type
            @test @inferred(spline(x, Derivative(2))) isa bspline_type
            @test @inferred(spline(x, Derivative(2), leftknot=leftknot)) isa bspline_type
            @test_throws ArgumentError splinevalue(spline, x, Derivative(2), leftknot=leftknot-1)
            @test_throws ArgumentError splinevalue(spline, x, Derivative(2), leftknot=leftknot+1)
            @test_throws ArgumentError splinevalue(spline, x, Derivative(2), leftknot=nothing)
            @test_throws ArgumentError spline(x, Derivative(2), leftknot=leftknot-1)
            @test_throws ArgumentError spline(x, Derivative(2), leftknot=leftknot+1)
            @test_throws ArgumentError spline(x, Derivative(2), leftknot=nothing)
            @test @inferred(splinevalue(spline, xtype(100), leftknot=nothing)) isa bspline_type
            @test @inferred(splinevalue(spline, xtype(100), Derivative(2), leftknot=nothing)) isa bspline_type
            @test @inferred(spline(xtype(100), leftknot=nothing)) isa bspline_type
            @test @inferred(spline(xtype(100), Derivative(2), leftknot=nothing)) isa bspline_type
            @test_throws ArgumentError splinevalue(spline, xtype(100), leftknot=leftknot)
            @test_throws ArgumentError splinevalue(spline, xtype(100), Derivative(2), leftknot=leftknot)
            @test_throws ArgumentError spline(xtype(100), leftknot=leftknot)
            @test_throws ArgumentError spline(xtype(100), Derivative(2), leftknot=leftknot)
        end

        # b1 = BSplineBasis(1, 0:5)
        b1_bspline_types = [(Int,Float64), (Float64,Float64), (Float32,Float32),
                            (Rational{Int},Rational{Int}), (BigInt,BigFloat),
                            (BigFloat,BigFloat), (Rational{BigInt},Rational{BigInt})]
        # b2 = BSplineBasis(2, big(0):1:big(5))
        b2_bspline_types = [(Int,BigFloat), (Float64,BigFloat), (Float32,BigFloat),
                            (Rational{Int},Rational{BigInt}), (BigInt,BigFloat),
                            (BigFloat,BigFloat), (Rational{BigInt},Rational{BigInt})]
        # b3 = BSplineBasis(3, 0//1:5//1)
        b3_bspline_types = [(Int,Rational{Int}), (Float64,Float64), (Float32,Float32),
                            (Rational{Int},Rational{Int}), (BigInt,Rational{BigInt}),
                            (BigFloat,BigFloat), (Rational{BigInt},Rational{BigInt})]
        # b4 = BSplineBasis(4, Int128[0,1,2,3,4,5])
        b4_bspline_types = [(Int,Float64), (Float64,Float64), (Float32,Float32),
                            (Rational{Int},Rational{Int128}), (BigInt,BigFloat),
                            (BigFloat,BigFloat), (Rational{BigInt},Rational{BigInt})]
        # b5 = BSplineBasis(3, [-5.0, -3.5, -2.2, 1.0, 2.0, 3.6])
        b5_bspline_types = [(Int,Float64), (Float64,Float64), (Float32,Float64),
                            (Rational{Int},Float64), (BigInt,BigFloat),
                            (BigFloat,BigFloat), (Rational{BigInt},BigFloat)]
        # b6 = BSplineBasis(4, Rational{BigInt}[-5//1, -7//2, -11//5, 1//1, 2//1, 18//5])
        b6_bspline_types = [(Int,Rational{BigInt}), (Float64,BigFloat), (Float32,BigFloat),
                            (Rational{Int},Rational{BigInt}), (BigInt,Rational{BigInt}),
                            (BigFloat,BigFloat), (Rational{BigInt},Rational{BigInt})]
        # b7 = BSplineBasis(3, Int8[-6,-4,-2,0,0,1,2])
        b7_bspline_types = [(Int,Float64), (Float64,Float64), (Float32,Float32),
                            (Rational{Int},Rational{Int}), (BigInt,BigFloat),
                            (BigFloat,BigFloat), (Rational{BigInt},Rational{BigInt})]
        for (basis, types) = [(b1, b1_bspline_types), (b2, b2_bspline_types), (b3, b3_bspline_types),
                              (b4, b4_bspline_types), (b5, b5_bspline_types), (b6, b6_bspline_types),
                              (b7, b7_bspline_types)]
            bspline = basis[1]
            for (xtype, bspline_type) in types
                test_splinevalue_returntype(bspline, xtype, bspline_type)
            end
        end

        # b1 = BSplineBasis(1, 0:5)
        b1_coeffs_and_types = [(Float32.(1:5), Int, Float32), (1:5, Int, Float64), (1:5, BigInt, BigFloat),
                               (Float64.(1:5), Int, Float64), (1:5, Float32, Float32), (1:5, Rational{Int}, Rational{Int}),
                               (1//2:9//2, BigInt, Rational{BigInt}), (1//2:9//2, Float64, Float64),
                               (Float16.(1:5), BigFloat, BigFloat), (big(1//2):big(9//2), Float64, BigFloat)]
        # b2 = BSplineBasis(2, big(0):1:big(5))
        b2_coeffs_and_types = [(Float64.(1:6), Int, BigFloat), (Float32.(1:6), Float32, BigFloat),
                               (1//1:6//1, Float64, BigFloat), (1//1:6//1, Int8, Rational{BigInt}), (1:6, Int, BigFloat),
                               (1:6, Rational{Int}, Rational{BigInt}), (1:6, Int, BigFloat)]
        # b3 = BSplineBasis(3, 0//1:5//1)
        b3_coeffs_and_types = [(1:7, Int, Rational{Int}), (1:7, Float32, Float32), (1:7, Float64, Float64),
                               (Float32.(1:7), Rational{Int}, Float32), (Float64.(1:7), BigFloat, BigFloat),
                               (1:7, BigInt, Rational{BigInt}), (Float32.(1:7), Float64, Float64),
                               (Float64.(1:7), Float32, Float64), (BigFloat.(1:7), Rational{Int}, BigFloat)]
        # b4 = BSplineBasis(4, Int128[0,1,2,3,4,5])
        b4_coeffs_and_types = [(Float32.(1:8), Int, Float32), (Float64.(1:8), Rational{Int}, Float64),
                               (BigFloat.(1:8), Int, BigFloat), (BigFloat.(1:8), Float64, BigFloat),
                               (1:8, Rational{Int}, Rational{Int128}), (1:8, BigInt, BigFloat)]
        # b5 = BSplineBasis(3, [-5.0, -3.5, -2.2, 1.0, 2.0, 3.6])
        b5_coeffs_and_types = [(1:7, Int, Float64), (Float32.(1:7), Float32, Float64), (1//1:7//1, BigInt, BigFloat),
                               (Float64.(1:7), BigFloat, BigFloat), (big(1):big(7), Float64, BigFloat)]
        # b6 = BSplineBasis(4, Rational{BigInt}[-5//1, -7//2, -11//5, 1//1, 2//1, 18//5])
        b6_coeffs_and_types = [(1:8, Int, Rational{BigInt}), (Float32.(1:8), Float32, BigFloat), (1:8, Float64, BigFloat),
                               (1//1:8//1, Int, Rational{BigInt}), (Float64.(1:8), Float64, BigFloat)]
        # b7 = BSplineBasis(3, Int8[-6,-4,-2,0,0,1,2])
        b7_coeffs_and_types = [(Float32.(1:8), Int, Float32), (Float64.(1:8), Rational{Int}, Float64),
                               (BigFloat.(1:8), Int, BigFloat), (BigFloat.(1:8), Float64, BigFloat),
                               (1:8, Rational{Int}, Rational{Int}), (1:8, BigInt, BigFloat)]
        for (basis, coeffs_and_types) = [(b1, b1_coeffs_and_types), (b2, b2_coeffs_and_types), (b3, b3_coeffs_and_types),
                                         (b4, b4_coeffs_and_types), (b5, b5_coeffs_and_types), (b6, b6_coeffs_and_types),
                                         (b7, b7_coeffs_and_types)]
            for (coeffs, xtype, bspline_type) in coeffs_and_types
                spline = Spline(basis, coeffs)
                test_splinevalue_returntype(spline, xtype, bspline_type)
            end
        end
    end

    splinevalue_exact_bspline(splines, i, x, nderiv=0) =
        nderiv ≥ size(splines,2) ? 0 : (splines[i,nderiv+1])(maybebig(x))

    maybebig(x) = x
    maybebig(x::AbstractFloat) = big(x)

    function test_splinevalue_zero(spline::Spline, x, isexact)
        spl = splinevalue(spline, x)
        isexact ? (@test spl == 0) : (@test spl ≈ₑₗ 0)
        @test splinevalue(spline, x, leftknot=nothing) == spl
        @test_throws ArgumentError splinevalue(spline, x, leftknot=1)
        for N = 0:5
            spl = splinevalue(spline, x, Derivative(N))
            isexact ? (@test spl == 0) : (@test spl ≈ₑₗ 0)
            @test splinevalue(spline, x, Derivative(N), leftknot=nothing) == spl
            @test_throws ArgumentError splinevalue(spline, x, Derivative(N), leftknot=1)
        end
    end

    function test_splinevalue_bsplines_nonzero(bspline::BSpline, basis_exact, i, x, isexact)
        leftknot = intervalindex(basis(bspline), x)
        leftknot === nothing && error("this function is only suitable for x for which the B-spline is non-zero")
        spl = splinevalue(bspline, x)
        if isexact
            @test spl == splinevalue_exact_bspline(basis_exact, i, x)
        else
            @test spl ≈ₑₗ splinevalue_exact_bspline(basis_exact, i, x)
        end
        @test splinevalue(bspline, x, leftknot=leftknot) == spl
        @test_throws ArgumentError splinevalue(bspline, x, leftknot=nothing)
        for N = 0:5
            spl = splinevalue(bspline, x, Derivative(N))
            if isexact
                @test spl == splinevalue_exact_bspline(basis_exact, i, x, N)
            else
                @test spl ≈ₑₗ splinevalue_exact_bspline(basis_exact, i, x, N)
            end
            @test splinevalue(bspline, x, Derivative(N), leftknot=leftknot) == spl
            @test_throws ArgumentError splinevalue(bspline, x, Derivative(N), leftknot=nothing)
        end
    end

    function test_splinevalue_bsplines(basis::AbstractBSplineBasis, basis_exact, xtype, isexact)
        for (i, bspline) = enumerate(basis)
            a, b = support(bspline)
            for x = xrange(xtype, a, b, 10)
                test_splinevalue_bsplines_nonzero(bspline, basis_exact, i, x, isexact)
            end
            a′, b′ = support(basis)
            for x = xrange(xtype, a′, b′, 5)
                test_splinevalue_bsplines_nonzero(bspline, basis_exact, i, x, isexact)
            end
            test_splinevalue_zero(bspline, xtype(floor(a′-1)), isexact)
            test_splinevalue_zero(bspline, xtype(ceil(b′+1)), isexact)
        end
    end

    function xrange(xtype, a, b, length)
        if xtype <: Integer
            xs = ceil(xtype, a):floor(xtype, b)
        elseif xtype <: AbstractFloat
            xa = xtype(a)
            xb = xtype(b)
            xa < a && (xa = nextfloat(xa))
            xb > b && (xb = prevfloat(xb))
            xs = range(xa, stop=xb, length=length)
        else
            xs = range(xtype(a), stop=xtype(b), length=length)
        end
        xs
    end

    # b1 = BSplineBasis(1, 0:5)
    test_splinevalue_bsplines(b1, b1_exact, Int,              true)
    test_splinevalue_bsplines(b1, b1_exact, BigInt,           true)
    test_splinevalue_bsplines(b1, b1_exact, Float32,          true)
    test_splinevalue_bsplines(b1, b1_exact, Rational{Int},    true)
    # b2 = BSplineBasis(2, big(0):1:big(5))
    test_splinevalue_bsplines(b2, b2_exact, BigInt,           false)
    test_splinevalue_bsplines(b2, b2_exact, BigFloat,         false)
    test_splinevalue_bsplines(b2, b2_exact, Rational{Int},    true)
    test_splinevalue_bsplines(b2, b2_exact, Rational{BigInt}, true)
    # b3 = BSplineBasis(3, 0//1:5//1)
    test_splinevalue_bsplines(b3, b3_exact, Int,              true)
    test_splinevalue_bsplines(b3, b3_exact, Float32,          false)
    test_splinevalue_bsplines(b3, b3_exact, Float64,          false)
    test_splinevalue_bsplines(b3, b3_exact, Rational{BigInt}, true)
    # b4 = BSplineBasis(4, Int128[0,1,2,3,4,5])
    test_splinevalue_bsplines(b4, b4_exact, BigInt,           false)
    test_splinevalue_bsplines(b4, b4_exact, Float64,          false)
    test_splinevalue_bsplines(b4, b4_exact, Rational{Int},    true)
    test_splinevalue_bsplines(b4, b4_exact, Rational{BigInt}, true)
    # b5 = BSplineBasis(3, [-5.0, -3.5, -2.2, 1.0, 2.0, 3.6])
    test_splinevalue_bsplines(b5, b5_exact, Int,              false)
    test_splinevalue_bsplines(b5, b5_exact, Float32,          false)
    test_splinevalue_bsplines(b5, b5_exact, Float64,          false)
    test_splinevalue_bsplines(b5, b5_exact, BigFloat,         false)
    # b6 = BSplineBasis(4, Rational{BigInt}[-5//1, -7//2, -11//5, 1//1, 2//1, 18//5])
    test_splinevalue_bsplines(b6, b6_exact, Int,              true)
    test_splinevalue_bsplines(b6, b6_exact, Float32,          false)
    test_splinevalue_bsplines(b6, b6_exact, Float64,          false)
    test_splinevalue_bsplines(b6, b6_exact, Rational{Int},    true)
    # b7 = BSplineBasis(3, Int8[-6,-4,-2,0,0,1,2])
    test_splinevalue_bsplines(b4, b4_exact, Int32,            false)
    test_splinevalue_bsplines(b4, b4_exact, BigFloat,         false)
    test_splinevalue_bsplines(b4, b4_exact, Rational{Int16},  true)
    test_splinevalue_bsplines(b4, b4_exact, Rational{Int128}, true)

    function test_splinevalue_Inf(spline::Spline)
        for T = (Float16, Float32, Float64, BigFloat, Rational{Int}, Rational{BigInt})
            test_splinevalue_zero(spline, T(Inf), true)
            test_splinevalue_zero(spline, T(-Inf), true)
        end
    end

    function test_splinevalue_NaN(spline::Spline)
        for T = (Float16, Float32, Float64, BigFloat)
            x = T(NaN)
            @test isnan(splinevalue(spline, x))
            @test isnan(splinevalue(spline, x, leftknot=nothing))
            @test_throws ArgumentError splinevalue(spline, x, leftknot=1)
            for N = 0:5
                @test isnan(splinevalue(spline, x, Derivative(N)))
                @test isnan(splinevalue(spline, x, Derivative(N), leftknot=nothing))
                @test_throws ArgumentError splinevalue(spline, x, Derivative(N), leftknot=1)
            end
        end
    end

    for basis = (b1, b2, b3, b4, b5, b6, b7)
        test_splinevalue_Inf(basis[1])
        test_splinevalue_NaN(basis[1])
    end

    function splinevalue_exact(splines, coeffs, x, nderiv=0)
        if nderiv ≥ size(splines,2)
            0
        else
            bigx = maybebig(x)
            sum(coeffs[i]*(splines[i,nderiv+1])(bigx) for i=eachindex(coeffs))
        end
    end

    function test_splinevalue_nonzero(spline::Spline, basis_exact, x, isexact)
        leftknot = intervalindex(basis(spline), x)
        leftknot === nothing && error("this function is only suitable for x for which the B-spline is non-zero")
        c = coeffs(spline)
        spl = splinevalue(spline, x)
        if isexact
            @test spl == splinevalue_exact(basis_exact, c, x)
        else
            @test spl ≈ₑₗ splinevalue_exact(basis_exact, c, x)
        end
        @test splinevalue(spline, x, leftknot=leftknot) == spl
        @test_throws ArgumentError splinevalue(spline, x, leftknot=nothing)
        for N = 0:5
            spl = splinevalue(spline, x, Derivative(N))
            if isexact
                @test spl == splinevalue_exact(basis_exact, c, x, N)
            else
                @test spl ≈ₑₗ splinevalue_exact(basis_exact, c, x, N)
            end
            @test splinevalue(spline, x, Derivative(N), leftknot=leftknot) == spl
            @test_throws ArgumentError splinevalue(spline, x, Derivative(N), leftknot=nothing)
        end
    end

    function test_splinevalue(spline::Spline, basis_exact, xtype, isexact)
        a, b = support(spline)
        for x = xrange(xtype, a, b, 10)
            test_splinevalue_nonzero(spline, basis_exact, x, isexact)
        end
        test_splinevalue_zero(spline, xtype(floor(a-1)), isexact)
        test_splinevalue_zero(spline, xtype(ceil(b+1)), isexact)
    end

    # b1 = BSplineBasis(1, 0:5)
    s1 = Spline(b1, 1.0:5.0)
    test_splinevalue(s1, b1_exact, Int,              true)
    test_splinevalue(s1, b1_exact, BigInt,           true)
    test_splinevalue(s1, b1_exact, Float64,          true)
    test_splinevalue(s1, b1_exact, Rational{BigInt}, true)
    # b2 = BSplineBasis(2, big(0):1:big(5))
    s2 = Spline(b2, 6.0:-1:1.0)
    test_splinevalue(s2, b2_exact, BigInt,           false)
    test_splinevalue(s2, b2_exact, Float32,          false)
    test_splinevalue(s2, b2_exact, BigFloat,         false)
    test_splinevalue(s2, b2_exact, Rational{Int},    false)
    # b3 = BSplineBasis(3, 0//1:5//1)
    s3 = Spline(b3, -30:10:30)
    test_splinevalue(s3, b3_exact, Int,              true)
    test_splinevalue(s3, b3_exact, Float32,          false)
    test_splinevalue(s3, b3_exact, BigFloat,         false)
    test_splinevalue(s3, b3_exact, Rational{BigInt}, true)
    # b4 = BSplineBasis(4, Int128[0,1,2,3,4,5])
    s4 = Spline(b4, 1:8)
    test_splinevalue(s4, b4_exact, BigInt,           false)
    test_splinevalue(s4, b4_exact, Float64,          false)
    test_splinevalue(s4, b4_exact, Rational{Int},    true)
    test_splinevalue(s4, b4_exact, Rational{BigInt}, true)
    # b5 = BSplineBasis(3, [-5.0, -3.5, -2.2, 1.0, 2.0, 3.6])
    s5 = Spline(b5, ones(Float64, 7))
    test_splinevalue(s5, b5_exact, Int,              false)
    test_splinevalue(s5, b5_exact, Float32,          false)
    test_splinevalue(s5, b5_exact, Float64,          false)
    test_splinevalue(s5, b5_exact, BigFloat,         false)
    # b6 = BSplineBasis(4, Rational{BigInt}[-5//1, -7//2, -11//5, 1//1, 2//1, 18//5])
    s6 = Spline(b6, [1//i for i=eachindex(b6)])
    test_splinevalue(s6, b6_exact, Int,              true)
    test_splinevalue(s6, b6_exact, Float32,          false)
    test_splinevalue(s6, b6_exact, Float64,          false)
    test_splinevalue(s6, b6_exact, Rational{Int},    true)
    # b7 = BSplineBasis(3, Int8[-6,-4,-2,0,0,1,2])
    s7 = Spline(b7, view(1:10, 2:9))
    test_splinevalue(s7, b7_exact, Int16,            false)
    test_splinevalue(s7, b7_exact, BigFloat,         false)
    test_splinevalue(s7, b7_exact, Rational{Int16},  true)
    test_splinevalue(s7, b7_exact, Rational{Int128}, true)

    for spline = (s1, s2, s3, s4, s5, s6, s7)
        test_splinevalue_Inf(spline)
        test_splinevalue_NaN(spline)
        leftknot = intervalindex(basis(spline), 1.5)
        spl = splinevalue(spline, 1.5, leftknot=leftknot)
        drv = splinevalue(spline, 1.5, Derivative(1), leftknot=leftknot)
        @test spline(1.5) == spl
        @test spline(1.5, leftknot=leftknot) == spl
        @test spline(1.5, Derivative(1)) == drv
        @test spline(1.5, Derivative(1), leftknot=leftknot) == drv
        @test_throws ArgumentError spline(1.5, leftknot=leftknot+1)
        @test_throws ArgumentError spline(1.5, leftknot=leftknot-1)
        @test_throws ArgumentError spline(1.5, leftknot=nothing)
        @test_throws ArgumentError spline(1.5, Derivative(1), leftknot=leftknot+1)
        @test_throws ArgumentError spline(1.5, Derivative(1), leftknot=leftknot-1)
        @test_throws ArgumentError spline(1.5, Derivative(1), leftknot=nothing)
    end
end
