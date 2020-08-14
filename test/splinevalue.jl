@time @testset "splinevalue" begin

    # b1 = BSplineBasis(1, 0:5)
    s1 = Spline(b1, 1.0:5.0)
    # b2 = BSplineBasis(2, big(0):1:big(5))
    s2 = Spline(b2, 6.0:-1:1.0)
    # b3 = BSplineBasis(3, 0//1:5//1)
    s3 = Spline(b3, -30:10:30)
    # b4 = BSplineBasis(4, Int128[0,1,2,3,4,5])
    s4 = Spline(b4, 1:8)
    # b5 = BSplineBasis(3, [-5.0, -3.5, -2.2, 1.0, 2.0, 3.6])
    s5 = Spline(b5, ones(Float64, 7))
    # b6 = BSplineBasis(4, Rational{Int64}[-5//1, -7//2, -11//5, 1//1, 2//1, 18//5])
    s6 = Spline(b6, [1//i for i=eachindex(b6)])
    # b7 = BSplineBasis(3, Int8[-6,-4,-2,0,0,1,2])
    s7 = Spline(b7, view(1:10, 2:9))

    @testset "Return types" begin
        function test_splinevalue_returntype(spline, xtype, bspline_type, leftknot=intervalindex(basis(spline), 1))
            x = one(xtype)
            @test @inferred(splinevalue(spline, x, leftknot=leftknot)) isa bspline_type
            @test @inferred(splinevalue(spline, x, Derivative(1), leftknot=leftknot)) isa bspline_type
        end
        function test_splinevalue_returntype_workspace(spline, leftknot=intervalindex(basis(spline), 1))
            for (T,N) = ((Float32,0), (Float64,1), (BigFloat,2), (Rational{Int64},3), (Rational{BigInt},4))
                work = Vector{T}(undef, order(spline))
                @test @inferred(splinevalue(spline, 1, Derivative(N), leftknot=leftknot, workspace=work)) isa T
            end
        end

        # b1 = BSplineBasis(1, 0:5)
        b1_bspline_types = [(Int,Float64), (BigInt,BigFloat),
                            (Float32,Float32), (Float64,Float64), (BigFloat,BigFloat),
                            (Rational{Int},Rational{Int}), (Rational{BigInt},Rational{BigInt})]
        # b2 = BSplineBasis(2, big(0):1:big(5))
        b2_bspline_types = [(Int,BigFloat), (Float32,BigFloat), (Rational{Int},Rational{BigInt})]
        # b3 = BSplineBasis(3, 0//1:5//1)
        b3_bspline_types = [(BigInt,Rational{BigInt}), (Float64,Float64), (Rational{BigInt},Rational{BigInt})]
        # b4 = BSplineBasis(4, Int128[0,1,2,3,4,5])
        b4_bspline_types = [(Int,Float64), (BigFloat,BigFloat), (Rational{Int},Rational{Int128})]
        # b5 = BSplineBasis(3, [-5.0, -3.5, -2.2, 1.0, 2.0, 3.6])
        b5_bspline_types = [(BigInt,BigFloat), (Float32,Float64), (Rational{BigInt},BigFloat)]
        # b6 = BSplineBasis(4, Rational{Int64}[-5//1, -7//2, -11//5, 1//1, 2//1, 18//5])
        b6_bspline_types = [(Int,Rational{Int64}), (Float64,Float64), (Rational{Int},Rational{Int64})]
        # b7 = BSplineBasis(3, Int8[-6,-4,-2,0,0,1,2])
        b7_bspline_types = [(BigInt,BigFloat), (BigFloat,BigFloat), (Rational{Int},Rational{Int})]
        for (basis, types) = [(b1, b1_bspline_types), (b2, b2_bspline_types),
                              (b3, b3_bspline_types), (b4, b4_bspline_types),
                              (b5, b5_bspline_types), (b6, b6_bspline_types),
                              (b7, b7_bspline_types)]
            bspline = basis[1]
            leftknot = intervalindex(basis, 1)
            test_splinevalue_returntype_workspace(bspline, leftknot)
            for (xtype, bspline_type) in types
                test_splinevalue_returntype(bspline, xtype, bspline_type, leftknot)
            end
        end

        # b1 = BSplineBasis(1, 0:5)
        b1_coeffs_and_types = [(1:5,Int,Float64),
                               (1:5,Float32,Float32)]
        # b2 = BSplineBasis(2, big(0):1:big(5))
        b2_coeffs_and_types = [(1//1:6//1,Int8,Rational{BigInt}),
                               (1:6,Int,BigFloat)]
        # b3 = BSplineBasis(3, 0//1:5//1)
        b3_coeffs_and_types = [(1:7,Float32,Float32),
                               (Float32.(1:7),Float64,Float64)]
        # b4 = BSplineBasis(4, Int128[0,1,2,3,4,5])
        b4_coeffs_and_types = [(Float64.(1:8),Rational{Int},Float64),
                               (1:8,Rational{Int},Rational{Int128})]
        # b5 = BSplineBasis(3, [-5.0, -3.5, -2.2, 1.0, 2.0, 3.6])
        b5_coeffs_and_types = [(Float32.(1:7),Float32,Float64),
                               (big(1):big(7),Float64,BigFloat)]
        # b6 = BSplineBasis(4, Rational{Int64}[-5//1, -7//2, -11//5, 1//1, 2//1, 18//5])
        b6_coeffs_and_types = [(1:8,Int,Rational{Int64}),
                               (1:8,Float64,Float64)]
        # b7 = BSplineBasis(3, Int8[-6,-4,-2,0,0,1,2])
        b7_coeffs_and_types = [(Float64.(1:8),Rational{Int},Float64),
                               (BigFloat.(1:8),Int,BigFloat)]
        for (basis, coeffs_and_types) = [(b1, b1_coeffs_and_types), (b2, b2_coeffs_and_types),
                                         (b3, b3_coeffs_and_types), (b4, b4_coeffs_and_types),
                                         (b5, b5_coeffs_and_types), (b6, b6_coeffs_and_types),
                                         (b7, b7_coeffs_and_types)]
            leftknot = intervalindex(basis, 1)
            for (coeffs, xtype, bspline_type) in coeffs_and_types
                spline = Spline(basis, coeffs)
                test_splinevalue_returntype(spline, xtype, bspline_type, leftknot)
            end
        end
        for spline = [s1, s2, s3, s4, s5, s6, s7]
            test_splinevalue_returntype_workspace(spline)
        end
    end

    function test_splinevalue_zero(spline::Spline, x)
        work = zeros(order(spline))
        for N = 0:3
            @test splinevalue(spline, x, Derivative(N)) == 0
            @test splinevalue(spline, x, Derivative(N), workspace=work) == 0
            @test splinevalue(spline, x, Derivative(N), leftknot=nothing) == 0
            @test splinevalue(spline, x, Derivative(N), leftknot=nothing, workspace=work) == 0
            @test_throws ArgumentError splinevalue(spline, x, Derivative(N), leftknot=1)
            @test_throws ArgumentError splinevalue(spline, x, Derivative(N), leftknot=1, workspace=work)
        end
    end
    test_splinevalue_zero(b1[1], Inf16)
    test_splinevalue_zero(b2[2], -Inf32)
    test_splinevalue_zero(b3[3], Inf64)
    test_splinevalue_zero(b4[4], 1//0)
    test_splinevalue_zero(b5[5], -1//0)
    test_splinevalue_zero(b6[6], -10.0)
    test_splinevalue_zero(b7[7], 100)

    function test_splinevalue_nan(spline::Spline, x)
        work = zeros(order(spline))
        for N = 0:3
            @test isnan(splinevalue(spline, x, Derivative(N)))
            @test isnan(splinevalue(spline, x, Derivative(N), workspace=work))
        end
        if isnan(x)
            for N = 0:3
                @test isnan(splinevalue(spline, x, Derivative(N), leftknot=nothing))
                @test isnan(splinevalue(spline, x, Derivative(N), leftknot=nothing, workspace=work))
                @test_throws ArgumentError splinevalue(spline, x, Derivative(N), leftknot=1)
                @test_throws ArgumentError splinevalue(spline, x, Derivative(N), leftknot=1, workspace=work)
            end
        end
    end
    test_splinevalue_nan(b1[1], NaN16)
    test_splinevalue_nan(b2[2], NaN32)
    test_splinevalue_nan(b3[3], NaN64)
    test_splinevalue_nan(b4[4], big(NaN))

    function splinevalue_exact(splines, coeffs, x, nderiv=0)
        if nderiv ≥ size(splines,2)
            0
        elseif coeffs isa BSplines.StandardBasisVector
            (splines[coeffs.index,nderiv+1])(maybebig(x))
        else
            bigx = maybebig(x)
            sum(coeffs[i]*(splines[i,nderiv+1])(bigx) for i=eachindex(coeffs))
        end
    end
    # Enlarge x for calculating exact B-splines:
    # * for Integer/Rational: use 64-bit to avoid overflow in rational arithmetic
    # * for floating-point: use BigFloat to minimize error
    maybebig(x::Integer) = Int64(x)
    maybebig(x::Rational) = Rational{Int64}(x)
    maybebig(x::AbstractFloat) = big(x)

    function test_splinevalue_nonzero(spline::Spline, basis_exact, x, isexact)
        leftknot = intervalindex(basis(spline), x)
        leftknot === nothing && error("this function is only suitable for x which are inside the support of the basis")
        work = Array{isexact ? Rational{Int64} : Float64}(undef, order(spline))
        for N = 0:3
            spl = splinevalue(spline, x, Derivative(N), leftknot=leftknot)
            exact = splinevalue_exact(basis_exact, coeffs(spline), x, N)
            if isexact
                @test splinevalue(spline, x, Derivative(N), leftknot=leftknot) == exact
                @test splinevalue(spline, x, Derivative(N), leftknot=leftknot, workspace=work) == exact
            else
                @test splinevalue(spline, x, Derivative(N), leftknot=leftknot) ≈ₑₗ exact
                @test splinevalue(spline, x, Derivative(N), leftknot=leftknot, workspace=work) ≈ₑₗ exact
            end
        end
    end

    function test_splinevalue_bsplines(basis::BSplineBasis, basis_exact, xtype, isexact)
        for bspline = basis
            a, b = support(bspline)
            for x = xrange(xtype, a, b, 5)
                test_splinevalue_nonzero(bspline, basis_exact, x, isexact)
            end
            a′, b′ = support(basis)
            for x = xrange(xtype, a′, b′, 5)
                test_splinevalue_nonzero(bspline, basis_exact, x, isexact)
            end
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
    test_splinevalue_bsplines(b1, b1_exact, BigInt,           true)
    test_splinevalue_bsplines(b1, b1_exact, Float32,          true)
    # b2 = BSplineBasis(2, big(0):1:big(5))
    test_splinevalue_bsplines(b2, b2_exact, BigFloat,         false)
    test_splinevalue_bsplines(b2, b2_exact, Rational{Int32},  true)
    # b3 = BSplineBasis(3, 0//1:5//1)
    test_splinevalue_bsplines(b3, b3_exact, Float64,          false)
    test_splinevalue_bsplines(b3, b3_exact, Rational{BigInt}, true)
    # b4 = BSplineBasis(4, Int128[0,1,2,3,4,5])
    test_splinevalue_bsplines(b4, b4_exact, BigInt,           false)
    test_splinevalue_bsplines(b4, b4_exact, Rational{Int64},  true)
    # b5 = BSplineBasis(3, [-5.0, -3.5, -2.2, 1.0, 2.0, 3.6])
    test_splinevalue_bsplines(b5, b5_exact, Int,              false)
    test_splinevalue_bsplines(b5, b5_exact, Float32,          false)
    # b6 = BSplineBasis(4, Rational{Int64}[-5//1, -7//2, -11//5, 1//1, 2//1, 18//5])
    test_splinevalue_bsplines(b6, b6_exact, Int,              true)
    test_splinevalue_bsplines(b6, b6_exact, Float64,          false)
    # b7 = BSplineBasis(3, Int8[-6,-4,-2,0,0,1,2])
    test_splinevalue_bsplines(b7, b7_exact, BigFloat,         false)
    test_splinevalue_bsplines(b7, b7_exact, Rational{Int16},  true)

    function test_splinevalue(spline::Spline, basis_exact, xtype, isexact)
        a, b = support(spline)
        for x = xrange(xtype, a, b, 10)
            test_splinevalue_nonzero(spline, basis_exact, x, isexact)
        end
    end

    # b1 = BSplineBasis(1, 0:5)
    # s1 = Spline(b1, 1.0:5.0)
    test_splinevalue(s1, b1_exact, BigInt,           true)
    test_splinevalue(s1, b1_exact, Float64,          true)
    # b2 = BSplineBasis(2, big(0):1:big(5))
    # s2 = Spline(b2, 6.0:-1:1.0)
    test_splinevalue(s2, b2_exact, BigFloat,         false)
    test_splinevalue(s2, b2_exact, Rational{Int32},  false)
    # b3 = BSplineBasis(3, 0//1:5//1)
    # s3 = Spline(b3, -30:10:30)
    test_splinevalue(s3, b3_exact, Float32,          false)
    test_splinevalue(s3, b3_exact, Rational{BigInt}, true)
    # b4 = BSplineBasis(4, Int128[0,1,2,3,4,5])
    # s4 = Spline(b4, 1:8)
    test_splinevalue(s4, b4_exact, BigInt,           false)
    test_splinevalue(s4, b4_exact, Rational{Int64},  true)
    # b5 = BSplineBasis(3, [-5.0, -3.5, -2.2, 1.0, 2.0, 3.6])
    # s5 = Spline(b5, ones(Float64, 7))
    test_splinevalue(s5, b5_exact, Int,              false)
    test_splinevalue(s5, b5_exact, Float32,          false)
    # b6 = BSplineBasis(4, Rational{Int64}[-5//1, -7//2, -11//5, 1//1, 2//1, 18//5])
    # s6 = Spline(b6, [1//i for i=eachindex(b6)])
    test_splinevalue(s6, b6_exact, Int,              true)
    test_splinevalue(s6, b6_exact, Float64,          false)
    # b7 = BSplineBasis(3, Int8[-6,-4,-2,0,0,1,2])
    # s7 = Spline(b7, view(1:10, 2:9))
    test_splinevalue(s7, b7_exact, BigFloat,         false)
    test_splinevalue(s7, b7_exact, Rational{Int16},  true)

    test_splinevalue_zero(s1, Inf16)
    test_splinevalue_zero(s2, -Inf32)
    test_splinevalue_zero(s3, Inf64)
    test_splinevalue_zero(s4, 1//0)
    test_splinevalue_zero(s5, -1//0)
    test_splinevalue_zero(s6, -10.0)
    test_splinevalue_zero(s7, 100)
    test_splinevalue_nan(s1, NaN16)
    test_splinevalue_nan(s2, NaN32)
    test_splinevalue_nan(s3, NaN64)
    test_splinevalue_nan(s4, big(NaN))
    test_splinevalue_nan(Spline(b4, fill(NaN, length(b4))), 1.5)

    k = order(s5);
    @test_throws DimensionMismatch splinevalue(s5, 0, workspace=zeros(k+1))
    @test_throws DimensionMismatch splinevalue(s5, 0, workspace=zeros(k-1))
    @test_throws DimensionMismatch splinevalue(s5, 0, workspace=zeros(1, k))
    @test_throws DimensionMismatch splinevalue(s5, 0, workspace=zeros(k, 1))
    @test_throws DimensionMismatch splinevalue(s5, 0, workspace=zeros(k, k))
    leftknot = intervalindex(basis(s5), 0)
    @test_throws ArgumentError splinevalue(s5, 0, leftknot=leftknot-1)
    @test_throws ArgumentError splinevalue(s5, 0, leftknot=leftknot+1)
    @test_throws ArgumentError splinevalue(s5, 0, leftknot=nothing)
    @test_throws ArgumentError splinevalue(s5, 0, Derivative(0), leftknot=leftknot-1)
    @test_throws ArgumentError splinevalue(s5, 0, Derivative(1), leftknot=leftknot+1)
    @test_throws ArgumentError splinevalue(s5, 0, Derivative(2), leftknot=nothing)
    @test s5(0) == s5(0, leftknot=leftknot) == splinevalue(s5, 0, leftknot=leftknot)
    @test s5(0, Derivative(1)) == s5(0, Derivative(1), leftknot=leftknot) == splinevalue(s5, 0, Derivative(1), leftknot=leftknot)
end
