using Test
using BSplines
using BSplines: KnotVector, StandardBasisVector, check_intervalindex
using LinearAlgebra: ldiv!, lmul!, rdiv!, rmul!
using OffsetArrays: OffsetArray

include("knotvector.jl")
include("standardbasisvector.jl")

@testset "intervalindex" begin
    function test_intervalindex(basis::BSplineBasis, x, leftknot)
        @test intervalindex(basis, x) === leftknot
        for start = eachindex(knots(basis))
            @test intervalindex(basis, x, start) === leftknot
        end
        @test_throws BoundsError intervalindex(basis, x, firstindex(knots(basis))-1)
        @test_throws BoundsError intervalindex(basis, x, lastindex(knots(basis))+1)
    end
    function test_intervalindex(vec, x, leftknot)
        @test intervalindex(vec, x) === leftknot
        for start = eachindex(vec)
            @test intervalindex(vec, x, start) === leftknot
        end
        @test_throws BoundsError intervalindex(vec, x, firstindex(vec)-1)
        @test_throws BoundsError intervalindex(vec, x, lastindex(vec)+1)
    end
    test_intervalindex(-10:10, -Inf32,         nothing)
    test_intervalindex(-10:10, -11,            nothing)
    test_intervalindex(-10:10, -10.0,          1)
    test_intervalindex(-10:10, -9.5f0,         1)
    test_intervalindex(-10:10, big(π),         14)
    test_intervalindex(-10:10, prevfloat(9.0), 19)
    test_intervalindex(-10:10, 9,              20)
    test_intervalindex(-10:10, 10.0,           20)
    test_intervalindex(-10:10, 10.1,           nothing)
    test_intervalindex(-10:10, NaN16,          nothing)
    test_intervalindex(-10:10, big(Inf),       nothing)
    test_intervalindex(-5//1:1//2:5//1, -Inf16,         nothing)
    test_intervalindex(-5//1:1//2:5//1, -5.1,           nothing)
    test_intervalindex(-5//1:1//2:5//1, -5.0,           1)
    test_intervalindex(-5//1:1//2:5//1, -4.5f0,         2)
    test_intervalindex(-5//1:1//2:5//1, 1.6,            14)
    test_intervalindex(-5//1:1//2:5//1, prevfloat(4.5), 19)
    test_intervalindex(-5//1:1//2:5//1, 9//2,           20)
    test_intervalindex(-5//1:1//2:5//1, 5.0,            20)
    test_intervalindex(-5//1:1//2:5//1, 5.1,            nothing)
    test_intervalindex(-5//1:1//2:5//1, NaN32,          nothing)
    test_intervalindex(-5//1:1//2:5//1, Inf64,          nothing)
    test_intervalindex(10_000_000.1:0.001:10_000_000.2, big(-Inf),                 nothing)
    test_intervalindex(10_000_000.1:0.001:10_000_000.2, 10_000_000,                nothing)
    test_intervalindex(10_000_000.1:0.001:10_000_000.2, 10_000_000.1,              1)
    test_intervalindex(10_000_000.1:0.001:10_000_000.2, prevfloat(10_000_000.123), 23)
    test_intervalindex(10_000_000.1:0.001:10_000_000.2, 10_000_000.123,            24)
    test_intervalindex(10_000_000.1:0.001:10_000_000.2, 10_000_000.15,             51)
    test_intervalindex(10_000_000.1:0.001:10_000_000.2, 10_000_000.1987,           99)
    test_intervalindex(10_000_000.1:0.001:10_000_000.2, 10_000_000.199,            100)
    test_intervalindex(10_000_000.1:0.001:10_000_000.2, 10_000_000.19999,          100)
    test_intervalindex(10_000_000.1:0.001:10_000_000.2, 10_000_000.2,              100)
    test_intervalindex(10_000_000.1:0.001:10_000_000.2, nextfloat(10_000_000.2),   nothing)
    test_intervalindex(10_000_000.1:0.001:10_000_000.2, big(NaN),                  nothing)
    test_intervalindex(10_000_000.1:0.001:10_000_000.2, Inf16,                     nothing)
    test_intervalindex([0, 1, 2, 3, 4, 4, 5, 6, 6], -Inf16, nothing)
    test_intervalindex([0, 1, 2, 3, 4, 4, 5, 6, 6], -1,     nothing)
    test_intervalindex([0, 1, 2, 3, 4, 4, 5, 6, 6], 0.0,    1)
    test_intervalindex([0, 1, 2, 3, 4, 4, 5, 6, 6], 0.5,    1)
    test_intervalindex([0, 1, 2, 3, 4, 4, 5, 6, 6], 1.0,    2)
    test_intervalindex([0, 1, 2, 3, 4, 4, 5, 6, 6], 4.0,    6)
    test_intervalindex([0, 1, 2, 3, 4, 4, 5, 6, 6], 4.5,    6)
    test_intervalindex([0, 1, 2, 3, 4, 4, 5, 6, 6], 6.0,    7)
    test_intervalindex([0, 1, 2, 3, 4, 4, 5, 6, 6], 6.1,    nothing)
    test_intervalindex([0, 1, 2, 3, 4, 4, 5, 6, 6], Inf64,  nothing)
    test_intervalindex([0, 1, 2, 3, 4, 4, 5, 6, 6], NaN,    nothing)
    test_intervalindex([5, 5, 5, 5, 5], 5, nothing)
    test_intervalindex(KnotVector(-5:5, 2), -Inf32,         nothing)
    test_intervalindex(KnotVector(-5:5, 2), -10,            nothing)
    test_intervalindex(KnotVector(-5:5, 2), -5,             3)
    test_intervalindex(KnotVector(-5:5, 2), -4.5,           3)
    test_intervalindex(KnotVector(-5:5, 2), 0,              8)
    test_intervalindex(KnotVector(-5:5, 2), 4.5,            12)
    test_intervalindex(KnotVector(-5:5, 2), 5,              12)
    test_intervalindex(KnotVector(-5:5, 2), nextfloat(5.0), nothing)
    test_intervalindex(KnotVector(-5:5, 2), NaN32,          nothing)
    test_intervalindex(KnotVector([0, 1, 2, 3, 4, 4, 5, 6], 4, 6), -1,  nothing)
    test_intervalindex(KnotVector([0, 1, 2, 3, 4, 4, 5, 6], 4, 6), 0.0, 5)
    test_intervalindex(KnotVector([0, 1, 2, 3, 4, 4, 5, 6], 4, 6), 0.5, 5)
    test_intervalindex(KnotVector([0, 1, 2, 3, 4, 4, 5, 6], 4, 6), 1.0, 6)
    test_intervalindex(KnotVector([0, 1, 2, 3, 4, 4, 5, 6], 4, 6), 4.0, 10)
    test_intervalindex(KnotVector([0, 1, 2, 3, 4, 4, 5, 6], 4, 6), 4.5, 10)
    test_intervalindex(KnotVector([0, 1, 2, 3, 4, 4, 5, 6], 4, 6), 6.0, 11)
    test_intervalindex(KnotVector([0, 1, 2, 3, 4, 4, 5, 6], 4, 6), 6.1, nothing)
    test_intervalindex(BSplineBasis(3, 0:5), -Inf, nothing)
    test_intervalindex(BSplineBasis(3, 0:5), -0.5, nothing)
    test_intervalindex(BSplineBasis(3, 0:5), 0//1, 3)
    test_intervalindex(BSplineBasis(3, 0:5), 0.5,  3)
    test_intervalindex(BSplineBasis(3, 0:5), 1,    4)
    test_intervalindex(BSplineBasis(3, 0:5), 4.0,  7)
    test_intervalindex(BSplineBasis(3, 0:5), 4.5,  7)
    test_intervalindex(BSplineBasis(3, 0:5), 5.0,  7)
    test_intervalindex(BSplineBasis(3, 0:5), 5.1,  nothing)
    test_intervalindex(BSplineBasis(3, 0:5), Inf,  nothing)
    test_intervalindex(BSplineBasis(3, 0:5), NaN,  nothing)
    test_intervalindex(BSplineBasis(4, [1,2,3,3,4,5,6]), 0.5,  nothing)
    test_intervalindex(BSplineBasis(4, [1,2,3,3,4,5,6]), 1,    4)
    test_intervalindex(BSplineBasis(4, [1,2,3,3,4,5,6]), 3,    7)
    test_intervalindex(BSplineBasis(4, [1,2,3,3,4,5,6]), 3.5,  7)
    test_intervalindex(BSplineBasis(4, [1,2,3,3,4,5,6]), 4.5,  8)
    test_intervalindex(BSplineBasis(4, [1,2,3,3,4,5,6]), 5.0,  9)
    test_intervalindex(BSplineBasis(4, [1,2,3,3,4,5,6]), 6.0,  9)
    test_intervalindex(BSplineBasis(4, [1,2,3,3,4,5,6]), 6.1,  nothing)
    test_intervalindex(BSplineBasis(4, [1,2,3,3,4,5,6]), NaN,  nothing)
    test_intervalindex(BSplineBasis(4, [1,2,3,3,4,5,6]), Inf,  nothing)
    test_intervalindex(BSplineBasis(4, [1,2,3,3,4,5,6]), -Inf, nothing)
    # check_intervalindex
    for vec = [[0, 1, 2, 3, 4, 4, 5, 6, 6], BSplineBasis(3, 0:5), BSplineBasis(4, [1,2,3,3,4,5,6])]
        for x = Real[-1, -0.5, 0, 0.5, 1, 3, 3.5, 4, 9//2, 5.0, 5.1, 6.0, 6.1, -Inf, Inf, NaN]
            for index = [-1:20; nothing]
                if intervalindex(vec, x) == index
                    @test check_intervalindex(Bool, vec, x, index)
                    @test check_intervalindex(vec, x, index) === nothing
                else
                    @test !check_intervalindex(Bool, vec, x, index)
                    @test_throws ArgumentError check_intervalindex(vec, x, index)
                end
            end
        end
    end
end

@testset "intervalindices" begin
    support = [(1, 3:3), (2, 3:4), (3, 3:5), (4, 4:6), (5, 5:7), (6, 6:7), (7, 7:7)]
    for basis = [BSplineBasis(3, 0:5), BSplineBasis(3, 1//2:1//2:3//1), BSplineBasis(3, 0.0:0.2:1.0)]
        for itr = [intervalindices(basis), intervalindices(basis, :), intervalindices(basis, 1:7),
                   intervalindices(basis, 1), intervalindices(basis, 1, 2)]
            @test itr isa UnitRange{Int}
        end
        @test intervalindices(basis)      == 3:7
        @test intervalindices(basis, :)   == 3:7
        @test intervalindices(basis, 1:7) == 3:7
        @test intervalindices(basis, 1:6) == 3:7
        @test intervalindices(basis, 3:7) == 3:7
        @test intervalindices(basis, 1:5) == 3:7
        @test intervalindices(basis, 2:4) == 3:6
        @test isempty(intervalindices(basis, 1:0))
        @test isempty(intervalindices(basis, 4:3))
        @test isempty(intervalindices(basis, 100:99))
        for (i, irange) = support
            @test intervalindices(basis, i:i) == irange
            @test intervalindices(basis, i) == irange
            for (j, jrange) = support
                @test intervalindices(basis, i, j) == irange ∩ jrange
                for (k, krange) = support
                    @test intervalindices(basis, i, j, k) == irange ∩ jrange ∩ krange
                end
            end
        end
        @test_throws BoundsError intervalindices(basis, 0:1)
        @test_throws BoundsError intervalindices(basis, 1:8)
        @test_throws BoundsError intervalindices(basis, 0)
        @test_throws BoundsError intervalindices(basis, 8)
        @test_throws BoundsError intervalindices(basis, 7, 0)
        @test_throws BoundsError intervalindices(basis, 1, 2, 8)
    end
    support = [(1, [4]), (2, [4,5]), (3, [4,5]), (4, [4,5,7]), (5, [5,7,8]),
               (6, [7,8,9]), (7, [7,8,9]), (8, [8,9]), (9, [9])]
    basis = BSplineBasis(4, [1,2,3,3,4,5,6])
    for itr = [intervalindices(basis), intervalindices(basis, :), intervalindices(basis, 1:7),
               intervalindices(basis, 1), intervalindices(basis, 1, 2)]
        @test eltype(typeof(itr)) === Int
        @test Base.IteratorEltype(typeof(itr)) === Base.HasEltype()
        @test Base.IteratorSize(typeof(itr)) === Base.SizeUnknown()
    end
    @test collect(intervalindices(basis)) == [4,5,7,8,9]
    @test collect(intervalindices(basis, :)) == [4,5,7,8,9]
    @test collect(intervalindices(basis, 1:9)) == [4,5,7,8,9]
    @test collect(intervalindices(basis, 1:8)) == [4,5,7,8,9]
    @test collect(intervalindices(basis, 4:8)) == [4,5,7,8,9]
    @test collect(intervalindices(basis, 1:6)) == [4,5,7,8,9]
    @test collect(intervalindices(basis, 5:6)) == [5,7,8,9]
    @test isempty(intervalindices(basis, 1:0))
    @test isempty(intervalindices(basis, 7:6))
    @test isempty(intervalindices(basis, -100:-101))
    for (i, irange) = support
        @test collect(intervalindices(basis, i:i)) == irange
        @test collect(intervalindices(basis, i)) == irange
        for (j, jrange) = support
            @test collect(intervalindices(basis, i, j)) == irange ∩ jrange
            for (k, krange) = support
                @test collect(intervalindices(basis, i, j, k)) == irange ∩ jrange ∩ krange
            end
        end
    end

    @testset "Iterators.reverse" begin
        support = [(1, 3:3), (2, 3:4), (3, 3:5), (4, 4:6), (5, 5:7), (6, 6:7), (7, 7:7)]
        for basis = [BSplineBasis(3, 0:5), BSplineBasis(3, 1//2:1//2:3//1), BSplineBasis(3, 0.0:0.2:1.0)]
            for itr = Iterators.reverse.([intervalindices(basis), intervalindices(basis, :),
                                          intervalindices(basis, 1:7), intervalindices(basis, 1), 
                                          intervalindices(basis, 1, 2)])
                @test itr isa StepRange{Int}
            end
            @test Iterators.reverse(intervalindices(basis))      == 7:-1:3
            @test Iterators.reverse(intervalindices(basis, :))   == 7:-1:3
            @test Iterators.reverse(intervalindices(basis, 1:7)) == 7:-1:3
            @test Iterators.reverse(intervalindices(basis, 1:6)) == 7:-1:3
            @test Iterators.reverse(intervalindices(basis, 3:7)) == 7:-1:3
            @test Iterators.reverse(intervalindices(basis, 1:5)) == 7:-1:3
            @test Iterators.reverse(intervalindices(basis, 2:4)) == 6:-1:3
            @test isempty(Iterators.reverse(intervalindices(basis, 1:0)))
            @test isempty(Iterators.reverse(intervalindices(basis, 4:3)))
            @test isempty(Iterators.reverse(intervalindices(basis, 100:99)))
            for (i, irange) = support
                @test Iterators.reverse(intervalindices(basis, i:i)) == reverse(irange)
                @test Iterators.reverse(intervalindices(basis, i)) == reverse(irange)
                for (j, jrange) = support
                    @test Iterators.reverse(intervalindices(basis, i, j)) == reverse(irange ∩ jrange)
                    for (k, krange) = support
                        @test Iterators.reverse(intervalindices(basis, i, j, k)) == reverse(irange ∩ jrange ∩ krange)
                    end
                end
            end
        end
        support = [(1, [4]), (2, [4,5]), (3, [4,5]), (4, [4,5,7]), (5, [5,7,8]),
                   (6, [7,8,9]), (7, [7,8,9]), (8, [8,9]), (9, [9])]
        basis = BSplineBasis(4, [1,2,3,3,4,5,6])
        for itr = Iterators.reverse.([intervalindices(basis), intervalindices(basis, :),
                                      intervalindices(basis, 1:7), intervalindices(basis, 1),
                                      intervalindices(basis, 1, 2)])
            @test eltype(typeof(itr)) === Int
            @test Base.IteratorEltype(typeof(itr)) === Base.HasEltype()
            @test Base.IteratorSize(typeof(itr)) === Base.SizeUnknown()
        end
        @test collect(Iterators.reverse(intervalindices(basis))) == [9,8,7,5,4]
        @test collect(Iterators.reverse(intervalindices(basis, :))) == [9,8,7,5,4]
        @test collect(Iterators.reverse(intervalindices(basis, 1:9))) == [9,8,7,5,4]
        @test collect(Iterators.reverse(intervalindices(basis, 1:8))) == [9,8,7,5,4]
        @test collect(Iterators.reverse(intervalindices(basis, 4:8))) == [9,8,7,5,4]
        @test collect(Iterators.reverse(intervalindices(basis, 1:6))) == [9,8,7,5,4]
        @test collect(Iterators.reverse(intervalindices(basis, 5:6))) == [9,8,7,5]
        @test isempty(Iterators.reverse(intervalindices(basis, 1:0)))
        @test isempty(Iterators.reverse(intervalindices(basis, 7:6)))
        @test isempty(Iterators.reverse(intervalindices(basis, -100:-101)))
        for (i, irange) = support
            @test collect(Iterators.reverse(intervalindices(basis, i:i))) == reverse(irange)
            @test collect(Iterators.reverse(intervalindices(basis, i))) == reverse(irange)
            for (j, jrange) = support
                @test collect(Iterators.reverse(intervalindices(basis, i, j))) == reverse(irange ∩ jrange)
                for (k, krange) = support
                    @test collect(Iterators.reverse(intervalindices(basis, i, j, k))) == reverse(irange ∩ jrange ∩ krange)
                end
            end
        end
    end
end

@testset "support" begin
    # Basis
    @test support(BSplineBasis(3, 0:5)) === (0, 5)
    @test support(BSplineBasis(4, 0//1:1//10:5//1)) === (0//1, 5//1)
    @test support(BSplineBasis(5, [-5.0, -3.5, -2.2, 1.0, 2.0, 3.6])) === (-5.0, 3.6)
    # BSpline
    @test support(BSplineBasis(3, 0:5)[1]) === (0, 1)
    @test support(BSplineBasis(3, 0:5)[2]) === (0, 2)
    @test support(BSplineBasis(3, 0:5)[3]) === (0, 3)
    @test support(BSplineBasis(3, 0:5)[4]) === (1, 4)
    @test support(BSplineBasis(3, 0:5)[5]) === (2, 5)
    @test support(BSplineBasis(3, 0:5)[6]) === (3, 5)
    @test support(BSplineBasis(3, 0:5)[7]) === (4, 5)
    @test support(BSplineBasis(5, [-5.0, -3.5, -2.2, 1.0, 2.0, 3.6])[1]) === (-5.0, -3.5)
    @test support(BSplineBasis(5, [-5.0, -3.5, -2.2, 1.0, 2.0, 3.6])[2]) === (-5.0, -2.2)
    @test support(BSplineBasis(5, [-5.0, -3.5, -2.2, 1.0, 2.0, 3.6])[3]) === (-5.0, 1.0)
    @test support(BSplineBasis(5, [-5.0, -3.5, -2.2, 1.0, 2.0, 3.6])[4]) === (-5.0, 2.0)
    @test support(BSplineBasis(5, [-5.0, -3.5, -2.2, 1.0, 2.0, 3.6])[5]) === (-5.0, 3.6)
    @test support(BSplineBasis(5, [-5.0, -3.5, -2.2, 1.0, 2.0, 3.6])[6]) === (-3.5, 3.6)
    @test support(BSplineBasis(5, [-5.0, -3.5, -2.2, 1.0, 2.0, 3.6])[7]) === (-2.2, 3.6)
    @test support(BSplineBasis(5, [-5.0, -3.5, -2.2, 1.0, 2.0, 3.6])[8]) === (1.0, 3.6)
    @test support(BSplineBasis(5, [-5.0, -3.5, -2.2, 1.0, 2.0, 3.6])[9]) === (2.0, 3.6)
    # Spline
    @test support(Spline(BSplineBasis(3, 0:5), ones(7))) === (0, 5)
    @test support(Spline(BSplineBasis(4, 0//1:1//10:5//1), zeros(53))) === (0//1, 5//1)
    @test support(Spline(BSplineBasis(5, [-5.0, -3.5, -2.2, 1.0, 2.0, 3.6]), 1:9)) === (-5.0, 3.6)
end

@testset "Equality and hashing" begin
    bases = [BSplineBasis(3, 0:5),
             BSplineBasis(3, [0.0,1.0,2.0,3.0,4.0,5.0]),
             BSplineBasis(3, 0//1:1//1:5//1),
             BSplineBasis(3, 0:0.5:5),
             BSplineBasis(4, 0:5),
             BSplineBasis(3, [1,2,3,3,4,5])]
    for b1 = bases, b2 = bases
        @test (b1 == b2) == (order(b1) == order(b2) && breakpoints(b1) == breakpoints(b2))
        @test (b1 == b2) == (hash(b1) == hash(b2))
    end
end

@testset "BSplineBasis" begin
    @testset "Construction" begin
        @test_throws DomainError BSplineBasis(0, 0:5) # order must be positive
        @test_throws InexactError BSplineBasis(5.5, 0:5) # order must be integer
        @test_throws TypeError BSplineBasis(5.5, [1+2im, 3+4im, 5+6im]) # breakpoints must be real
        @test_throws ArgumentError BSplineBasis(5, 1:0) # breakpoints must have length ≥ 2
        @test_throws ArgumentError BSplineBasis(5, [1]) # breakpoints must have length ≥ 2
        @test_throws ArgumentError BSplineBasis(5, OffsetArray(1:5, 1)) # breakpoints must not have offset axes
    end

    @testset "Dimensions" begin
        @test length(BSplineBasis(3, 0:5)) === 7
        @test length(BSplineBasis(3, big(0):big(5))) === 7
        @test length(BSplineBasis(5, [1,2,3.5,6,10])) === 8
    end

    @testset "Indices" begin
        @test firstindex(BSplineBasis(3, 0:5)) == 1
        @test firstindex(BSplineBasis(5, [1,2,3.5,6,10])) == 1
        @test lastindex(BSplineBasis(3, 0:5)) == 7
        @test lastindex(BSplineBasis(5, [1,2,3.5,6,10])) == 8
        @test eachindex(BSplineBasis(3, 0:5)) == 1:7
        @test eachindex(BSplineBasis(5, [1,2,3.5,6,10])) == 1:8
        @test keys(BSplineBasis(3, 0:5)) == 1:7
        @test keys(BSplineBasis(5, [1,2,3.5,6,10])) == 1:8
    end

    @testset "checkbounds" begin
        @test checkbounds(BSplineBasis(3, 0:5), :) === nothing
        @test checkbounds(BSplineBasis(3, 0:5), 2) === nothing
        @test checkbounds(BSplineBasis(3, 0:5), 1:7) === nothing
        @test checkbounds(BSplineBasis(3, 0:5), 2:7) === nothing
        @test_throws BoundsError checkbounds(BSplineBasis(3, 0:5), 0)
        @test_throws BoundsError checkbounds(BSplineBasis(3, 0:5), 8)
        @test_throws BoundsError checkbounds(BSplineBasis(3, 0:5), 2:8)
        @test checkbounds(BSplineBasis(5, [1,2,3.5,6,10]), :) === nothing
        @test checkbounds(BSplineBasis(5, [1,2,3.5,6,10]), 2) === nothing
        @test checkbounds(BSplineBasis(5, [1,2,3.5,6,10]), 1:8) === nothing
        @test checkbounds(BSplineBasis(5, [1,2,3.5,6,10]), 2:8) === nothing
        @test_throws BoundsError checkbounds(BSplineBasis(5, [1,2,3.5,6,10]), 0)
        @test_throws BoundsError checkbounds(BSplineBasis(5, [1,2,3.5,6,10]), 9)
        @test_throws BoundsError checkbounds(BSplineBasis(5, [1,2,3.5,6,10]), 2:9)
        @test checkbounds(Bool, BSplineBasis(3, 0:5), :)
        @test checkbounds(Bool, BSplineBasis(3, 0:5), 2)
        @test checkbounds(Bool, BSplineBasis(3, 0:5), 1:7)
        @test checkbounds(Bool, BSplineBasis(3, 0:5), 2:7)
        @test checkbounds(Bool, BSplineBasis(3, 0:5), 0) === false
        @test checkbounds(Bool, BSplineBasis(3, 0:5), 8) === false
        @test checkbounds(Bool, BSplineBasis(3, 0:5), 2:8) === false
        @test checkbounds(Bool, BSplineBasis(5, [1,2,3.5,6,10]), :)
        @test checkbounds(Bool, BSplineBasis(5, [1,2,3.5,6,10]), 2)
        @test checkbounds(Bool, BSplineBasis(5, [1,2,3.5,6,10]), 1:8)
        @test checkbounds(Bool, BSplineBasis(5, [1,2,3.5,6,10]), 2:8)
        @test checkbounds(Bool, BSplineBasis(5, [1,2,3.5,6,10]), 0) === false
        @test checkbounds(Bool, BSplineBasis(5, [1,2,3.5,6,10]), 9) === false
        @test checkbounds(Bool, BSplineBasis(5, [1,2,3.5,6,10]), 2:9) === false
    end

    @testset "Indexing and iteration" begin
        b1 = BSplineBasis(5, 0:5)
        @test eltype(typeof(b1)) === BSpline{BSplineBasis{UnitRange{Int}}}
        @test Base.IteratorEltype(typeof(b1)) === Base.HasEltype()
        @test Base.IteratorSize(typeof(b1)) === Base.HasLength()
        @test collect(b1) == [b1[i] for i=1:9]
        @test @inferred(first(b1)) isa BSpline{BSplineBasis{UnitRange{Int}}}
        @test @inferred(last(b1)) isa BSpline{BSplineBasis{UnitRange{Int}}}
        @test @inferred(b1[1]) isa BSpline{BSplineBasis{UnitRange{Int}}}
        @test first(b1) == Spline(b1, [1,0,0,0,0,0,0,0,0])
        @test last(b1) == Spline(b1, [0,0,0,0,0,0,0,0,1])
        @test b1[5] == Spline(b1, [0,0,0,0,1,0,0,0,0])
        @test b1[end] == Spline(b1, [0,0,0,0,0,0,0,0,1])
        @test_throws BoundsError b1[0]
        @test_throws BoundsError b1[10]
        b2 = BSplineBasis(3, [1,2,3.5,6,10])
        @test eltype(typeof(b2)) === BSpline{BSplineBasis{Vector{Float64}}}
        @test Base.IteratorEltype(typeof(b2)) === Base.HasEltype()
        @test Base.IteratorSize(typeof(b2)) === Base.HasLength()
        @test collect(b2) == [b2[i] for i=1:6]
        @test @inferred(first(b2)) isa BSpline{BSplineBasis{Vector{Float64}}}
        @test @inferred(last(b2)) isa BSpline{BSplineBasis{Vector{Float64}}}
        @test @inferred(b2[1]) isa BSpline{BSplineBasis{Vector{Float64}}}
        @test first(b2) == Spline(b2, [1,0,0,0,0,0])
        @test last(b2) == Spline(b2, [0,0,0,0,0,1])
        @test b2[5] == Spline(b2, [0,0,0,0,1,0])
        @test b2[end] == Spline(b2, [0,0,0,0,0,1])
        @test_throws BoundsError b2[0]
        @test_throws BoundsError b2[7]

        @testset "Iterators.reverse" begin
            b1 = BSplineBasis(5, 0:5)
            rb1 = Iterators.reverse(b1)
            @test eltype(typeof(rb1)) === BSpline{BSplineBasis{UnitRange{Int}}}
            @test Base.IteratorEltype(typeof(rb1)) === Base.HasEltype()
            @test Base.IteratorSize(typeof(rb1)) === Base.HasLength()
            @test collect(rb1) == [b1[i] for i=9:-1:1]
            @test @inferred(first(rb1)) isa BSpline{BSplineBasis{UnitRange{Int}}}
            @test @inferred(last(rb1)) isa BSpline{BSplineBasis{UnitRange{Int}}}
            @test first(rb1) == Spline(b1, [0,0,0,0,0,0,0,0,1])
            @test last(rb1) == Spline(b1, [1,0,0,0,0,0,0,0,0])
            b2 = BSplineBasis(3, [1,2,3.5,6,10])
            rb2 = Iterators.reverse(b2)
            @test eltype(typeof(rb2)) === BSpline{BSplineBasis{Vector{Float64}}}
            @test Base.IteratorEltype(typeof(rb2)) === Base.HasEltype()
            @test Base.IteratorSize(typeof(rb2)) === Base.HasLength()
            @test collect(rb2) == [b2[i] for i=6:-1:1]
            @test @inferred(first(rb2)) isa BSpline{BSplineBasis{Vector{Float64}}}
            @test @inferred(last(rb2)) isa BSpline{BSplineBasis{Vector{Float64}}}
            @test first(rb2) == Spline(b2, [0,0,0,0,0,1])
            @test last(rb2) == Spline(b2, [1,0,0,0,0,0])
        end
    end

    @testset "breakpoints" begin
        @test breakpoints(BSplineBasis(5, 0:5)) == 0:5
        @test breakpoints(BSplineBasis(3, [1,2,3,4,4,5])) == [1,2,3,4,4,5]
    end

    @testset "knots" begin
        @test knots(BSplineBasis(3, 0:5)) == [0, 0, 0, 1, 2, 3, 4, 5, 5, 5]
        @test knots(BSplineBasis(2, [1.5, 4.0])) == [1.5, 1.5, 4.0, 4.0]
    end

    @testset "order" begin
        @test order(BSplineBasis(5, 0:5)) == 5
        @test order(BSplineBasis(3, 0:1//2:10)) == 3
        @test order(BSplineBasis(5, 0:5)[1]) == 5
        @test order(BSplineBasis(3, 0:1//2:10)[2]) == 3
    end

    @testset "Printing" begin
        @test summary(BSplineBasis(5, 0:5)) == "9-element $(BSplineBasis{UnitRange{Int}})"
        @test summary(BSplineBasis(8, [1.0:0.1:3.0;])) == "27-element $(BSplineBasis{Vector{Float64}})"
        @test summary(Spline(BSplineBasis(5, 0:5), zeros(9))) == string(Spline{BSplineBasis{UnitRange{Int}},Vector{Float64}})
        @test summary(Spline(BSplineBasis(8, [1.0:0.1:3.0;]), 1:27)) == string(Spline{BSplineBasis{Vector{Float64}},UnitRange{Int}})
        @test summary(BSpline(BSplineBasis(5, 0:5), 1)) == "BSpline{$(BSplineBasis{UnitRange{Int}})}"
        @test summary(BSpline(BSplineBasis(8, [1.0:0.1:3.0;]), 1)) == "BSpline{$(BSplineBasis{Vector{Float64}})}"
        @test repr(MIME"text/plain"(), BSplineBasis(5, [1,π])) ==
            """
            5-element $(BSplineBasis{Vector{Float64}}):
             order: 5
             breakpoints: [1.0, 3.14159]"""
        @test repr(MIME"text/plain"(), Spline(BSplineBasis(5, [1,π]), sqrt.(0:4))) ==
            """
            $(Spline{BSplineBasis{Vector{Float64}},Vector{Float64}}):
             basis: 5-element $(BSplineBasis{Vector{Float64}}):
              order: 5
              breakpoints: [1.0, 3.14159]
             coeffs: [0.0, 1.0, 1.41421, 1.73205, 2.0]"""
        @test repr(MIME"text/plain"(), BSpline(BSplineBasis(5, [1,π]), 2)) ==
            """
            BSpline{$(BSplineBasis{Vector{Float64}})}:
             basis: 5-element $(BSplineBasis{Vector{Float64}}):
              order: 5
              breakpoints: [1.0, 3.14159]
             index: 2 (knots: [1.0, 1.0, 1.0, 1.0, 3.14159, 3.14159])"""
    end
end

@testset "Spline" begin
    @testset "Construction" begin
        @test_throws DimensionMismatch Spline(BSplineBasis(3, 0:5), zeros(8))
        @test_throws DimensionMismatch Spline(BSplineBasis(3, 0:5), OffsetArray{Int}(undef, 2:8))
        @test_throws DomainError BSpline(BSplineBasis(3, 0:5), 0)
        @test_throws DomainError BSpline(BSplineBasis(3, 0:5), 8)
        @test BSplineBasis(3, 0:5)[4] === BSpline(BSplineBasis(3, 0:5), 4)
    end

    @testset "basis" begin
        @test (b = BSplineBasis(4, [-1, 0, 1.5, 5, 10]); basis(Spline(b, ones(7))) === b)
        @test (b = BSplineBasis(3, 0:5); basis(Spline(b, 1:7)) === b)
        @test (b = BSplineBasis(5, 0:5); basis(Spline(b, zeros(9))) === b)
        @test (b = BSplineBasis(5, 0:5); basis(BSpline(b, 5)) === b)
    end

    @testset "coeffs" begin
        @test (c = ones(7); coeffs(Spline(BSplineBasis(4, [-1, 0, 1.5, 5, 10]), c)) === c)
        @test (c = 1:7; coeffs(Spline(BSplineBasis(3, 0:5), c)) === c)
        @test coeffs(BSpline(BSplineBasis(5, 0:10), 5)) == [0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    end

    @testset "order" begin
        @test order(Spline(BSplineBasis(5, 0:5), 1:9)) == 5
        @test order(Spline(BSplineBasis(3, 0:1//2:10), zeros(22))) == 3
        @test order(BSpline(BSplineBasis(4, [1:10;]), 3)) == 4
    end

    @testset "Equality and hashing" begin
        splines = [Spline(BSplineBasis(3, 0:5), Float64[1:7;]),
                   Spline(BSplineBasis(3, Float64[0:5;]), 1:7),
                   BSplineBasis(4, 0:4)[1],
                   BSpline(BSplineBasis(4, 0:4), 1),
                   Spline(BSplineBasis(4, [0,1,2,3,4]), [1,0,0,0,0,0,0]),
                   BSplineBasis(4, 0:5)[1],
                   BSplineBasis(4, Float64[0:5;])[2]]
        for s1 = splines, s2 = splines
            @test (s1 == s2) == (basis(s1) == basis(s2) && coeffs(s1) == coeffs(s2))
            @test (s1 == s2) == (hash(s1) == hash(s2))
        end
    end

    @testset "Broadcasting" begin
        # Splines act as scalars for broadcasting
        @test size(Spline(BSplineBasis(3, 0:5), 1:7) .== BSplineBasis(3, 0:5)[1]) == ()
    end

    @testset "Arithmetic" begin
        @testset "Unary +/-" begin
            @test +(BSplineBasis(3, 0:5)[4]) == Spline(BSplineBasis(3, 0:5), [0,0,0,1,0,0,0])
            @test +(Spline(BSplineBasis(4, [1:10;]), 1:12)) == Spline(BSplineBasis(4, 1:10), 1:12)
            @test +(Spline(BSplineBasis(4, [1:10;]), [1:12;])) == Spline(BSplineBasis(4, 1:10), 1:12)
            @test -(BSplineBasis(3, 0:5)[4]) == Spline(BSplineBasis(3, 0:5), [0,0,0,-1,0,0,0])
            @test -(Spline(BSplineBasis(4, [1:10;]), 1:12)) == Spline(BSplineBasis(4, 1:10), -1:-1:-12)
            @test -(Spline(BSplineBasis(4, [1:10;]), [1:12;])) == Spline(BSplineBasis(4, 1:10), -1:-1:-12)
        end

        @testset "Binary +/-" begin
            # Addition
            @test (BSplineBasis(3, 0:5)[1]) + (BSplineBasis(3, 0:5)[4]) == Spline(BSplineBasis(3, 0:5), [1,0,0,1,0,0,0])
            @test (BSplineBasis(3, 0:5)[2]) + (BSplineBasis(3, [0:5;])[2]) == Spline(BSplineBasis(3, 0:5), [0,2,0,0,0,0,0])
            @test (BSplineBasis(4, [1:10;])[1]) + Spline(BSplineBasis(4, 1:10), ones(12)) ==
                Spline(BSplineBasis(4, 1:10), [2,1,1,1,1,1,1,1,1,1,1,1])
            @test Spline(BSplineBasis(4, [1:10;]), -11:0) + (BSplineBasis(4, [1:10;])[12]) ==
                Spline(BSplineBasis(4, 1:10), [-11,-10,-9,-8,-7,-6,-5,-4,-3,-2,-1,1])
            @test Spline(BSplineBasis(5, -6:2:6), sqrt.(1:10)) + Spline(BSplineBasis(5, -6:2:6), log.(1:10)) ==
                Spline(BSplineBasis(5, -6:2:6), [sqrt(i)+log(i) for i=1:10])
            @test Spline(BSplineBasis(5, -3:1:3), 1.0:10.0) + Spline(BSplineBasis(5, -3:3), -10.0:-1.0) ==
                Spline(BSplineBasis(5, -3:3), -9.0:2.0:9.0)
            # Subtraction
            @test (BSplineBasis(3, 0:5)[1]) - (BSplineBasis(3, 0:5)[4]) == Spline(BSplineBasis(3, 0:5), [1,0,0,-1,0,0,0])
            @test (BSplineBasis(3, 0:5)[2]) - (BSplineBasis(3, [0:5;])[2]) == Spline(BSplineBasis(3, 0:5), zeros(7))
            @test (BSplineBasis(4, [1:10;])[1]) - Spline(BSplineBasis(4, 1:10), ones(12)) ==
                Spline(BSplineBasis(4, 1:10), [0,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1])
            @test Spline(BSplineBasis(4, [1:10;]), -11:0) - (BSplineBasis(4, [1:10;])[12]) ==
                Spline(BSplineBasis(4, 1:10), [-11,-10,-9,-8,-7,-6,-5,-4,-3,-2,-1,-1])
            @test Spline(BSplineBasis(5, -6:2:6), sqrt.(1:10)) - Spline(BSplineBasis(5, -6:2:6), log.(1:10)) ==
                Spline(BSplineBasis(5, -6:2:6), [sqrt(i)-log(i) for i=1:10])
            @test Spline(BSplineBasis(5, -3:1:3), 1.0:10.0) - Spline(BSplineBasis(5, -3:3), -10.0:-1.0) ==
                Spline(BSplineBasis(5, -3:3), fill(11.0, 10))
            # Different bases
            @test_throws ArgumentError (BSplineBasis(3, 0:5)[2]) + (BSplineBasis(4, 0:4)[2])
            @test_throws ArgumentError (BSplineBasis(3, 0:5)[2]) + Spline(BSplineBasis(4, 0:5), ones(8))
            @test_throws ArgumentError Spline(BSplineBasis(3, 0:5), zeros(7)) + Spline(BSplineBasis(4, 1:10), ones(12))
            @test_throws ArgumentError (BSplineBasis(3, 0:5)[2]) - (BSplineBasis(4, 0:4)[2])
            @test_throws ArgumentError (BSplineBasis(3, 0:5)[2]) - Spline(BSplineBasis(4, 0:5), ones(8))
            @test_throws ArgumentError Spline(BSplineBasis(3, 0:5), zeros(7)) - Spline(BSplineBasis(4, 1:10), ones(12))
        end

        @testset "Multiplication/division" begin
            # Base.*
            @test (BSplineBasis(3, 0:5)[5]) * 2 == Spline(BSplineBasis(3, 0:5), [0,0,0,0,2,0,0])
            @test Spline(BSplineBasis(4, [1:10;]), [sqrt(i) for i=1:12]) * 0.1 ==
                Spline(BSplineBasis(4, [1:10;]), [sqrt(i)*0.1 for i=1:12])
            @test 3.0 * (BSplineBasis(3, 0:5)[5])  == Spline(BSplineBasis(3, 0:5), [0,0,0,0,3,0,0])
            @test (5//2) * Spline(BSplineBasis(4, [1:10;]), [sqrt(i) for i=1:12]) ==
                Spline(BSplineBasis(4, [1:10;]), [(5//2)*sqrt(i) for i=1:12])
            # Base./
            @test (BSplineBasis(3, 0:5)[5]) / 2 == Spline(BSplineBasis(3, 0:5), [0,0,0,0,.5,0,0])
            @test Spline(BSplineBasis(4, [1:10;]), [sqrt(i) for i=1:12]) / 3.0 ==
                Spline(BSplineBasis(4, [1:10;]), [sqrt(i)/3.0 for i=1:12])
            # Base.\
            @test 10 \ (BSplineBasis(3, 0:5)[5])  == Spline(BSplineBasis(3, 0:5), [0,0,0,0,1/10,0,0])
            @test (5//2) \ Spline(BSplineBasis(4, [1:10;]), [sqrt(i) for i=1:12]) ==
                Spline(BSplineBasis(4, [1:10;]), [(5//2)\sqrt(i) for i=1:12])
            # LinearAlgebra.lmul!
            @test begin
                spl = Spline(BSplineBasis(4, [1:10;]), [sqrt(i) for i=1:12])
                lmul!(5//2, spl) === spl && spl == Spline(BSplineBasis(4, [1:10;]), [(5//2)*sqrt(i) for i=1:12])
            end
            # LinearAlgebra.rmul!
            @test begin
                spl = Spline(BSplineBasis(3, 0:5), [0,0,0,0,1,0,0])
                rmul!(spl, 2) === spl && spl == Spline(BSplineBasis(3, 0:5), [0,0,0,0,2,0,0])
            end
            @static if VERSION ≥ v"1.2"
                # LinearAlgebra.ldiv!
                @test begin
                    spl = Spline(BSplineBasis(3, 0:5), [0,0,0,0,1//1,0,0])
                    ldiv!(10, spl) === spl && spl == Spline(BSplineBasis(3, 0:5), [0,0,0,0,1//10,0,0])
                end
                # LinearAlgebra.rdiv!
                @test begin
                    spl = Spline(BSplineBasis(4, [1:10;]), [sqrt(i) for i=1:12])
                    rdiv!(spl, 3.0) === spl && spl ==  Spline(BSplineBasis(4, [1:10;]), [sqrt(i)/3.0 for i=1:12])
                end
            end
        end
    end
end

@testset "Derivatives" begin
    # Derivative
    @test Derivative{1}() isa Derivative
    @test_throws DomainError Derivative{-1}()
    @test_throws ArgumentError Derivative{UInt(1)}()
    @test Derivative(2) === Derivative{2}()
    @test Derivative(Int8(3)) === Derivative{3}()
    @test Derivative(UInt(4)) === Derivative{4}()
    @test_throws DomainError Derivative(-1)
    @test_throws DomainError Derivative(Int8(-2))
    @test Derivative(2) == Derivative(2)
    @test Derivative(2) != Derivative(3)
    # AllDerivatives
    @test AllDerivatives{2}() isa AllDerivatives
    @test_throws DomainError AllDerivatives{0}()
    @test_throws DomainError AllDerivatives{-1}()
    @test_throws ArgumentError AllDerivatives{UInt(1)}()
    @test AllDerivatives(2) === AllDerivatives{2}()
    @test AllDerivatives(Int8(3)) === AllDerivatives{3}()
    @test AllDerivatives(UInt(4)) === AllDerivatives{4}()
    @test_throws DomainError AllDerivatives(0)
    @test_throws DomainError AllDerivatives(-1)
    @test_throws DomainError AllDerivatives(Int8(-2))
    @test_throws DomainError AllDerivatives(UInt(0))
    @test AllDerivatives(2) == AllDerivatives(2)
    @test AllDerivatives(2) != AllDerivatives(3)
    # Derivative/AllDerivatives act as scalars for broadcasting
    @test size(Derivative(1) .== Derivative(1)) == ()
    @test size(AllDerivatives(4) .== AllDerivatives(4)) == ()
end

include("bases.jl")
include("bsplines.jl")
include("splinevalue.jl")

@testset "SplineFunctionWrapper" begin
    function test_splinefunctionwrapper(spline)
        a,b = support(spline)
        f  = Function(spline)
        f⁰ = Function(spline, Derivative(0))
        f¹ = Function(spline, Derivative(1))
        f² = Function(spline, Derivative(2))
        g  = Function(spline, false)
        g⁰ = Function(spline, Derivative(0), false)
        g¹ = Function(spline, Derivative(1), false)
        g² = Function(spline, Derivative(2), false)
        xs = range(a, stop=b, length=10)
        @test all(f.(xs) .== spline.(xs))
        @test all(f⁰.(xs) .== spline.(xs))
        @test all(f¹.(xs) .== spline.(xs, Derivative(1)))
        @test all(f².(xs) .== spline.(xs, Derivative(2)))
        @test all(g.(xs) .== spline.(xs))
        @test all(g⁰.(xs) .== spline.(xs))
        @test all(g¹.(xs) .== spline.(xs, Derivative(1)))
        @test all(g².(xs) .== spline.(xs, Derivative(2)))
        for x = (a-1, b+1)
            for fun = (f, f⁰, f¹, f²)
                @test isnan(fun(x))
            end
            for fun = (g, g⁰, g¹, g²)
                @test iszero(fun(x))
            end
        end
    end
    for spline = [BSpline(BSplineBasis(4, 0:2:20), 8),
                  Spline(BSplineBasis(3, 0//1:5//1), 1:7),
                  Spline(BSplineBasis(5, [0:0.1:1;]), ones(14))]
        test_splinefunctionwrapper(spline)
    end
end

@testset "basismatrix" begin
    basis = BSplineBasis(3, 0:5)
    xvalues = [0.3, 1.5, 3.2]
    @test basismatrix(basis, xvalues) == [basis[i](x) for x=xvalues, i=eachindex(basis)]
    @test basismatrix(basis, xvalues, indices=2:5) == [basis[i](x) for x=xvalues, i=2:5]
    @test_throws ArgumentError basismatrix(basis, xvalues, indices=0:5)
end

@testset "basismatrix!" begin
    basis = BSplineBasis(3, 0:5)
    xvalues = [0.3, 1.5, 3.2]
    @test begin
        dest = Matrix{Float64}(undef, 3, length(basis))
        basismatrix!(dest, basis, xvalues)
        dest == [basis[i](x) for x=xvalues, i=eachindex(basis)]
    end
    @test begin
        dest = Matrix{Float64}(undef, 3, 4)
        basismatrix!(dest, basis, xvalues, indices=2:5)
        dest == [basis[i](x) for x=xvalues, i=2:5]
    end
    @test_throws ArgumentError basismatrix!(Array{Float64}(undef, 3, 6), basis, xvalues, indices=0:5)
    @test_throws ArgumentError basismatrix!(OffsetArray(Array{Float64}(undef, 3, length(basis)),1,1), basis, xvalues)
    @test_throws ArgumentError basismatrix!(OffsetArray(Array{Float64}(undef, 3, 4),1,1), basis, xvalues, indices=2:5)
    @test_throws DimensionMismatch basismatrix!(Array{Float64}(undef, 3, 4), basis, xvalues)
    @test_throws DimensionMismatch basismatrix!(Array{Float64}(undef, 2, 4), basis, xvalues, indices=2:5)
end

@testset "approximate" begin
    basis = BSplineBasis(3, 0:5)
    knotavgs = [0.0, 0.5, 1.5, 2.5, 3.5, 4.5, 5.0]
    @test begin
        spline = @inferred(approximate(sqrt, basis))
        all(spline(x) ≈ sqrt(x) for x=knotavgs)
    end
    @test @inferred(approximate(sqrt, basis, indices=:)) == approximate(sqrt, basis)
    @test @inferred(approximate(sqrt, basis, indices=1:7)) == approximate(sqrt, basis)
    @test_throws ArgumentError approximate(sqrt, basis, indices=0:7)

    @test begin
        spline = approximate(sqrt, basis, indices=2:7)
        all(spline(x) ≈ sqrt(x) for x=knotavgs) && iszero(coeffs(spline)[1])
    end

    basis = BSplineBasis(5, [-5.0, -3.5, -2.2, 1.0, 2.0, 3.6])
    knotavgs = [-5.0, -4.625, -3.925, -2.425, -0.675, 1.1, 2.55, 3.2, 3.6]
    # basis = BSplineBasis(5, [-5, -7//2, -11//5, 1, 2, 18//5])
    # knotavgs = [-5, -37//8, -157//40, -97//40, -27//40, 11//10, 51//20, 16//5, 18//5]
    @test begin
        spline = @inferred(approximate(sin, basis))
        all(spline(x) ≈ sin(x) for x = knotavgs)
    end
    @test @inferred(approximate(sin, basis, indices=:)) == approximate(sin, basis)
    @test @inferred(approximate(sin, basis, indices=1:9)) == approximate(sin, basis)
    @test_throws ArgumentError approximate(sin, basis, indices=2:10)

    @test begin
        spline = approximate(one, BSplineBasis(3, 0:5))
        all(spline(x) ≈ 1 for x=range(0, stop=5, length=100))
    end
    @test begin
        basis = BSplineBasis(9, 0:0.1:1)
        spline = approximate(sinpi, basis, indices=2:length(basis)-1)
        all(spline(x) ≈ sinpi(x) for x=range(0, stop=1, length=100))
    end
end

@testset "interpolate" begin
    basis = BSplineBasis(3, 0:5)
    xvalues = [0.0, 0.5, 1.5, 2.5, 3.5, 4.5, 5.0]
    yvalues = sqrt.(xvalues)
    @test begin
        spline = @inferred(interpolate(basis, xvalues, yvalues))
        all(spline(xvalues[i]) ≈ yvalues[i] for i=eachindex(xvalues))
    end
    @test @inferred(interpolate(basis, xvalues, yvalues, indices=:)) == interpolate(basis, xvalues, yvalues)
    @test @inferred(interpolate(basis, xvalues, yvalues, indices=1:7)) == interpolate(basis, xvalues, yvalues)
    @test_throws ArgumentError interpolate(basis, xvalues, yvalues, indices=0:7)
    @test_throws DimensionMismatch interpolate(basis, xvalues, [yvalues; 0])
    @test_throws DimensionMismatch interpolate(basis, [xvalues; 0], yvalues, indices=:)
    @test_throws DimensionMismatch interpolate(basis, [xvalues; 0], yvalues, indices=1:7)

    @test begin
        spline = interpolate(basis, xvalues, yvalues, indices=2:7)
        all(spline(xvalues[i]) ≈ yvalues[i] for i=eachindex(xvalues)) && iszero(coeffs(spline)[1])
    end

    xvalues = xvalues[2:end]
    yvalues = yvalues[2:end]
    @test begin
        spline = interpolate(basis, xvalues, yvalues, indices=2:7)
        all(spline(xvalues[i]) ≈ yvalues[i] for i=eachindex(xvalues)) && iszero(coeffs(spline)[1])
    end

    basis = BSplineBasis(5, [-5.0, -3.5, -2.2, 1.0, 2.0, 3.6])
    xvalues = [-5.0, -4.625, -3.925, -2.425, -0.675, 1.1, 2.55, 3.2, 3.6]
    # basis = BSplineBasis(5, [-5, -7//2, -11//5, 1, 2, 18//5])
    # xvalues = [-5, -37//8, -157//40, -97//40, -27//40, 11//10, 51//20, 16//5, 18//5]
    yvalues = sin.(xvalues)
    @test begin
        spline = @inferred(interpolate(basis, xvalues, yvalues))
        all(spline(xvalues[i]) ≈ yvalues[i] for i=eachindex(xvalues))
    end
    @test @inferred(interpolate(basis, xvalues, yvalues, indices=:)) == interpolate(basis, xvalues, yvalues)
    @test @inferred(interpolate(basis, xvalues, yvalues, indices=1:9)) == interpolate(basis, xvalues, yvalues)
    @test_throws ArgumentError interpolate(basis, xvalues, yvalues, indices=2:10)
    @test_throws DimensionMismatch interpolate(basis, xvalues, [yvalues; 0])
    @test_throws DimensionMismatch interpolate(basis, [xvalues; 0], yvalues, indices=:)
    @test_throws DimensionMismatch interpolate(basis, [xvalues; 0], yvalues, indices=1:9)

    yvalues[end] = 0
    @test begin
        spline = interpolate(basis, xvalues, yvalues, indices=1:8)
        all(spline(xvalues[i]) ≈ yvalues[i] for i=eachindex(xvalues)) && iszero(coeffs(spline)[end])
    end

    xvalues = xvalues[1:end-1]
    yvalues = yvalues[1:end-1]
    @test begin
        spline = interpolate(basis, xvalues, yvalues, indices=1:8)
        all(spline(xvalues[i]) ≈ yvalues[i] for i=eachindex(xvalues)) && iszero(coeffs(spline)[end])
    end
end

@testset "averagebasis" begin
    ≈ₜₑₛₜ(x::BSplineBasis, y::BSplineBasis; kwargs...) =
        order(x) == order(y) && isapprox(breakpoints(x), breakpoints(y); kwargs...)

    @test averagebasis(3, [0.0, 0.5, 1.5, 2.5, 3.5, 4.5, 5.0]) ≈ₜₑₛₜ BSplineBasis(3, 0:5)
    @test averagebasis(3, [0//1, 1//2, 3//2, 5//2, 7//2, 9//2, 5//1]) == BSplineBasis(3, 0:5)
    @test averagebasis(4, [0.0, 1/3, 1.0, 2.0, 3.0, 4.0, 14/3, 5.0]) ≈ₜₑₛₜ
        BSplineBasis(4, [0.0, 10/9, 2.0, 3.0, 35/9, 5])
    @test averagebasis(4, [0//1, 1//3, 1//1, 2//1, 3//1, 4//1, 14//3, 5//1]) ==
        BSplineBasis(4, [0//1, 10//9, 2//1, 3//1, 35//9, 5//1])
    @test averagebasis(5, 1:10) ≈ₜₑₛₜ BSplineBasis(5, [1.0, 3.5, 4.5, 5.5, 6.5, 7.5, 10.0])
    @test averagebasis(5, 1//1:10//1) == BSplineBasis(5, [1//1, 7//2, 9//2, 11//2, 13//2, 15//2, 10//1])
    @test averagebasis(6, 1:1:10) ≈ₜₑₛₜ BSplineBasis(6, [1, 4, 5, 6, 7, 10])
    @test averagebasis(6, 1//1:1//1:10//1) == BSplineBasis(6, [1, 4, 5, 6, 7, 10])
    @test averagebasis(7, 1:7) ≈ₜₑₛₜ BSplineBasis(7, [1, 7])
    @test averagebasis(7, 1//1:7//1) == BSplineBasis(7, [1, 7])
    @test_throws ArgumentError averagebasis(8, 1:7) # too few datapoints
end

@testset "knotaverages" begin
    @test knotaverages(BSplineBasis(3, 0//1:5//1))              == [0, 1//2, 3//2, 5//2, 7//2, 9//2, 5]
    @test knotaverages(BSplineBasis(3, 0//1:5//1), indices=:)   == [0, 1//2, 3//2, 5//2, 7//2, 9//2, 5]
    @test knotaverages(BSplineBasis(3, 0//1:5//1), indices=2:6) == [1//2, 3//2, 5//2, 7//2, 9//2]
    @test_throws ArgumentError knotaverages(BSplineBasis(3, 0//1:5//1), indices=0:3)
    @test_throws ArgumentError knotaverages(BSplineBasis(3, 0//1:5//1), indices=1:8)
    @test knotaverages(BSplineBasis(4, 0:1:5))              ≈ [0, 1/3, 1, 2, 3, 4, 14/3, 5]
    @test knotaverages(BSplineBasis(4, 0:1:5), indices=:)   ≈ [0, 1/3, 1, 2, 3, 4, 14/3, 5]
    @test knotaverages(BSplineBasis(4, 0:1:5), indices=2:8) ≈ [1/3, 1, 2, 3, 4, 14/3, 5]
    @test knotaverages(BSplineBasis(4, 0.0:5.0))              ≈ [0, 1/3, 1, 2, 3, 4, 14/3, 5]
    @test knotaverages(BSplineBasis(4, 0.0:5.0), indices=:)   ≈ [0, 1/3, 1, 2, 3, 4, 14/3, 5]
    @test knotaverages(BSplineBasis(4, 0.0:5.0), indices=2:8) ≈ [1/3, 1, 2, 3, 4, 14/3, 5]
    @test knotaverages(BSplineBasis(5, [-5.0, -3.5, -2.2, 1.0, 2.0, 3.6])) ≈
        [-5.0, -4.625, -3.925, -2.425, -0.675, 1.1, 2.55, 3.2, 3.6]
    @test knotaverages(BSplineBasis(5, [-5.0, -3.5, -2.2, 1.0, 2.0, 3.6]), indices=:) ≈
        [-5.0, -4.625, -3.925, -2.425, -0.675, 1.1, 2.55, 3.2, 3.6]
    @test knotaverages(BSplineBasis(5, [-5.0, -3.5, -2.2, 1.0, 2.0, 3.6]), indices=1:8) ≈
        [-5.0, -4.625, -3.925, -2.425, -0.675, 1.1, 2.55, 3.2]
    @test knotaverages(BSplineBasis(5, [-5//1, -7//2, -11//5, 1//1, 2//1, 18//5])) ==
        [-5//1, -37//8, -157//40, -97//40, -27//40, 11//10, 51//20, 16//5, 18//5]
    @test knotaverages(BSplineBasis(5, [-5//1, -7//2, -11//5, 1//1, 2//1, 18//5]), indices=:) ==
        [-5//1, -37//8, -157//40, -97//40, -27//40, 11//10, 51//20, 16//5, 18//5]
    @test knotaverages(BSplineBasis(5, [-5//1, -7//2, -11//5, 1//1, 2//1, 18//5]), indices=1:8) ==
        [-5//1, -37//8, -157//40, -97//40, -27//40, 11//10, 51//20, 16//5]
end

@testset "knotaverages!" begin
    @test begin
        dest = Vector{Float64}(undef, 7)
        knotaverages!(dest, BSplineBasis(3, 0//1:5//1))
        dest ≈ [0, 1/2, 3/2, 5/2, 7/2, 9/2, 5]
    end
    @test begin
        dest = Vector{Float64}(undef, 7)
        knotaverages!(dest, BSplineBasis(3, 0//1:5//1))
        dest ≈ [0, 1/2, 3/2, 5/2, 7/2, 9/2, 5]
    end
    @test begin
        dest = Vector{Rational{Int}}(undef, 5)
        knotaverages!(dest, BSplineBasis(3, 0//1:5//1), indices=2:6)
        dest == [1//2, 3//2, 5//2, 7//2, 9//2]
    end
    @test_throws ArgumentError knotaverages!(zeros(4), BSplineBasis(3, 0//1:5//1), indices=0:3)
    @test_throws ArgumentError knotaverages!(zeros(8), BSplineBasis(3, 0//1:5//1), indices=1:8)
    @test_throws ArgumentError knotaverages!(OffsetArray{Float64}(undef, 2:8), BSplineBasis(3, 0//1:5//1))
    @test_throws ArgumentError knotaverages!(OffsetArray{Float64}(undef, 2:8), BSplineBasis(3, 0//1:5//1), indices=:)
    @test_throws ArgumentError knotaverages!(OffsetArray{Float64}(undef, 0:5), BSplineBasis(3, 0//1:5//1), indices=2:7)
    @test_throws DimensionMismatch knotaverages!(zeros(6), BSplineBasis(3, 0//1:5//1))
    @test_throws DimensionMismatch knotaverages!(zeros(6), BSplineBasis(3, 0//1:5//1), indices=:)
    @test_throws DimensionMismatch knotaverages!(zeros(7), BSplineBasis(3, 0//1:5//1), indices=2:7)
    @test begin
        dest = Vector{Float32}(undef, 8)
        knotaverages!(dest, BSplineBasis(4, 0:1:5))
        dest ≈ [0, 1/3, 1, 2, 3, 4, 14/3, 5]
    end
    @test begin
        dest = Vector{Float64}(undef, 8)
        knotaverages!(dest, BSplineBasis(4, 0:1:5), indices=:)
        dest ≈ [0, 1/3, 1, 2, 3, 4, 14/3, 5]
    end
    @test begin
        dest = Vector{Rational{Int}}(undef, 5)
        knotaverages!(dest, BSplineBasis(4, 0:1:5), indices=2:6)
        dest == [1//3, 1, 2, 3, 4]
    end
    @test begin
        dest = Vector{Float32}(undef, 9)
        knotaverages!(dest, BSplineBasis(5, [-5.0, -3.5, -2.2, 1.0, 2.0, 3.6]))
        dest ≈ [-5.0, -4.625, -3.925, -2.425, -0.675, 1.1, 2.55, 3.2, 3.6]
    end
    @test begin
        dest = Vector{Float64}(undef, 9)
        knotaverages!(dest, BSplineBasis(5, [-5.0, -3.5, -2.2, 1.0, 2.0, 3.6]), indices=:)
        dest ≈ [-5.0, -4.625, -3.925, -2.425, -0.675, 1.1, 2.55, 3.2, 3.6]
    end
    @test begin
        dest = Vector{BigFloat}(undef, 8)
        knotaverages!(dest, BSplineBasis(5, [-5.0, -3.5, -2.2, 1.0, 2.0, 3.6]), indices=1:8)
        dest ≈ [-5.0, -4.625, -3.925, -2.425, -0.675, 1.1, 2.55, 3.2]
    end
end
