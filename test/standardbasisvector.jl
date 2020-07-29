@testset "StandardBasisVector" begin
    @testset "Construction" begin
        @test StandardBasisVector(5, 2) === StandardBasisVector(Bool, 5, 2)
        @test StandardBasisVector(5//1, 2.0) === StandardBasisVector(Bool, 5, 2)
        @test_throws DomainError StandardBasisVector(5, 0)
        @test_throws DomainError StandardBasisVector(5, 6)
        @test_throws DomainError StandardBasisVector(0, 0)
        @test_throws DomainError StandardBasisVector(5, -2)
        @test_throws InexactError StandardBasisVector(5, 2.5)
        for T in (:Bool, :Int8, :Int, :UInt, :BigInt, :Float64, :BigFloat, :(Rational{Int}))
            @eval @test StandardBasisVector($T, 5, 2) isa AbstractVector{$T}
            @eval @test StandardBasisVector($T, 5//1, 2.0) isa AbstractVector{$T}
            @eval @test StandardBasisVector($T, 5, 2) == [0,1,0,0,0]
            @eval @test_throws DomainError StandardBasisVector($T, 5, 0)
            @eval @test_throws DomainError StandardBasisVector($T, 5, 6)
            @eval @test_throws DomainError StandardBasisVector($T, 5, -2)
            @eval @test_throws InexactError StandardBasisVector($T, 5, 2.5)
        end
    end

    @testset "copy" begin
        @test copy(StandardBasisVector(5, 2)) === StandardBasisVector(5, 2)
        @test copy(StandardBasisVector(Float64, 5, 2)) === StandardBasisVector(Float64, 5, 2)
    end

    @testset "Dimensions" begin
        @test axes(StandardBasisVector(10, 3)) == (1:10,)
        @test axes(StandardBasisVector(10, 3), 1) == 1:10
        @test axes(StandardBasisVector(10, 3), 2) == 1:1
        @test length(StandardBasisVector(10, 3)) == 10
        @test size(StandardBasisVector(10, 3)) == (10,)
        @test size(StandardBasisVector(10, 3), 1) == 10
        @test size(StandardBasisVector(10, 3), 2) == 1
    end

    @testset "Indices" begin
        @test firstindex(StandardBasisVector(8, 1)) == 1
        @test firstindex(StandardBasisVector(8, 2), 1) == 1
        @test firstindex(StandardBasisVector(8, 3), 2) == 1
        @test lastindex(StandardBasisVector(8, 4)) == 8
        @test lastindex(StandardBasisVector(8, 5), 1) == 8
        @test lastindex(StandardBasisVector(8, 6), 2) == 1
        @test eachindex(StandardBasisVector(8, 7)) == 1:8
        @test keys(StandardBasisVector(8, 8)) == 1:8
    end

    @testset "Indexing and iteration" begin
        for T in (:Bool, :Int8, :Int, :UInt, :BigInt, :Float64, :BigFloat, :(Rational{Int}))
            @eval @test collect(StandardBasisVector($T, 6, 3)) isa AbstractVector{$T}
            @eval @test collect(StandardBasisVector($T, 6, 3)) == [0,0,1,0,0,0]
            @eval @test @inferred(first(StandardBasisVector($T, 6, 3))) isa $T
            @eval @test @inferred(last(StandardBasisVector($T, 6, 3))) isa $T
            @eval @test @inferred(StandardBasisVector($T, 6, 5)[1]) isa $T
            @eval @test @inferred(StandardBasisVector($T, 6, 5)[2,1]) isa $T
            @eval @test first(StandardBasisVector($T, 6, 3)) == 0
            @eval @test first(StandardBasisVector($T, 6, 1)) == 1
            @eval @test last(StandardBasisVector($T, 6, 3)) == 0
            @eval @test last(StandardBasisVector($T, 6, 6)) == 1
            @eval @test StandardBasisVector($T, 6, 5)[:] === StandardBasisVector($T, 6, 5)
            @eval @test StandardBasisVector($T, 6, 5)[1] == 0
            @eval @test StandardBasisVector($T, 6, 5)[5] == 1
            @eval @test StandardBasisVector($T, 6, 5)[end] == 0
            @eval @test StandardBasisVector($T, 6, 5)[1:5] == [0,0,0,0,1]
            @eval @test StandardBasisVector($T, 6, 5)[2,1] == 0
            @eval @test StandardBasisVector($T, 6, 5)[1:5,1] == [0,0,0,0,1]
            @eval @test_throws BoundsError StandardBasisVector($T, 6, 5)[0]
            @eval @test_throws BoundsError StandardBasisVector($T, 6, 5)[7]
            @eval @test_throws BoundsError StandardBasisVector($T, 6, 5)[3,2]
            @eval @test_throws BoundsError StandardBasisVector($T, 6, 5)[1:5,2]
            @eval @test_throws BoundsError StandardBasisVector($T, 6, 5)[1:8,1]
        end
    end

    @testset "Equality" begin
        @test StandardBasisVector(8, 3) != StandardBasisVector(8, 4)
        @test StandardBasisVector(8, 3) == StandardBasisVector(8, 3)
        @test StandardBasisVector(UInt8, 8, 3) == StandardBasisVector(Int64, 8, 3)
        @test StandardBasisVector(UInt8, 7, 3) != StandardBasisVector(Int64, 8, 3)
    end
end
