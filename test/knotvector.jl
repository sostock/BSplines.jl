@testset "KnotVector" begin
    @testset "Construction" begin
        @test KnotVector([1.0, 1.5, 2.5], 3) == [1.0, 1.0, 1.0, 1.0, 1.5, 2.5, 2.5, 2.5, 2.5]
        @test KnotVector(1:5, 2, 3) == [1, 1, 1, 2, 3, 4, 5, 5, 5, 5]
        @test KnotVector(big(-1):big(1), 2) == [-1, -1, -1, 0, 1, 1, 1]
        @test KnotVector(1//2:5//2, 2, 0) == [1//2, 1//2, 1//2, 3//2, 5//2]
        @test KnotVector([true, false, true], 0) == [true, false, true]
        @test KnotVector(["test", "string"], 0, 1) == ["test", "string", "string"]
        @test KnotVector(KnotVector(UInt8(0):UInt8(3), 2, 1), 2) == [0, 0, 0, 0, 0, 1, 2, 3, 3, 3, 3]
        @test KnotVector([1.0, 1.5, 2.5], 3) isa AbstractVector{Float64}
        @test KnotVector(1:5, 2, 3) isa AbstractVector{Int}
        @test KnotVector(big(-1):big(1), 2) isa AbstractVector{BigInt}
        @test KnotVector(1//2:5//2, 2, 0) isa AbstractVector{Rational{Int}}
        @test KnotVector([true, false, true], 0) isa AbstractVector{Bool}
        @test KnotVector(["test", "string"], 0, 1) isa AbstractVector{String}
        @test KnotVector(KnotVector(UInt8(0):UInt8(3), 2, 1), 2) isa AbstractVector{UInt8}
        @test_throws ArgumentError KnotVector(1:0, 2)
        @test_throws ArgumentError KnotVector(Float64[], 2)
        @test_throws DomainError KnotVector(1:5, -1)
        @test_throws DomainError KnotVector(1:5, -1, 2)
        @test_throws DomainError KnotVector(1:5, 2, -1)
        @test_throws ArgumentError KnotVector(OffsetArray(1:5, 1), 2)
    end

    @testset "Dimensions" begin
        @test axes(KnotVector(1:5, 3)) == (1:11,)
        @test axes(KnotVector(1:5, 3), 1) == 1:11
        @test axes(KnotVector(1:5, 3), 2) == 1:1
        @test axes(KnotVector([1.0, 4.0], 2, 1)) == (1:5,)
        @test axes(KnotVector([1.0, 4.0], 2, 1), 1) == 1:5
        @test axes(KnotVector([1.0, 4.0], 2, 1), 2) == 1:1
        @test length(KnotVector(1:5, 3)) == 11
        @test length(KnotVector([1.0, 4.0], 2, 1)) == 5
        @test size(KnotVector(1:5, 3)) == (11,)
        @test size(KnotVector(1:5, 3), 1) == 11
        @test size(KnotVector(1:5, 3), 2) == 1
        @test size(KnotVector([1.0, 4.0], 2, 1)) == (5,)
        @test size(KnotVector([1.0, 4.0], 2, 1), 1) == 5
        @test size(KnotVector([1.0, 4.0], 2, 1), 2) == 1
    end

    @testset "Indices" begin
        @test firstindex(KnotVector(1:5, 3)) == 1
        @test firstindex(KnotVector(1:5, 3), 1) == 1
        @test firstindex(KnotVector(1:5, 3), 2) == 1
        @test firstindex(KnotVector([1.0, 4.0], 2, 1)) == 1
        @test firstindex(KnotVector([1.0, 4.0], 2, 1), 1) == 1
        @test firstindex(KnotVector([1.0, 4.0], 2, 1), 2) == 1
        @test lastindex(KnotVector(1:5, 3)) == 11
        @test lastindex(KnotVector(1:5, 3), 1) == 11
        @test lastindex(KnotVector(1:5, 3), 2) == 1
        @test lastindex(KnotVector([1.0, 4.0], 2, 1)) == 5
        @test lastindex(KnotVector([1.0, 4.0], 2, 1), 1) == 5
        @test lastindex(KnotVector([1.0, 4.0], 2, 1), 2) == 1
        @test eachindex(KnotVector(1:5, 3)) == 1:11
        @test eachindex(KnotVector([1.0, 4.0], 2, 1)) == 1:5
        @test keys(KnotVector(1:5, 3)) == 1:11
        @test keys(KnotVector([1.0, 4.0], 2, 1)) == 1:5
    end

    @testset "Indexing and iteration" begin
        @test collect(KnotVector(1:5, 3)) == [1, 1, 1, 1, 2, 3, 4, 5, 5, 5, 5]
        @test collect(KnotVector([1.0, 4.0], 2, 1)) == [1.0, 1.0, 1.0, 4.0, 4.0]
        @test first(KnotVector(1:5, 3)) == 1
        @test first(KnotVector([1.0, 4.0], 2, 1)) == 1.0
        @test last(KnotVector(1:5, 3)) == 5
        @test last(KnotVector([1.0, 4.0], 2, 1)) == 4.0
        @test KnotVector(1:5, 3)[:] == [1, 1, 1, 1, 2, 3, 4, 5, 5, 5, 5]
        @test KnotVector(1:5, 3)[1] == 1
        @test KnotVector(1:5, 3)[5] == 2
        @test KnotVector(1:5, 3)[end] == 5
        @test KnotVector(1:5, 3)[1:5] == [1, 1, 1, 1, 2]
        @test KnotVector(1:5, 3)[2,1] == 1
        @test KnotVector(1:5, 3)[1:5,1] == [1, 1, 1, 1, 2]
        @test_throws BoundsError KnotVector(1:5, 3)[0]
        @test_throws BoundsError KnotVector(1:5, 3)[12]
        @test_throws BoundsError KnotVector(1:5, 3)[1,2]
        @test_throws BoundsError KnotVector(1:5, 3)[1:5,2]
        @test_throws BoundsError KnotVector(1:5, 3)[1:20,1]
        @test KnotVector([1.0, 4.0], 2, 1)[:] == [1.0, 1.0, 1.0, 4.0, 4.0]
        @test KnotVector([1.0, 4.0], 2, 1)[1] == 1.0
        @test KnotVector([1.0, 4.0], 2, 1)[4] == 4.0
        @test KnotVector([1.0, 4.0], 2, 1)[end] == 4.0
        @test KnotVector([1.0, 4.0], 2, 1)[1:4] == [1.0, 1.0, 1.0, 4.0]
        @test KnotVector([1.0, 4.0], 2, 1)[4,1] == 4.0
        @test KnotVector([1.0, 4.0], 2, 1)[1:4,1] == [1.0, 1.0, 1.0, 4.0]
        @test_throws BoundsError KnotVector([1.0, 4.0], 2, 1)[0]
        @test_throws BoundsError KnotVector([1.0, 4.0], 2, 1)[6]
        @test_throws BoundsError KnotVector([1.0, 4.0], 2, 1)[1,2]
        @test_throws BoundsError KnotVector([1.0, 4.0], 2, 1)[1:4,2]
        @test_throws BoundsError KnotVector([1.0, 4.0], 2, 1)[1:10,1]
    end

    @testset "issorted" begin
        @test issorted(KnotVector(1:10, 5))
        @test issorted(KnotVector([1.0, 2.5, 5.0], 3, 5))
        @test !issorted(KnotVector(5:-1:1, 5))
        @test !issorted(KnotVector([1.0, 3.0, 2.0], 2, 1))
    end

    @testset "parent" begin
        @test (vec = [1.0, 1.5, 2.5]; parent(KnotVector(vec, 3)) === vec)
        @test (vec = 1:5; parent(KnotVector(vec, 2, 3)) === vec)
        @test (vec = big(-1):big(1); parent(KnotVector(vec, 2)) === vec)
        @test (vec = UInt8(0):UInt8(3); parent(KnotVector(KnotVector(vec, 2, 1), 2)) === vec)
    end

    @testset "unique/allunique" begin
        @test unique(KnotVector(1:3, 5)) == [1,2,3]
        @test unique(KnotVector([0], 2, 3)) == [0]
        @test unique(KnotVector([1,1,2,1], 0)) == [1,2]
        @test unique(KnotVector([1,1,2,1], 3)) == [1,2]
        @test allunique(KnotVector(1:5, 0))
        @test !allunique(KnotVector(1:5, 1))
        @test !allunique(KnotVector(1:5, 0, 1))
        @test !allunique(KnotVector(1:5, 2, 0))
        @test !allunique(KnotVector([1,1], 0))
    end
end
