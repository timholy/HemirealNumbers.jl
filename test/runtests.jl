using HemirealNumbers
using Test

@testset "HemirealNumbers" begin
    @testset "ambiguities" begin
        @test isempty(detect_ambiguities(HemirealNumbers))
    end

    p1 = @inferred(PureHemi(1, 2))
    p2 = PureHemi(1, 2)
    p3 = @inferred(PureHemi(1.0, 2))
    @test p1 == p2
    @test p1 == p3

    @test @inferred(convert(PureHemi{Float64}, p1)) === PureHemi(1.0, 2.0)
    @test @inferred(convert(PureHemi{Float64}, 0)) == PureHemi(0.0, 0.0)
    @test_throws DomainError convert(PureHemi{Float64}, 1)
    @test convert(Float64, PureHemi(0, 0)) === 0.0
    @test_throws DomainError convert(Float64, PureHemi(1, 0))

    @test @inferred(-p1) == PureHemi(-1, -2)
    @test @inferred(+p1) == p1

    @test @inferred(PureHemi(1, 2) + PureHemi(5, 3)) === PureHemi(6, 5)
    @test @inferred((μ + 2ν) + (5μ + 3ν)) === 6μ + 5ν
    @test PureHemi(1, 2) - PureHemi(5, 3) == PureHemi(-4, -1)
    @test PureHemi(0.0, 0.0) == PureHemi(-0.0, -0.0)
    @test PureHemi(NaN, 1.0) != PureHemi(NaN, 1.0)
    @test isequal(PureHemi(NaN, 1.0), PureHemi(NaN, 1.0))
    @test μ * μ == 0
    @test ν * ν == 0
    @test ν * μ == 1
    @test μ * ν == 1
    @test 2 * PureHemi(3, 4) == PureHemi(6, 8)
    @test PureHemi(3, 4) * (-1) == PureHemi(-3, -4)
    @test false * PureHemi(3, 4) == 0
    @test PureHemi(3, 4) * true == 3μ + 4ν
    @test @inferred((3μ + 4ν) / 2) === 1.5μ + 2ν
    @test 2 \ (3μ + 4ν) === 1.5μ + 2ν
    @test 18 / PureHemi(3, 3) == PureHemi(3, 3)
    @test PureHemi(3, 3) \ 7 === PureHemi(7 / 6, 7 / 6)
    @test @inferred((4μ + 2ν) / (2μ + ν)) === 2.0
    @test_throws ArgumentError (4μ + 2ν) / (2μ + 3ν)
    @test @inferred((4μ + 2ν) \ (2μ + ν)) === 0.5
    @test PureHemi(3, 4)^2 == 24
    @test PureHemi(3, 4)^(1 // 4) == 24^(1 // 8)
    @test PureHemi(3, 4)^2.0 == 24.0
    @test real(PureHemi(0.2, 0.3)) == 0.0
    @test real(PureHemi{Bool}) == Bool
    @test mu(PureHemi(0.2, 0.3)) == 0.2
    @test nu(PureHemi(0.2, 0.3)) == 0.3
    @test zero(PureHemi(0.2, 0.3)) == PureHemi(0, 0)
    @test zero(PureHemi{Bool}) == PureHemi(false, false)
    @test conj(3μ + 4ν) == 3μ + 4ν
    @test isfinite(3μ + 4ν)
    @test !isfinite(Inf * μ + 3.2 * ν)
    @test PureHemi(1, 2) ≈ PureHemi(nextfloat(1.0), nextfloat(2.0))

    @testset "PureHemi predicates" begin
        @test iszero(PureHemi(0, 0))
        @test !iszero(PureHemi(1, 0))
        @test isnan(PureHemi(NaN, 0.0))
        @test isnan(PureHemi(0.0, NaN))
        @test !isnan(PureHemi(1.0, 2.0))
        @test isinf(PureHemi(Inf, 0.0))
        @test isinf(PureHemi(0.0, Inf))
        @test !isinf(PureHemi(1.0, 2.0))
        @test isreal(PureHemi(0, 0))
        @test !isreal(PureHemi(1, 0))
    end

    @testset "PureHemi float/widen/big/bswap" begin
        @test float(PureHemi(1, 2)) === PureHemi(1.0, 2.0)
        @test float(PureHemi(1.0, 2.0)) === PureHemi(1.0, 2.0)
        @test float(PureHemi{Int}) == PureHemi{Float64}
        @test float(PureHemi{Float32}) == PureHemi{Float32}
        @test widen(PureHemi{Float32}) == PureHemi{Float64}
        @test widen(PureHemi{Int32}) == PureHemi{Int64}
        @test big(PureHemi{Int}) == PureHemi{BigInt}
        @test big(PureHemi(1, 2)) == PureHemi(BigInt(1), BigInt(2))
        @test bswap(PureHemi(1, 2)) == PureHemi(bswap(1), bswap(2))
    end

    @testset "PureHemi hash" begin
        @test hash(PureHemi(0, 0)) == hash(0)
        @test hash(PureHemi(0.0, 0.0)) == hash(0)
        @test hash(PureHemi(1, 2)) == hash(PureHemi(1, 2))
        @test hash(PureHemi(1, 2)) != hash(PureHemi(2, 1))
    end

    @test @inferred(Hemireal(1, 2, 3)) === 1 + 2μ + 3ν
    @test @inferred(Hemireal(1.0, 2, 3)) === 1.0 + 2.0μ + 3.0ν
    @test @inferred(Hemireal(2μ + 3ν)) === Hemireal(0, 2, 3)
    @test @inferred(Hemireal(7.0, 2μ + 3ν)) === Hemireal(7.0, 2.0, 3.0)
    @test @inferred(Hemireal(7)) === Hemireal(7, 0, 0)

    @test convert(Hemireal{Float64}, 1) === Hemireal(1.0, 0.0, 0.0)
    @test convert(Hemireal{Float64}, 2μ + 3ν) === Hemireal(0.0, 2.0, 3.0)
    @test convert(Hemireal{Float64}, Hemireal(1, 2, 3)) === Hemireal(1.0, 2.0, 3.0)
    @test convert(Float64, Hemireal(1, 0, 0)) === 1.0
    @test_throws DomainError convert(Float64, Hemireal(1, 2, 3))
    @test convert(PureHemi{Float64}, Hemireal(0, 2, 3)) === PureHemi(2.0, 3.0)
    @test_throws DomainError convert(PureHemi{Float64}, Hemireal(1, 2, 3))

    @test @inferred(-(7 + 2μ)) === Hemireal(-7, -2, 0)
    @test @inferred(+(7 + 2μ)) === Hemireal(7, 2, 0)
    @test @inferred((1 + 2μ + 3ν) + (5.5 + 2.1μ + 3.2ν)) === Hemireal(6.5, 4.1, 6.2)
    @test (1 + 2μ + 3ν) - (5.5 + 2.1μ + 3.2ν) ≈ Hemireal(-4.5, -0.1, -0.2)
    @test (1 + 2μ + 3ν) * (5.5 + 2.1μ + 3.2ν) ≈ Hemireal(5.5 + 6.4 + 6.3, 13.1, 19.7)

    @test @inferred((1 + 2μ + 3ν) + PureHemi(2.1, 3.2)) === Hemireal(1, 4.1, 6.2)
    @test @inferred(PureHemi(7, 9) + (2 + μ - 6ν)) === Hemireal(2, 8, 3)
    @test (1 + 2μ + 3ν) - PureHemi(2.1, 3.2) ≈ Hemireal(1.0, -0.1, -0.2)
    @test @inferred(PureHemi(7, 9) - (2 + μ - 6ν)) === Hemireal(-2, 6, 15)
    @test (1 + 2μ + 3ν) * PureHemi(2.1, 3.2) ≈ Hemireal(6.4 + 6.3, 2.1, 3.2)
    @test @inferred(PureHemi(7, 9) * (2 + μ - 6ν)) === Hemireal(-33, 14, 18)

    @test @inferred((1 + 2μ + 3ν) + 5.5) === Hemireal(6.5, 2.0, 3.0)
    @test @inferred(3 + (2 + μ - 6ν)) === Hemireal(5, 1, -6)
    @test @inferred((1 + 2μ + 3ν) - 5.5) === Hemireal(-4.5, 2.0, 3.0)
    @test @inferred(3 - (2 + μ - 6ν)) === Hemireal(1, -1, 6)
    @test @inferred(false * (2 + μ - 6ν)) === Hemireal(0, 0, 0)
    @test @inferred((2 + μ - 6ν) * true) === Hemireal(2, 1, -6)
    @test @inferred(1.5 * (2 + μ - 6ν)) === Hemireal(3.0, 1.5, -9.0)
    @test @inferred((2 + μ - 6ν) * 1.5) === Hemireal(3.0, 1.5, -9.0)
    @test @inferred((2 + μ - 6ν) / 1.5) === Hemireal(2 / 1.5, 1 / 1.5, -6 / 1.5)
    @test @inferred(1.5 \ (2 + μ - 6ν)) === Hemireal(2 / 1.5, 1 / 1.5, -6 / 1.5)

    x = 3.2 + μ - ν
    @test real(x) === 3.2
    @test real(Hemireal{Float16}) == Float16
    @test mu(x) == 1
    @test nu(x) == -1
    @test mu([x]) == [1]
    @test nu([x]) == [-1]
    @test conj(x) === x
    @test zero(x) === Hemireal(0.0, 0.0, 0.0)
    @test zero(Hemireal{Float16}) === Hemireal{Float16}(0, 0, 0)

    @testset "Hemireal predicates" begin
        @test iszero(Hemireal(0, 0, 0))
        @test !iszero(Hemireal(1, 0, 0))
        @test !iszero(Hemireal(0, 1, 0))
        @test isnan(Hemireal(NaN, 0.0, 0.0))
        @test isnan(Hemireal(0.0, NaN, 0.0))
        @test isnan(Hemireal(0.0, 0.0, NaN))
        @test !isnan(Hemireal(1.0, 2.0, 3.0))
        @test isinf(Hemireal(Inf, 0.0, 0.0))
        @test isinf(Hemireal(0.0, Inf, 0.0))
        @test isinf(Hemireal(0.0, 0.0, Inf))
        @test !isinf(Hemireal(1.0, 2.0, 3.0))
        @test isreal(Hemireal(1.0, 0.0, 0.0))
        @test !isreal(Hemireal(1.0, 1.0, 0.0))
        @test isinteger(Hemireal(2, 0, 0))
        @test !isinteger(Hemireal(2, 1, 0))
        @test !isinteger(Hemireal(2.5, 0, 0))
        @test isone(Hemireal(1, 0, 0))
        @test !isone(Hemireal(2, 0, 0))
        @test !isone(Hemireal(1, 1, 0))
    end

    @testset "Hemireal one" begin
        @test one(Hemireal{Int}) === Hemireal(1, 0, 0)
        @test one(Hemireal(3.0, 1.0, 2.0)) === Hemireal(1.0, 0.0, 0.0)
        @test isone(one(Hemireal{Float64}))
    end

    @testset "Hemireal equality" begin
        @test Hemireal(1, 2, 3) == Hemireal(1, 2, 3)
        @test Hemireal(1, 0, 0) == 1
        @test 1 == Hemireal(1, 0, 0)
        @test Hemireal(1, 1, 0) != 1
        @test isequal(Hemireal(1, 2, 3), Hemireal(1, 2, 3))
        @test isequal(Hemireal(NaN, 0.0, 0.0), Hemireal(NaN, 0.0, 0.0))
    end

    @testset "Hemireal float/widen/big/bswap" begin
        @test float(Hemireal(1, 2, 3)) === Hemireal(1.0, 2.0, 3.0)
        @test float(Hemireal(1.0, 2.0, 3.0)) === Hemireal(1.0, 2.0, 3.0)
        @test float(Hemireal{Int}) == Hemireal{Float64}
        @test float(Hemireal{Float32}) == Hemireal{Float32}
        @test widen(Hemireal{Float32}) == Hemireal{Float64}
        @test widen(Hemireal{Int32}) == Hemireal{Int64}
        @test big(Hemireal{Int}) == Hemireal{BigInt}
        @test big(Hemireal(1, 2, 3)) == Hemireal(BigInt(1), BigInt(2), BigInt(3))
        @test bswap(Hemireal(1, 2, 3)) == Hemireal(bswap(1), bswap(2), bswap(3))
    end

    @testset "Hemireal hash" begin
        @test hash(Hemireal(0, 0, 0)) == hash(0)
        @test hash(Hemireal(3, 0, 0)) == hash(3)
        @test hash(Hemireal(3.0, 0.0, 0.0)) == hash(3)
        @test hash(Hemireal(1, 2, 3)) == hash(Hemireal(1, 2, 3))
        @test hash(Hemireal(1, 2, 3)) != hash(Hemireal(1, 3, 2))
        # Consistent with PureHemi when real part is zero
        @test hash(Hemireal(0, 2, 3)) == hash(PureHemi(2, 3))
    end

    @testset "binary I/O" begin
        buf = IOBuffer()
        ph = PureHemi(1.5, -2.5)
        write(buf, ph)
        seekstart(buf)
        @test read(buf, PureHemi{Float64}) === ph

        buf = IOBuffer()
        hr = Hemireal(3.0, 1.5, -2.5)
        write(buf, hr)
        seekstart(buf)
        @test read(buf, Hemireal{Float64}) === hr
    end

    a = [μ + 2ν, 5μ - ν]
    @test mu(a) == [1, 5]
    @test nu(a) == [2, -1]
    @test @inferred(a * a') == [4 9; 9 -10]
    @test isa(a * a', Matrix{Int})
    @test @inferred(a' * a) === -6
    A = [μ 0; 2μ + 3ν 5μ - 4ν]
    @test @inferred(A * a) == [2, -18]
    @test isa(A * a, Vector{Hemireal{Int}})
    z = zero(PureHemi{Int})
    A = [μ z; 2μ + 3ν 5μ - 4ν]
    @test @inferred(A * a) == [2, -18]
    @test isa(A * a, Vector{Int})
    @test @inferred(A * A) == [0 0; 10 -40]
    @test isa(A * A, Matrix{Int})
    @test @inferred(A * A') == [0 3;  3 -28]
    @test isa(A * A', Matrix{Int})
    @test @inferred(A * transpose(A)) == [0 3;  3 -28]
    @test isa(A * transpose(A), Matrix{Int})
    @test @inferred(A' * A) == [12 7; 7 -40]
    @test isa(A' * A, Matrix{Int})
    @test @inferred(transpose(A) * A) == [12 7; 7 -40]
    @test isa(transpose(A) * A, Matrix{Int})
    @test @inferred(A' * A') == [0 10; 0 -40]
    @test isa(A' * A', Matrix{Int})
    @test @inferred(transpose(A) * transpose(A)) == [0 10; 0 -40]
    @test isa(transpose(A) * transpose(A), Matrix{Int})

    @test @inferred(A * [1, 2]) == [μ, 12μ - 5ν]
    @test isa(A * [1, 2], Vector{PureHemi{Int}})

    @test μ * a == [2, -1]
    @test isa(μ * a, Vector{Int})
    @test 3 * a == [3μ + 6ν, 15μ - 3ν]
    @test isa(3 * a, Vector{PureHemi{Int}})

    # Promotion
    @test isa([μ, false], Vector{Hemireal{Bool}})
    @test isa([2μ, false], Vector{Hemireal{Int}})
    @test isa([μ, π], Vector{Hemireal{Float64}})
    @test isa([1 + μ, false], Vector{Hemireal{Int}})
    @test isa([1 + μ, π], Vector{Hemireal{Float64}})
    @test isa([1 + μ, 2.0 + ν], Vector{Hemireal{Float64}})

    @testset "show" begin
        @test sprint(show, 1μ + 2ν) == "1μ + 2ν"
        @test sprint(show, 1μ - 2ν) == "1μ - 2ν"
        @test sprint(show, 3 - 4.0μ + 2ν) == "3.0 - 4.0μ + 2.0ν"
        @test sprint(show, 3 - 4.0μ - 2ν) == "3.0 - 4.0μ - 2.0ν"
        # compact mode
        @test sprint(show, 1μ + 2ν; context=:compact => true) == "1μ+2ν"
        @test sprint(show, 1μ - 2ν; context=:compact => true) == "1μ-2ν"
        @test sprint(show, 3.0 + 4.0μ + 2.0ν; context=:compact => true) == "3.0+4.0μ+2.0ν"
    end
end
