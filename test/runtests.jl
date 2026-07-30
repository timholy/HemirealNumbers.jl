using HemiplexNumbers
using Test

@testset "HemiplexNumbers" begin
    @testset "ambiguities" begin
        @test isempty(detect_ambiguities(HemiplexNumbers))
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
    @test ν * μ == -0.5
    @test μ * ν == -0.5
    @test 2 * PureHemi(3, 4) == PureHemi(6, 8)
    @test PureHemi(3, 4) * (-1) == PureHemi(-3, -4)
    @test false * PureHemi(3, 4) == 0
    @test PureHemi(3, 4) * true == 3μ + 4ν
    @test @inferred((3μ + 4ν) / 2) === 1.5μ + 2ν
    @test 2 \ (3μ + 4ν) === 1.5μ + 2ν
    @test 18 / PureHemi(3, 3) == PureHemi(-6, -6)
    @test (18 / PureHemi(4, 2)) * PureHemi(4, 2) == 18
    @test PureHemi(3, 3) \ 7 === PureHemi(-7 / 3, -7 / 3)
    @test (PureHemi(4, 2) \ 18) * PureHemi(4, 2) == 18
    @test @inferred((4μ + 2ν) / (2μ + ν)) === 2.0
    @test_throws ArgumentError (4μ + 2ν) / (2μ + 3ν)
    @test @inferred((4μ + 2ν) \ (2μ + ν)) === 0.5
    @test PureHemi(3, 4)^2 == -12
    @test PureHemi(3, -4)^(1 // 4) == 12^(1 // 8)
    @test PureHemi(3, 4)^2.0 == -12.0
    @test real(PureHemi(0.2, 0.3)) == 0.0
    @test real(PureHemi{Bool}) == Bool
    @test mu(PureHemi(0.2, 0.3)) == 0.2
    @test nu(PureHemi(0.2, 0.3)) == 0.3
    @test zero(PureHemi(0.2, 0.3)) == PureHemi(0, 0)
    @test zero(PureHemi{Bool}) == PureHemi(false, false)
    @test conj(3μ + 4ν) == -(3μ + 4ν)
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

    @test @inferred(Hemiplex(1, 2, 3)) === 1 + 2μ + 3ν
    @test @inferred(Hemiplex(1.0, 2, 3)) === 1.0 + 2.0μ + 3.0ν
    @test @inferred(Hemiplex(2μ + 3ν)) === Hemiplex(0, 2, 3)
    @test @inferred(Hemiplex(7.0, 2μ + 3ν)) === Hemiplex(7.0, 2.0, 3.0)
    @test @inferred(Hemiplex(7)) === Hemiplex(7, 0, 0)

    @test convert(Hemiplex{Float64}, 1) === Hemiplex(1.0, 0.0, 0.0)
    @test convert(Hemiplex{Float64}, 2μ + 3ν) === Hemiplex(0.0, 2.0, 3.0)
    @test convert(Hemiplex{Float64}, Hemiplex(1, 2, 3)) === Hemiplex(1.0, 2.0, 3.0)
    @test convert(Float64, Hemiplex(1, 0, 0)) === 1.0
    @test_throws DomainError convert(Float64, Hemiplex(1, 2, 3))
    @test convert(PureHemi{Float64}, Hemiplex(0, 2, 3)) === PureHemi(2.0, 3.0)
    @test_throws DomainError convert(PureHemi{Float64}, Hemiplex(1, 2, 3))

    @test @inferred(-(7 + 2μ)) === Hemiplex(-7, -2, 0)
    @test @inferred(+(7 + 2μ)) === Hemiplex(7, 2, 0)
    @test @inferred((1 + 2μ + 3ν) + (5.5 + 2.1μ + 3.2ν)) === Hemiplex(6.5, 4.1, 6.2)
    @test (1 + 2μ + 3ν) - (5.5 + 2.1μ + 3.2ν) ≈ Hemiplex(-4.5, -0.1, -0.2)
    @test (1 + 2μ + 3ν) * (5.5 + 2.1μ + 3.2ν) ≈ Hemiplex(5.5 - (6.4 + 6.3) / 2, 13.1, 19.7)

    @test @inferred((1 + 2μ + 3ν) + PureHemi(2.1, 3.2)) === Hemiplex(1, 4.1, 6.2)
    @test @inferred(PureHemi(7, 9) + (2 + μ - 6ν)) === Hemiplex(2, 8, 3)
    @test (1 + 2μ + 3ν) - PureHemi(2.1, 3.2) ≈ Hemiplex(1.0, -0.1, -0.2)
    @test @inferred(PureHemi(7, 9) - (2 + μ - 6ν)) === Hemiplex(-2, 6, 15)
    @test (1 + 2μ + 3ν) * PureHemi(2.1, 3.2) ≈ Hemiplex(-(6.4 + 6.3) / 2, 2.1, 3.2)
    @test @inferred(PureHemi(7, 9) * (2 + μ - 6ν)) === Hemiplex(16.5, 14, 18)

    @test @inferred((1 + 2μ + 3ν) + 5.5) === Hemiplex(6.5, 2.0, 3.0)
    @test @inferred(3 + (2 + μ - 6ν)) === Hemiplex(5, 1, -6)
    @test @inferred((1 + 2μ + 3ν) - 5.5) === Hemiplex(-4.5, 2.0, 3.0)
    @test @inferred(3 - (2 + μ - 6ν)) === Hemiplex(1, -1, 6)
    @test @inferred(false * (2 + μ - 6ν)) === Hemiplex(0, 0, 0)
    @test @inferred((2 + μ - 6ν) * true) === Hemiplex(2, 1, -6)
    @test @inferred(1.5 * (2 + μ - 6ν)) === Hemiplex(3.0, 1.5, -9.0)
    @test @inferred((2 + μ - 6ν) * 1.5) === Hemiplex(3.0, 1.5, -9.0)
    @test @inferred((2 + μ - 6ν) / 1.5) === Hemiplex(2 / 1.5, 1 / 1.5, -6 / 1.5)
    @test @inferred(1.5 \ (2 + μ - 6ν)) === Hemiplex(2 / 1.5, 1 / 1.5, -6 / 1.5)

    x = 3.2 + μ - ν
    @test real(x) === 3.2
    @test real(Hemiplex{Float16}) == Float16
    @test mu(x) == 1
    @test nu(x) == -1
    @test mu([x]) == [1]
    @test nu([x]) == [-1]
    @test conj(x) === 3.2 - μ + ν
    @test zero(x) === Hemiplex(0.0, 0.0, 0.0)
    @test zero(Hemiplex{Float16}) === Hemiplex{Float16}(0, 0, 0)

    @testset "Hemiplex predicates" begin
        @test iszero(Hemiplex(0, 0, 0))
        @test !iszero(Hemiplex(1, 0, 0))
        @test !iszero(Hemiplex(0, 1, 0))
        @test isnan(Hemiplex(NaN, 0.0, 0.0))
        @test isnan(Hemiplex(0.0, NaN, 0.0))
        @test isnan(Hemiplex(0.0, 0.0, NaN))
        @test !isnan(Hemiplex(1.0, 2.0, 3.0))
        @test isinf(Hemiplex(Inf, 0.0, 0.0))
        @test isinf(Hemiplex(0.0, Inf, 0.0))
        @test isinf(Hemiplex(0.0, 0.0, Inf))
        @test !isinf(Hemiplex(1.0, 2.0, 3.0))
        @test isreal(Hemiplex(1.0, 0.0, 0.0))
        @test !isreal(Hemiplex(1.0, 1.0, 0.0))
        @test isinteger(Hemiplex(2, 0, 0))
        @test !isinteger(Hemiplex(2, 1, 0))
        @test !isinteger(Hemiplex(2.5, 0, 0))
        @test isone(Hemiplex(1, 0, 0))
        @test !isone(Hemiplex(2, 0, 0))
        @test !isone(Hemiplex(1, 1, 0))
    end

    @testset "Hemiplex one" begin
        @test one(Hemiplex{Int}) === Hemiplex(1, 0, 0)
        @test one(Hemiplex(3.0, 1.0, 2.0)) === Hemiplex(1.0, 0.0, 0.0)
        @test isone(one(Hemiplex{Float64}))
    end

    @testset "Hemiplex equality" begin
        @test Hemiplex(1, 2, 3) == Hemiplex(1, 2, 3)
        @test Hemiplex(1, 0, 0) == 1
        @test 1 == Hemiplex(1, 0, 0)
        @test Hemiplex(1, 1, 0) != 1
        @test isequal(Hemiplex(1, 2, 3), Hemiplex(1, 2, 3))
        @test isequal(Hemiplex(NaN, 0.0, 0.0), Hemiplex(NaN, 0.0, 0.0))
    end

    @testset "Hemiplex float/widen/big/bswap" begin
        @test float(Hemiplex(1, 2, 3)) === Hemiplex(1.0, 2.0, 3.0)
        @test float(Hemiplex(1.0, 2.0, 3.0)) === Hemiplex(1.0, 2.0, 3.0)
        @test float(Hemiplex{Int}) == Hemiplex{Float64}
        @test float(Hemiplex{Float32}) == Hemiplex{Float32}
        @test widen(Hemiplex{Float32}) == Hemiplex{Float64}
        @test widen(Hemiplex{Int32}) == Hemiplex{Int64}
        @test big(Hemiplex{Int}) == Hemiplex{BigInt}
        @test big(Hemiplex(1, 2, 3)) == Hemiplex(BigInt(1), BigInt(2), BigInt(3))
        @test bswap(Hemiplex(1, 2, 3)) == Hemiplex(bswap(1), bswap(2), bswap(3))
    end

    @testset "Hemiplex hash" begin
        @test hash(Hemiplex(0, 0, 0)) == hash(0)
        @test hash(Hemiplex(3, 0, 0)) == hash(3)
        @test hash(Hemiplex(3.0, 0.0, 0.0)) == hash(3)
        @test hash(Hemiplex(1, 2, 3)) == hash(Hemiplex(1, 2, 3))
        @test hash(Hemiplex(1, 2, 3)) != hash(Hemiplex(1, 3, 2))
        # Consistent with PureHemi when real part is zero
        @test hash(Hemiplex(0, 2, 3)) == hash(PureHemi(2, 3))
    end

    @testset "binary I/O" begin
        buf = IOBuffer()
        ph = PureHemi(1.5, -2.5)
        write(buf, ph)
        seekstart(buf)
        @test read(buf, PureHemi{Float64}) === ph

        buf = IOBuffer()
        hr = Hemiplex(3.0, 1.5, -2.5)
        write(buf, hr)
        seekstart(buf)
        @test read(buf, Hemiplex{Float64}) === hr
    end

    a = [μ + 2ν, 5μ - ν]
    @test mu(a) == [1, 5]
    @test nu(a) == [2, -1]
    @test @inferred(a * a') == [2 4.5; 4.5 -5]
    @test isa(a * a', Matrix{Float64})
    @test @inferred(a' * a) === -3.0
    A = [μ 0; 2μ + 3ν 5μ - 4ν]
    @test @inferred(A * a) == [-1, 9]
    @test isa(A * a, Vector{Hemiplex{Float64}})
    z = zero(PureHemi{Int})
    A = [μ z; 2μ + 3ν 5μ - 4ν]
    @test @inferred(A * a) == [-1, 9]
    @test isa(A * a, Vector{Float64})
    @test @inferred(A * A) == [0 0; -5 20]
    @test isa(A * A, Matrix{Float64})
    @test @inferred(A * A') == [0 1.5; 1.5 -14]
    @test isa(A * A', Matrix{Float64})
    @test @inferred(A * transpose(A)) == [0 -1.5; -1.5 14]
    @test isa(A * transpose(A), Matrix{Float64})
    @test @inferred(A' * A) == [6 3.5; 3.5 -20]
    @test isa(A' * A, Matrix{Float64})
    @test @inferred(transpose(A) * A) == [-6 -3.5; -3.5 20]
    @test isa(transpose(A) * A, Matrix{Float64})
    @test @inferred(A' * A') == [0 -5; 0 20]
    @test isa(A' * A', Matrix{Float64})
    @test @inferred(transpose(A) * transpose(A)) == [0 -5; 0 20]
    @test isa(transpose(A) * transpose(A), Matrix{Float64})

    @test @inferred(A * [1, 2]) == [μ, 12μ - 5ν]
    @test isa(A * [1, 2], Vector{PureHemi{Int}})

    @test μ * a == [-1, 0.5]
    @test isa(μ * a, Vector{Float64})
    @test 3 * a == [3μ + 6ν, 15μ - 3ν]
    @test isa(3 * a, Vector{PureHemi{Int}})

    # Promotion
    @test isa([μ, false], Vector{Hemiplex{Bool}})
    @test isa([2μ, false], Vector{Hemiplex{Int}})
    @test isa([μ, π], Vector{Hemiplex{Float64}})
    @test isa([1 + μ, false], Vector{Hemiplex{Int}})
    @test isa([1 + μ, π], Vector{Hemiplex{Float64}})
    @test isa([1 + μ, 2.0 + ν], Vector{Hemiplex{Float64}})

    @testset "complex coefficients" begin
        z = @inferred(PureHemi(1.0 + 2.0im, 3.0 - 1.0im))
        @test z isa PureHemi{ComplexF64}
        @test @inferred(PureHemi(1.0 + 0.0im, 2)) === PureHemi(1.0 + 0.0im, 2.0 + 0.0im)
        @test PureHemi{ComplexF32}(1, 0) === PureHemi(1.0f0 + 0.0f0im, 0.0f0 + 0.0f0im)
        @test convert(PureHemi{ComplexF64}, PureHemi(1, 2)) === PureHemi(1.0 + 0.0im, 2.0 + 0.0im)
        @test @inferred(2 * z) === PureHemi(2.0 + 4.0im, 6.0 - 2.0im)
        @test @inferred((1 + 1im) * z) === PureHemi(-1.0 + 3.0im, 4.0 + 2.0im)
        @test @inferred(z + PureHemi(1, 1)) === PureHemi(2.0 + 2.0im, 4.0 - 1.0im)

        @test coefftype(z) === ComplexF64
        @test coefftype(PureHemi{ComplexF32}) === ComplexF32
        @test coefftype(Hemiplex{ComplexF64}) === ComplexF64
        @test coefftype(Hemiplex(1.0, 2.0, 3.0)) === Float64
        @test real(PureHemi{ComplexF64}) === Float64
        @test real(Hemiplex{ComplexF32}) === Float32
        @test real(z) === 0.0
        @test real(Hemiplex(1.0 + 2.0im, z)) === 1.0

        # Promotion with complex scalars
        @test promote_type(ComplexF64, PureHemi{Float64}) === Hemiplex{ComplexF64}
        @test isa([PureHemi(1.0, 2.0), 1.0 + 0im], Vector{Hemiplex{ComplexF64}})
        @test @inferred(z + (1 + 1im)) === Hemiplex(1.0 + 1.0im, 1.0 + 2.0im, 3.0 - 1.0im)
        @test convert(Hemiplex{ComplexF64}, 1 + 1im) === Hemiplex(1.0 + 1.0im, 0.0im, 0.0im)

        # Conjugation
        @test conj(z) === PureHemi(-1.0 + 2.0im, -3.0 - 1.0im)
        @test conj(conj(z)) === z
        w = PureHemi(0.5, 1.0 - 1.0im)
        @test (z * w)' == z' * w'
        @test conj(Hemiplex(1.0 + 2.0im, z)) === Hemiplex(1.0 - 2.0im, conj(z))
        L = [PureHemi(1.0 + 1.0im, 2.0 - 1.0im) PureHemi(0.5 + 0.0im, 0.0 + 1.0im);
             PureHemi(2.0 + 0.0im, 3.0 + 1.0im) PureHemi(0.0 + 1.0im, 1.0 + 0.0im)]
        Ladj = L'
        for i in axes(Ladj, 1), j in axes(Ladj, 2)
            @test Ladj[i, j] == conj(L[j, i])
        end

        # abs2 is a real quadratic form
        @test abs2(z) isa Float64
        @test abs2(z) == real(z.m * conj(z.n))
        @test z * conj(z) ≈ abs2(z)
        @test abs2(PureHemi(3, 4)) == 12

        # Symmetric division
        for c in (2.0 - 3.0im, 3.0)
            @test z * (c / z) ≈ c
            @test z * (z \ c) ≈ c
        end

        hc = Hemiplex(1.0 + 1.0im, 2.0 + 0.0im, 3.0 + 0.0im)
        @test_throws "not scalar-valued" abs2(hc)
        @test_throws "not scalar-valued" abs(hc)
        # inv uses the coefficient-linear hemi-conjugate
        @test hc * inv(hc) ≈ 1
        @test inv(inv(hc)) ≈ hc
        let zc = Hemiplex(0.0im, 2.0 + 1.0im, 1.0 - 3.0im)
            @test zc * inv(zc) ≈ 1
        end
        let zd = Hemiplex(1.0 + 0.0im, 1.0im, 1.0im)  # r² + mn = 0
            @test !isfinite(inv(zd))
        end
        @test inv(Hemiplex(1.0, 2.0, 3.0)) === Hemiplex(1.0, -2.0, -3.0) / 7.0

        # isapprox uses the coefficient sup-norm, so it works even where abs does not
        @test Hemiplex(0.0, 1.0, -1.0) ≈ Hemiplex(0.0, 1.0, -1.0)      # abs2 < 0
        @test PureHemi(1.0 + 1.0im, 2.0) ≈ PureHemi(1.0 + 1.0im, 2.0 + 1e-13)
        @test !(Hemiplex(1.0, 0.0, 0.0) ≈ Hemiplex(1.0, 0.1, 0.0))
        @test Hemiplex(1.0 + 0.0im, 0.0im, 0.0im) ≈ 1.0

        @test !isreal(Hemiplex(1.0 + 2.0im, 0.0im, 0.0im))
        @test isreal(Hemiplex(1.0 + 0.0im, 0.0im, 0.0im))
        @test convert(ComplexF64, Hemiplex(1.0 + 2.0im, 0.0im, 0.0im)) === 1.0 + 2.0im
        @test_throws InexactError convert(Float64, Hemiplex(1.0 + 2.0im, 0.0im, 0.0im))
        @test_throws DomainError convert(ComplexF64, Hemiplex(1.0 + 2.0im, 1.0im, 0.0im))

        # symmetric division by a real pure-hemi must not round-trip through abs2
        @test 2.0 / PureHemi(1e300, 1e300) === PureHemi(-2.0e-300, -2.0e-300)
        @test (3.0 + 1.0im) / PureHemi(1e300, 1e300) === PureHemi(-((3.0 + 1.0im) / 1e300), -((3.0 + 1.0im) / 1e300))

        @test mu([z]) == [1.0 + 2.0im]
        @test nu([z]) == [3.0 - 1.0im]
        @test sprint(show, z) == "(1.0 + 2.0im)μ + (3.0 - 1.0im)ν"
        @test sprint(show, Hemiplex(1.0 + 0.0im, z)) == "(1.0 + 0.0im) + (1.0 + 2.0im)μ + (3.0 - 1.0im)ν"

        buf = IOBuffer()
        write(buf, z)
        seekstart(buf)
        @test read(buf, PureHemi{ComplexF64}) === z
    end

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
