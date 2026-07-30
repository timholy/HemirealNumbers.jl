module HemiplexNumbers

import Base: +, -, *, /, \, ^, abs, abs2, big, bswap, conj, convert,
             float, hash, inv, isfinite, isinteger, isnan, isinf, isone, isreal,
             iszero, one, promote_rule, read, real, show, show_unquoted,
             widen, write, zero

export PureHemi, Hemiplex, coefftype, μ, ν, mu, nu

# The coefficient ring must be commutative and closed under `conj`.
const RealOrComplex = Union{Real, Complex}

struct PureHemi{T <: RealOrComplex} <: Number
    m::T
    n::T
end
PureHemi(m::RealOrComplex, n::RealOrComplex) = PureHemi(promote(m, n)...)

struct Hemiplex{T <: RealOrComplex} <: Number
    r::T
    h::PureHemi{T}
end

Hemiplex{T}(r::RealOrComplex, m::RealOrComplex, n::RealOrComplex) where {T <: RealOrComplex} = Hemiplex{T}(r, PureHemi{T}(m, n))
Hemiplex(r::T, m::T, n::T) where {T <: RealOrComplex} = Hemiplex{T}(r, PureHemi{T}(m, n))
Hemiplex(r::RealOrComplex, m::RealOrComplex, n::RealOrComplex) = Hemiplex(promote(r, m, n)...)

function Hemiplex(r::R, h::PureHemi{H}) where {R <: RealOrComplex, H <: RealOrComplex}
    T = promote_type(R, H)
    return Hemiplex{T}(r, PureHemi{T}(h.m, h.n))
end
Hemiplex(h::PureHemi{R}) where {R <: RealOrComplex} = Hemiplex(zero(R), h)
Hemiplex(r::R) where {R <: RealOrComplex} = Hemiplex(r, zero(PureHemi{R}))

const μ = PureHemi(true, false)
const ν = PureHemi(false, true)

"""
    coefftype(x) -> T
    coefftype(::Type{PureHemi{T}}) -> T

Return the type of the coefficients multiplying `μ` and `ν` (and, for [`Hemiplex`](@ref),
the scalar part). It may be real or complex; `real(typeof(x))` instead returns the real
type underlying those coefficients.
"""
coefftype(::Type{PureHemi{R}}) where {R <: RealOrComplex} = R
coefftype(::Type{Hemiplex{R}}) where {R <: RealOrComplex} = R
coefftype(x::Union{PureHemi, Hemiplex}) = coefftype(typeof(x))

## PureHemi implementation
convert(::Type{PureHemi{R}}, x::PureHemi) where {R} = PureHemi{R}(x.m, x.n)
convert(::Type{PureHemi{R}}, x::RealOrComplex) where {R} =
    iszero(x) ? zero(PureHemi{R}) : throw(DomainError(x, "Non-zero scalars cannot be converted to pure-hemi numbers"))

convert(::Type{R}, x::Hemiplex) where {R <: RealOrComplex} =
    iszero(x.h) ? R(x.r) : throw(DomainError(x, "Hemiplex numbers with non-zero hemi-part cannot be converted to scalars"))
convert(::Type{R}, x::PureHemi) where {R <: RealOrComplex} =
    iszero(x) ? zero(R) : throw(DomainError(x, "PureHemi numbers with non-zero components cannot be converted to scalars"))

iszero(x::PureHemi) = iszero(x.m) & iszero(x.n)
isnan(x::PureHemi) = isnan(x.m) | isnan(x.n)
isinf(x::PureHemi) = isinf(x.m) | isinf(x.n)
isfinite(x::PureHemi) = isfinite(x.m) & isfinite(x.n)
isreal(x::PureHemi) = iszero(x)

Base.:(==)(x::PureHemi, y::PureHemi) = (x.m == y.m) & (x.n == y.n)
Base.isequal(x::PureHemi, y::PureHemi) = isequal(x.m, y.m) & isequal(x.n, y.n)

(+)(x::PureHemi) = PureHemi(+x.m, +x.n)
(-)(x::PureHemi) = PureHemi(-x.m, -x.n)

(+)(x::PureHemi, y::PureHemi) = PureHemi(x.m + y.m, x.n + y.n)
(-)(x::PureHemi, y::PureHemi) = PureHemi(x.m - y.m, x.n - y.n)
(*)(x::PureHemi, y::PureHemi) = -(x.m * y.n + x.n * y.m) / 2
(*)(c::Bool, x::PureHemi) = PureHemi(c * x.m, c * x.n)
(*)(c::RealOrComplex, x::PureHemi) = PureHemi(c * x.m, c * x.n)
(*)(x::PureHemi, c::RealOrComplex) = c * x
(/)(x::PureHemi, c::RealOrComplex) = PureHemi(x.m / c, x.n / c)
(\)(c::RealOrComplex, x::PureHemi) = x / c
function (/)(x::PureHemi, y::PureHemi)
    x.m * y.n == x.n * y.m || throw(ArgumentError("$x is not a scalar multiple of $y"))
    return y.n == 0 ? x.m / y.m : x.n / y.n
end
(\)(x::PureHemi, y::PureHemi) = y / x
(^)(x::PureHemi, p::Integer) = (x * x)^(p / 2)  # note x*x is much faster than x^2
(^)(x::PureHemi, p::Rational) = (x * x)^(p / 2)
(^)(x::PureHemi, p::Real) = (x * x)^(p / 2)

# Symmetric division (there are other solutions y to y*x = c).
# For real x the coefficients cancel, avoiding spurious overflow of abs2(x).
(/)(c::RealOrComplex, x::PureHemi{<:Real}) = PureHemi(-c / x.n, -c / x.m)
(/)(c::RealOrComplex, x::PureHemi) = c * conj(x) / abs2(x)
(\)(x::PureHemi, c::RealOrComplex) = c / x

real(::Type{PureHemi{R}}) where {R <: RealOrComplex} = real(R)
real(x::PureHemi{R}) where {R <: RealOrComplex} = zero(real(R))
conj(x::PureHemi) = PureHemi(-conj(x.m), -conj(x.n))
mu(x::PureHemi) = x.m
nu(x::PureHemi) = x.n

zero(::Type{PureHemi{R}}) where {R <: RealOrComplex} = PureHemi{R}(0, 0)
zero(::PureHemi{R}) where {R <: RealOrComplex} = zero(PureHemi{R})

# (note that abs2(x) != x*x for PureHemi, and that abs2(x) may be negative)
abs2(x::PureHemi) = real(x.m * conj(x.n))
abs(x::PureHemi) = sqrt(abs2(x))

float(::Type{PureHemi{T}}) where {T <: AbstractFloat} = PureHemi{T}
float(::Type{PureHemi{T}}) where {T} = PureHemi{float(T)}
float(x::PureHemi{<:AbstractFloat}) = x
float(x::PureHemi) = PureHemi(float(x.m), float(x.n))

widen(::Type{PureHemi{T}}) where {T} = PureHemi{widen(T)}

big(::Type{PureHemi{T}}) where {T <: RealOrComplex} = PureHemi{big(T)}
big(x::PureHemi{T}) where {T <: RealOrComplex} = PureHemi{big(T)}(x.m, x.n)

bswap(x::PureHemi) = PureHemi(bswap(x.m), bswap(x.n))

function hash(x::PureHemi, h::UInt)
    # A zero PureHemi is equal to zero real
    iszero(x) && return hash(zero(x.m), h)
    hash(x.m, hash(x.n, h))
end

## Hemiplex implementation
convert(::Type{Hemiplex{R}}, x::Hemiplex) where {R <: RealOrComplex} = Hemiplex{R}(x.r, x.h)
convert(::Type{Hemiplex{R}}, x::PureHemi) where {R <: RealOrComplex} = Hemiplex{R}(0, x)
convert(::Type{Hemiplex{R}}, x::RealOrComplex) where {R <: RealOrComplex} = Hemiplex{R}(x, zero(PureHemi{R}))
convert(::Type{PureHemi{R}}, x::Hemiplex) where {R <: RealOrComplex} =
    iszero(x.r) ? convert(PureHemi{R}, x.h) : throw(DomainError(x, "Hemiplex numbers with non-zero real part cannot be converted to pure-hemi numbers"))

iszero(x::Hemiplex) = iszero(x.r) & iszero(x.h)
isnan(x::Hemiplex) = isnan(x.r) | isnan(x.h)
isinf(x::Hemiplex) = isinf(x.r) | isinf(x.h)
isfinite(x::Hemiplex) = isfinite(x.r) & isfinite(x.h)
isreal(x::Hemiplex) = isreal(x.r) & iszero(x.h)
isinteger(x::Hemiplex) = isreal(x) & isinteger(x.r)
isone(x::Hemiplex) = isone(x.r) & iszero(x.h)

Base.:(==)(x::Hemiplex, y::Hemiplex) = (x.r == y.r) & (x.h == y.h)
Base.isequal(x::Hemiplex, y::Hemiplex) = isequal(x.r, y.r) & isequal(x.h, y.h)

(+)(x::Hemiplex) = Hemiplex(+x.r, +x.h)
(-)(x::Hemiplex) = Hemiplex(-x.r, -x.h)

(+)(x::Hemiplex, y::Hemiplex) = Hemiplex(x.r + y.r, x.h + y.h)
(-)(x::Hemiplex, y::Hemiplex) = Hemiplex(x.r - y.r, x.h - y.h)
(*)(x::Hemiplex, y::Hemiplex) = Hemiplex(x.r * y.r + x.h * y.h, x.r * y.h + x.h * y.r)

(+)(x::Hemiplex, y::PureHemi) = x + Hemiplex(y)
(+)(x::PureHemi, y::Hemiplex) = y + x
(-)(x::Hemiplex, y::PureHemi) = x - Hemiplex(y)
(-)(x::PureHemi, y::Hemiplex) = Hemiplex(x) - y
(*)(x::Hemiplex, y::PureHemi) = x * Hemiplex(y)
(*)(x::PureHemi, y::Hemiplex) = y * x

(+)(x::Hemiplex, y::RealOrComplex) = x + Hemiplex(y)
(+)(x::RealOrComplex, y::Hemiplex) = y + x
(-)(x::Hemiplex, y::RealOrComplex) = x - Hemiplex(y)
(-)(x::RealOrComplex, y::Hemiplex) = Hemiplex(x) - y
(*)(c::Bool, x::Hemiplex) = Hemiplex(c * x.r, c * x.h)
(*)(c::RealOrComplex, x::Hemiplex) = Hemiplex(c * x.r, c * x.h)
(*)(x::Hemiplex, c::RealOrComplex) = c * x
(/)(x::Hemiplex, c::RealOrComplex) = Hemiplex(x.r / c, x.h / c)
(\)(c::RealOrComplex, x::Hemiplex) = x / c

real(::Type{Hemiplex{R}}) where {R <: RealOrComplex} = real(R)
real(x::Hemiplex{R}) where {R <: RealOrComplex} = real(x.r)
conj(x::Hemiplex) = Hemiplex(conj(x.r), conj(x.h))
mu(x::Hemiplex) = mu(x.h)
nu(x::Hemiplex) = nu(x.h)

zero(::Type{Hemiplex{R}}) where {R <: RealOrComplex} = Hemiplex{R}(0, PureHemi{R}(0, 0))
zero(::Hemiplex{R}) where {R <: RealOrComplex} = zero(Hemiplex{R})
one(::Type{Hemiplex{R}}) where {R <: RealOrComplex} = Hemiplex{R}(1, PureHemi{R}(0, 0))
one(::Hemiplex{R}) where {R <: RealOrComplex} = one(Hemiplex{R})

# (note that abs2(x) != x*x for the hemi part)
abs2(x::Hemiplex{<:Real}) = x.r * x.r + abs2(x.h)
abs(x::Hemiplex{<:Real}) = sqrt(abs2(x))

# With complex coefficients z*conj(z) retains a hemi-part, so these have no scalar answer.
abs2(x::Hemiplex{<:Complex}) =
    throw(ArgumentError("abs2(z) = z*conj(z) is not scalar-valued for complex-coefficient Hemiplex numbers"))
abs(x::Hemiplex{<:Complex}) =
    throw(ArgumentError("abs(z) = sqrt(z*conj(z)) is not defined for complex-coefficient Hemiplex numbers: z*conj(z) is not scalar-valued"))

# The coefficient-linear hemi-conjugate q = r - mμ - nν satisfies
# z*q = r² + mn even for complex coefficients.
inv(x::Hemiplex) = Hemiplex(x.r, -x.h) / (x.r * x.r + x.h.m * x.h.n)

# `abs` is not a norm on the hemiplexes (abs2 may be negative, or non-scalar for complex
# coefficients), so approximate equality uses the sup-norm of the coefficients instead.
_infnorm(x::Hemiplex) = max(abs(x.r), abs(x.h.m), abs(x.h.n))
function Base.isapprox(x::Hemiplex, y::Hemiplex;
                       atol::Real=0,
                       rtol::Real=Base.rtoldefault(real(coefftype(x)), real(coefftype(y)), atol))
    return _infnorm(x - y) <= max(atol, rtol * max(_infnorm(x), _infnorm(y)))
end
Base.isapprox(x::PureHemi, y::PureHemi; kwargs...) = isapprox(Hemiplex(x), Hemiplex(y); kwargs...)
Base.isapprox(x::Union{PureHemi, Hemiplex}, y::RealOrComplex; kwargs...) = isapprox(promote(x, y)...; kwargs...)
Base.isapprox(x::RealOrComplex, y::Union{PureHemi, Hemiplex}; kwargs...) = isapprox(y, x; kwargs...)

float(::Type{Hemiplex{T}}) where {T <: AbstractFloat} = Hemiplex{T}
float(::Type{Hemiplex{T}}) where {T} = Hemiplex{float(T)}
float(x::Hemiplex{<:AbstractFloat}) = x
float(x::Hemiplex) = Hemiplex(float(x.r), float(x.h))

widen(::Type{Hemiplex{T}}) where {T} = Hemiplex{widen(T)}

big(::Type{Hemiplex{T}}) where {T <: RealOrComplex} = Hemiplex{big(T)}
big(x::Hemiplex{T}) where {T <: RealOrComplex} = Hemiplex{big(T)}(x.r, x.h)

bswap(x::Hemiplex) = Hemiplex(bswap(x.r), bswap(x.h))

function hash(x::Hemiplex, h::UInt)
    # A real-valued Hemiplex hashes like its real part
    iszero(x.h) && return hash(x.r, h)
    # A pure-hemi Hemiplex hashes like its PureHemi
    iszero(x.r) && return hash(x.h, h)
    hash(x.r, hash(x.h, h))
end

for f in (:mu, :nu)
    @eval begin
        function $f(X::AbstractArray)
            out = Array{coefftype(eltype(X))}(undef, size(X))
            for I in eachindex(X)
                out[I] = $f(X[I])
            end
            return out
        end
    end
end

promote_rule(::Type{Bool}, ::Type{PureHemi{H}}) where {H <: RealOrComplex} = Hemiplex{promote_type(Bool, H)}
promote_rule(::Type{R}, ::Type{PureHemi{H}}) where {R <: AbstractIrrational, H <: RealOrComplex} = Hemiplex{promote_type(R, H)}
promote_rule(::Type{R}, ::Type{PureHemi{H}}) where {R <: RealOrComplex, H <: RealOrComplex} = Hemiplex{promote_type(R, H)}
promote_rule(::Type{PureHemi{H1}}, ::Type{PureHemi{H2}}) where {H1 <: RealOrComplex, H2 <: RealOrComplex} = PureHemi{promote_type(H1, H2)}
promote_rule(::Type{Bool}, ::Type{Hemiplex{H}}) where {H <: RealOrComplex} = Hemiplex{promote_type(Bool, H)}
promote_rule(::Type{R}, ::Type{Hemiplex{H}}) where {R <: AbstractIrrational, H <: RealOrComplex} = Hemiplex{promote_type(R, H)}
promote_rule(::Type{R}, ::Type{Hemiplex{H}}) where {R <: RealOrComplex, H <: RealOrComplex} = Hemiplex{promote_type(R, H)}
promote_rule(::Type{PureHemi{H1}}, ::Type{Hemiplex{H2}}) where {H1 <: RealOrComplex, H2 <: RealOrComplex} = Hemiplex{promote_type(H1, H2)}
promote_rule(::Type{Hemiplex{H1}}, ::Type{Hemiplex{H2}}) where {H1 <: RealOrComplex, H2 <: RealOrComplex} = Hemiplex{promote_type(H1, H2)}

## Show

# Use the IOBuffer trick (as in Base complex.jl) to robustly detect sign of a value
# without assuming anything about its string representation beyond the leading '-'.
# Complex coefficients carry internal signs, so they are parenthesized and always
# joined with '+'.

function show_coefficient(io::IO, c::Complex)
    print(io, "(")
    show(io, c)
    print(io, ")")
end
show_coefficient(io::IO, c::Real) = show(io, c)

function show_signed_coefficient(io::IO, c::Complex, compact::Bool)
    print(io, compact ? "+" : " + ")
    show_coefficient(io, c)
end
function show_signed_coefficient(io::IO, c::Real, compact::Bool)
    bufio = IOBuffer()
    show(IOContext(bufio, io), c)
    seekstart(bufio)
    if peek(bufio) === UInt8('-')
        seek(bufio, 1)
        print(io, compact ? "-" : " - ")
    else
        print(io, compact ? "+" : " + ")
    end
    write(io, bufio)
    return nothing
end

function show(io::IO, x::PureHemi)
    compact = get(io, :compact, false)::Bool
    show_coefficient(io, x.m)
    print(io, "μ")
    show_signed_coefficient(io, x.n, compact)
    print(io, "ν")
end

function show(io::IO, x::Hemiplex)
    compact = get(io, :compact, false)::Bool
    show_coefficient(io, x.r)
    show_signed_coefficient(io, x.h.m, compact)
    print(io, "μ")
    show_signed_coefficient(io, x.h.n, compact)
    print(io, "ν")
end

function show_unquoted(io::IO, x::Union{PureHemi, Hemiplex}, ::Int, prec::Int)
    if Base.operator_precedence(:+) <= prec
        print(io, "(")
        show(io, x)
        print(io, ")")
    else
        show(io, x)
    end
end

## Binary I/O

function read(s::IO, ::Type{PureHemi{T}}) where {T <: RealOrComplex}
    m = read(s, T)
    n = read(s, T)
    PureHemi{T}(m, n)
end
function write(s::IO, x::PureHemi)
    write(s, x.m, x.n)
end

function read(s::IO, ::Type{Hemiplex{T}}) where {T <: RealOrComplex}
    r = read(s, T)
    h = read(s, PureHemi{T})
    Hemiplex{T}(r, h)
end
function write(s::IO, x::Hemiplex)
    write(s, x.r, x.h)
end

end # module
