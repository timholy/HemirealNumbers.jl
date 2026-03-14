module HemirealNumbers

import Base: +, -, *, /, \, ^, abs, abs2, big, bswap, conj, convert,
             float, hash, isfinite, isinteger, isnan, isinf, isone, isreal,
             iszero, one, promote_rule, read, real, show, show_unquoted,
             widen, write, zero

export PureHemi, Hemireal, μ, ν, mu, nu

struct PureHemi{T <: Real} <: Number
    m::T
    n::T
end
PureHemi(m::Real, n::Real) = PureHemi(promote(m, n)...)

struct Hemireal{T <: Real} <: Number
    r::T
    h::PureHemi{T}
end

Hemireal{T}(r::Real, m::Real, n::Real) where {T <: Real} = Hemireal{T}(r, PureHemi{T}(m, n))
Hemireal(r::T, m::T, n::T) where {T <: Real} = Hemireal{T}(r, PureHemi{T}(m, n))
Hemireal(r::Real, m::Real, n::Real) = Hemireal(promote(r, m, n)...)

function Hemireal(r::R, h::PureHemi{H}) where {R <: Real, H <: Real}
    T = promote_type(R, H)
    return Hemireal{T}(r, PureHemi{T}(h.m, h.n))
end
Hemireal(h::PureHemi{R}) where {R <: Real} = Hemireal(zero(R), h)
Hemireal(r::R) where {R <: Real} = Hemireal(r, zero(PureHemi{R}))

const μ = PureHemi(true, false)
const ν = PureHemi(false, true)

## PureHemi implementation
convert(::Type{PureHemi{R}}, x::PureHemi) where {R} = PureHemi{R}(x.m, x.n)
convert(::Type{PureHemi{R}}, x::Real) where {R} =
    iszero(x) ? zero(PureHemi{R}) : throw(DomainError(x, "Non-zero reals cannot be converted to pure-hemi numbers"))

convert(::Type{R}, x::Hemireal) where {R <: Real} =
    isreal(x) ? R(x.r) : throw(DomainError(x, "Hemireal numbers with non-zero hemi-part cannot be converted to real numbers"))
convert(::Type{R}, x::PureHemi) where {R <: Real} =
    iszero(x) ? zero(R) : throw(DomainError(x, "PureHemi numbers with non-zero components cannot be converted to real numbers"))

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
(*)(x::PureHemi, y::PureHemi) = x.m * y.n + x.n * y.m
(*)(c::Bool, x::PureHemi) = PureHemi(c * x.m, c * x.n)
(*)(c::Real, x::PureHemi) = PureHemi(c * x.m, c * x.n)
(*)(x::PureHemi, c::Real) = c * x
(/)(x::PureHemi, c::Real) = PureHemi(x.m / c, x.n / c)
(\)(c::Real, x::PureHemi) = x / c
function (/)(x::PureHemi, y::PureHemi)
    x.m * y.n == x.n * y.m || throw(ArgumentError("$x is not a scalar multiple of $y"))
    return y.n == 0 ? x.m / y.m : x.n / y.n
end
(\)(x::PureHemi, y::PureHemi) = y / x
(^)(x::PureHemi, p::Integer) = (x * x)^(p / 2)  # note x*x is much faster than x^2
(^)(x::PureHemi, p::Rational) = (x * x)^(p / 2)
(^)(x::PureHemi, p::Real) = (x * x)^(p / 2)

# Symmetric division (there are other solutions y to y*x = c)
(/)(c::Real, x::PureHemi) = PureHemi(c / (2 * x.n), c / (2 * x.m))
(\)(x::PureHemi, c::Real) = c / x

real(::Type{PureHemi{R}}) where {R <: Real} = R
real(x::PureHemi{R}) where {R <: Real} = zero(R)
conj(x::PureHemi) = x
mu(x::PureHemi) = x.m
nu(x::PureHemi) = x.n

zero(::Type{PureHemi{R}}) where {R <: Real} = PureHemi{R}(0, 0)
zero(::PureHemi{R}) where {R <: Real} = zero(PureHemi{R})

# (note that abs2(x) != x*x for PureHemi)
abs2(x::PureHemi) = x.m * x.m + x.n * x.n
abs(x::PureHemi) = sqrt(abs2(x))

float(::Type{PureHemi{T}}) where {T <: AbstractFloat} = PureHemi{T}
float(::Type{PureHemi{T}}) where {T} = PureHemi{float(T)}
float(x::PureHemi{<:AbstractFloat}) = x
float(x::PureHemi) = PureHemi(float(x.m), float(x.n))

widen(::Type{PureHemi{T}}) where {T} = PureHemi{widen(T)}

big(::Type{PureHemi{T}}) where {T <: Real} = PureHemi{big(T)}
big(x::PureHemi{T}) where {T <: Real} = PureHemi{big(T)}(x.m, x.n)

bswap(x::PureHemi) = PureHemi(bswap(x.m), bswap(x.n))

function hash(x::PureHemi, h::UInt)
    # A zero PureHemi is equal to zero real
    iszero(x) && return hash(zero(x.m), h)
    hash(x.m, hash(x.n, h))
end

## Hemireal implementation
convert(::Type{Hemireal{R}}, x::Hemireal) where {R <: Real} = Hemireal{R}(x.r, x.h)
convert(::Type{Hemireal{R}}, x::PureHemi) where {R <: Real} = Hemireal{R}(0, x)
convert(::Type{Hemireal{R}}, x::Real) where {R <: Real} = Hemireal{R}(x, zero(PureHemi{R}))
convert(::Type{Hemireal{R}}, r::Real, m::Real, n::Real) where {R <: Real} = Hemireal{R}(r, PureHemi{R}(m, n))

iszero(x::Hemireal) = iszero(x.r) & iszero(x.h)
isnan(x::Hemireal) = isnan(x.r) | isnan(x.h)
isinf(x::Hemireal) = isinf(x.r) | isinf(x.h)
isfinite(x::Hemireal) = isfinite(x.r) & isfinite(x.h)
isreal(x::Hemireal) = iszero(x.h)
isinteger(x::Hemireal) = isreal(x) & isinteger(x.r)
isone(x::Hemireal) = isone(x.r) & iszero(x.h)

Base.:(==)(x::Hemireal, y::Hemireal) = (x.r == y.r) & (x.h == y.h)
Base.isequal(x::Hemireal, y::Hemireal) = isequal(x.r, y.r) & isequal(x.h, y.h)

(+)(x::Hemireal) = Hemireal(+x.r, +x.h)
(-)(x::Hemireal) = Hemireal(-x.r, -x.h)

(+)(x::Hemireal, y::Hemireal) = Hemireal(x.r + y.r, x.h + y.h)
(-)(x::Hemireal, y::Hemireal) = Hemireal(x.r - y.r, x.h - y.h)
(*)(x::Hemireal, y::Hemireal) = Hemireal(x.r * y.r + x.h * y.h, x.r * y.h + x.h * y.r)

(+)(x::Hemireal, y::PureHemi) = x + Hemireal(y)
(+)(x::PureHemi, y::Hemireal) = y + x
(-)(x::Hemireal, y::PureHemi) = x - Hemireal(y)
(-)(x::PureHemi, y::Hemireal) = Hemireal(x) - y
(*)(x::Hemireal, y::PureHemi) = x * Hemireal(y)
(*)(x::PureHemi, y::Hemireal) = y * x

(+)(x::Hemireal, y::Real) = x + Hemireal(y)
(+)(x::Real, y::Hemireal) = y + x
(-)(x::Hemireal, y::Real) = x - Hemireal(y)
(-)(x::Real, y::Hemireal) = Hemireal(x) - y
(*)(c::Bool, x::Hemireal) = Hemireal(c * x.r, c * x.h)
(*)(c::Real, x::Hemireal) = Hemireal(c * x.r, c * x.h)
(*)(x::Hemireal, c::Real) = c * x
(/)(x::Hemireal, c::Real) = Hemireal(x.r / c, x.h / c)
(\)(c::Real, x::Hemireal) = x / c

real(::Type{Hemireal{R}}) where {R <: Real} = R
real(x::Hemireal{R}) where {R <: Real} = x.r
conj(x::Hemireal) = x
mu(x::Hemireal) = mu(x.h)
nu(x::Hemireal) = nu(x.h)

zero(::Type{Hemireal{R}}) where {R <: Real} = Hemireal{R}(0, PureHemi{R}(0, 0))
zero(::Hemireal{R}) where {R <: Real} = zero(Hemireal{R})
one(::Type{Hemireal{R}}) where {R <: Real} = Hemireal{R}(1, PureHemi{R}(0, 0))
one(::Hemireal{R}) where {R <: Real} = one(Hemireal{R})

# (note that abs2(x) != x*x for the hemi part)
abs2(x::Hemireal) = x.r * x.r + abs2(x.h)
abs(x::Hemireal) = sqrt(abs2(x))

float(::Type{Hemireal{T}}) where {T <: AbstractFloat} = Hemireal{T}
float(::Type{Hemireal{T}}) where {T} = Hemireal{float(T)}
float(x::Hemireal{<:AbstractFloat}) = x
float(x::Hemireal) = Hemireal(float(x.r), float(x.h))

widen(::Type{Hemireal{T}}) where {T} = Hemireal{widen(T)}

big(::Type{Hemireal{T}}) where {T <: Real} = Hemireal{big(T)}
big(x::Hemireal{T}) where {T <: Real} = Hemireal{big(T)}(x.r, x.h)

bswap(x::Hemireal) = Hemireal(bswap(x.r), bswap(x.h))

function hash(x::Hemireal, h::UInt)
    # A real-valued Hemireal hashes like its real part
    iszero(x.h) && return hash(x.r, h)
    # A pure-hemi Hemireal hashes like its PureHemi
    iszero(x.r) && return hash(x.h, h)
    hash(x.r, hash(x.h, h))
end

for f in (:mu, :nu)
    @eval begin
        function $f(X::AbstractArray)
            out = Array{real(eltype(X))}(undef, size(X))
            for I in eachindex(X)
                out[I] = $f(X[I])
            end
            return out
        end
    end
end

promote_rule(::Type{Bool}, ::Type{PureHemi{H}}) where {H <: Real} = Hemireal{promote_type(Bool, H)}
promote_rule(::Type{R}, ::Type{PureHemi{H}}) where {R <: AbstractIrrational, H <: Real} = Hemireal{promote_type(R, H)}
promote_rule(::Type{R}, ::Type{PureHemi{H}}) where {R <: Real, H <: Real} = Hemireal{promote_type(R, H)}
promote_rule(::Type{PureHemi{H1}}, ::Type{PureHemi{H2}}) where {H1 <: Real, H2 <: Real} = PureHemi{promote_type(H1, H2)}
promote_rule(::Type{Bool}, ::Type{Hemireal{H}}) where {H <: Real} = Hemireal{promote_type(Bool, H)}
promote_rule(::Type{R}, ::Type{Hemireal{H}}) where {R <: AbstractIrrational, H <: Real} = Hemireal{promote_type(R, H)}
promote_rule(::Type{R}, ::Type{Hemireal{H}}) where {R <: Real, H <: Real} = Hemireal{promote_type(R, H)}
promote_rule(::Type{PureHemi{H1}}, ::Type{Hemireal{H2}}) where {H1 <: Real, H2 <: Real} = Hemireal{promote_type(H1, H2)}
promote_rule(::Type{Hemireal{H1}}, ::Type{Hemireal{H2}}) where {H1 <: Real, H2 <: Real} = Hemireal{promote_type(H1, H2)}

## Show

# Use the IOBuffer trick (as in Base complex.jl) to robustly detect sign of a value
# without assuming anything about its string representation beyond the leading '-'.

function show(io::IO, x::PureHemi)
    compact = get(io, :compact, false)::Bool
    show(io, x.m)
    print(io, "μ")
    bufio = IOBuffer()
    show(IOContext(bufio, io), x.n)
    seekstart(bufio)
    if peek(bufio) === UInt8('-')
        seek(bufio, 1)
        print(io, compact ? "-" : " - ")
    else
        print(io, compact ? "+" : " + ")
    end
    write(io, bufio)
    print(io, "ν")
end

function show(io::IO, x::Hemireal)
    compact = get(io, :compact, false)::Bool
    show(io, x.r)
    bufio = IOBuffer()
    show(IOContext(bufio, io), x.h.m)
    seekstart(bufio)
    if peek(bufio) === UInt8('-')
        seek(bufio, 1)
        print(io, compact ? "-" : " - ")
    else
        print(io, compact ? "+" : " + ")
    end
    write(io, bufio)
    print(io, "μ")
    bufio2 = IOBuffer()
    show(IOContext(bufio2, io), x.h.n)
    seekstart(bufio2)
    if peek(bufio2) === UInt8('-')
        seek(bufio2, 1)
        print(io, compact ? "-" : " - ")
    else
        print(io, compact ? "+" : " + ")
    end
    write(io, bufio2)
    print(io, "ν")
end

function show_unquoted(io::IO, x::Union{PureHemi, Hemireal}, ::Int, prec::Int)
    if Base.operator_precedence(:+) <= prec
        print(io, "(")
        show(io, x)
        print(io, ")")
    else
        show(io, x)
    end
end

## Binary I/O

function read(s::IO, ::Type{PureHemi{T}}) where {T <: Real}
    m = read(s, T)
    n = read(s, T)
    PureHemi{T}(m, n)
end
function write(s::IO, x::PureHemi)
    write(s, x.m, x.n)
end

function read(s::IO, ::Type{Hemireal{T}}) where {T <: Real}
    r = read(s, T)
    h = read(s, PureHemi{T})
    Hemireal{T}(r, h)
end
function write(s::IO, x::Hemireal)
    write(s, x.r, x.h)
end

end # module
