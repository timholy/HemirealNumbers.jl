module HemiplexNumbers

import Base: +, -, *, /, \, ^, abs, abs2, big, bswap, conj, convert,
             float, hash, isfinite, isinteger, isnan, isinf, isone, isreal,
             iszero, one, promote_rule, read, real, show, show_unquoted,
             widen, write, zero

export PureHemi, Hemiplex, μ, ν, mu, nu

struct PureHemi{T <: Real} <: Number
    m::T
    n::T
end
PureHemi(m::Real, n::Real) = PureHemi(promote(m, n)...)

struct Hemiplex{T <: Real} <: Number
    r::T
    h::PureHemi{T}
end

Hemiplex{T}(r::Real, m::Real, n::Real) where {T <: Real} = Hemiplex{T}(r, PureHemi{T}(m, n))
Hemiplex(r::T, m::T, n::T) where {T <: Real} = Hemiplex{T}(r, PureHemi{T}(m, n))
Hemiplex(r::Real, m::Real, n::Real) = Hemiplex(promote(r, m, n)...)

function Hemiplex(r::R, h::PureHemi{H}) where {R <: Real, H <: Real}
    T = promote_type(R, H)
    return Hemiplex{T}(r, PureHemi{T}(h.m, h.n))
end
Hemiplex(h::PureHemi{R}) where {R <: Real} = Hemiplex(zero(R), h)
Hemiplex(r::R) where {R <: Real} = Hemiplex(r, zero(PureHemi{R}))

const μ = PureHemi(true, false)
const ν = PureHemi(false, true)

## PureHemi implementation
convert(::Type{PureHemi{R}}, x::PureHemi) where {R} = PureHemi{R}(x.m, x.n)
convert(::Type{PureHemi{R}}, x::Real) where {R} =
    iszero(x) ? zero(PureHemi{R}) : throw(DomainError(x, "Non-zero reals cannot be converted to pure-hemi numbers"))

convert(::Type{R}, x::Hemiplex) where {R <: Real} =
    isreal(x) ? R(x.r) : throw(DomainError(x, "Hemiplex numbers with non-zero hemi-part cannot be converted to real numbers"))
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
(*)(x::PureHemi, y::PureHemi) = -(x.m * y.n + x.n * y.m) / 2
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
(/)(c::Real, x::PureHemi) = PureHemi(-c / x.n, -c / x.m)
(\)(x::PureHemi, c::Real) = c / x

real(::Type{PureHemi{R}}) where {R <: Real} = R
real(x::PureHemi{R}) where {R <: Real} = zero(R)
conj(x::PureHemi) = -x
mu(x::PureHemi) = x.m
nu(x::PureHemi) = x.n

zero(::Type{PureHemi{R}}) where {R <: Real} = PureHemi{R}(0, 0)
zero(::PureHemi{R}) where {R <: Real} = zero(PureHemi{R})

# (note that abs2(x) != x*x for PureHemi)
abs2(x::PureHemi) = x * conj(x)
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

## Hemiplex implementation
convert(::Type{Hemiplex{R}}, x::Hemiplex) where {R <: Real} = Hemiplex{R}(x.r, x.h)
convert(::Type{Hemiplex{R}}, x::PureHemi) where {R <: Real} = Hemiplex{R}(0, x)
convert(::Type{Hemiplex{R}}, x::Real) where {R <: Real} = Hemiplex{R}(x, zero(PureHemi{R}))
convert(::Type{PureHemi{R}}, x::Hemiplex) where {R <: Real} =
    iszero(x.r) ? convert(PureHemi{R}, x.h) : throw(DomainError(x, "Hemiplex numbers with non-zero real part cannot be converted to pure-hemi numbers"))

iszero(x::Hemiplex) = iszero(x.r) & iszero(x.h)
isnan(x::Hemiplex) = isnan(x.r) | isnan(x.h)
isinf(x::Hemiplex) = isinf(x.r) | isinf(x.h)
isfinite(x::Hemiplex) = isfinite(x.r) & isfinite(x.h)
isreal(x::Hemiplex) = iszero(x.h)
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

(+)(x::Hemiplex, y::Real) = x + Hemiplex(y)
(+)(x::Real, y::Hemiplex) = y + x
(-)(x::Hemiplex, y::Real) = x - Hemiplex(y)
(-)(x::Real, y::Hemiplex) = Hemiplex(x) - y
(*)(c::Bool, x::Hemiplex) = Hemiplex(c * x.r, c * x.h)
(*)(c::Real, x::Hemiplex) = Hemiplex(c * x.r, c * x.h)
(*)(x::Hemiplex, c::Real) = c * x
(/)(x::Hemiplex, c::Real) = Hemiplex(x.r / c, x.h / c)
(\)(c::Real, x::Hemiplex) = x / c

real(::Type{Hemiplex{R}}) where {R <: Real} = R
real(x::Hemiplex{R}) where {R <: Real} = x.r
conj(x::Hemiplex) = Hemiplex(x.r, -x.h)
mu(x::Hemiplex) = mu(x.h)
nu(x::Hemiplex) = nu(x.h)

zero(::Type{Hemiplex{R}}) where {R <: Real} = Hemiplex{R}(0, PureHemi{R}(0, 0))
zero(::Hemiplex{R}) where {R <: Real} = zero(Hemiplex{R})
one(::Type{Hemiplex{R}}) where {R <: Real} = Hemiplex{R}(1, PureHemi{R}(0, 0))
one(::Hemiplex{R}) where {R <: Real} = one(Hemiplex{R})

# (note that abs2(x) != x*x for the hemi part)
abs2(x::Hemiplex) = x.r * x.r + abs2(x.h)
abs(x::Hemiplex) = sqrt(abs2(x))

float(::Type{Hemiplex{T}}) where {T <: AbstractFloat} = Hemiplex{T}
float(::Type{Hemiplex{T}}) where {T} = Hemiplex{float(T)}
float(x::Hemiplex{<:AbstractFloat}) = x
float(x::Hemiplex) = Hemiplex(float(x.r), float(x.h))

widen(::Type{Hemiplex{T}}) where {T} = Hemiplex{widen(T)}

big(::Type{Hemiplex{T}}) where {T <: Real} = Hemiplex{big(T)}
big(x::Hemiplex{T}) where {T <: Real} = Hemiplex{big(T)}(x.r, x.h)

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
            out = Array{real(eltype(X))}(undef, size(X))
            for I in eachindex(X)
                out[I] = $f(X[I])
            end
            return out
        end
    end
end

promote_rule(::Type{Bool}, ::Type{PureHemi{H}}) where {H <: Real} = Hemiplex{promote_type(Bool, H)}
promote_rule(::Type{R}, ::Type{PureHemi{H}}) where {R <: AbstractIrrational, H <: Real} = Hemiplex{promote_type(R, H)}
promote_rule(::Type{R}, ::Type{PureHemi{H}}) where {R <: Real, H <: Real} = Hemiplex{promote_type(R, H)}
promote_rule(::Type{PureHemi{H1}}, ::Type{PureHemi{H2}}) where {H1 <: Real, H2 <: Real} = PureHemi{promote_type(H1, H2)}
promote_rule(::Type{Bool}, ::Type{Hemiplex{H}}) where {H <: Real} = Hemiplex{promote_type(Bool, H)}
promote_rule(::Type{R}, ::Type{Hemiplex{H}}) where {R <: AbstractIrrational, H <: Real} = Hemiplex{promote_type(R, H)}
promote_rule(::Type{R}, ::Type{Hemiplex{H}}) where {R <: Real, H <: Real} = Hemiplex{promote_type(R, H)}
promote_rule(::Type{PureHemi{H1}}, ::Type{Hemiplex{H2}}) where {H1 <: Real, H2 <: Real} = Hemiplex{promote_type(H1, H2)}
promote_rule(::Type{Hemiplex{H1}}, ::Type{Hemiplex{H2}}) where {H1 <: Real, H2 <: Real} = Hemiplex{promote_type(H1, H2)}

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

function show(io::IO, x::Hemiplex)
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

function read(s::IO, ::Type{PureHemi{T}}) where {T <: Real}
    m = read(s, T)
    n = read(s, T)
    PureHemi{T}(m, n)
end
function write(s::IO, x::PureHemi)
    write(s, x.m, x.n)
end

function read(s::IO, ::Type{Hemiplex{T}}) where {T <: Real}
    r = read(s, T)
    h = read(s, PureHemi{T})
    Hemiplex{T}(r, h)
end
function write(s::IO, x::Hemiplex)
    write(s, x.r, x.h)
end

end # module
