module HemirealNumbers

import Base: +, -, *, /, \, ^, abs, abs2, conj, convert, isfinite, promote_op, promote_rule, real, show, zero

export PureHemi, Hemireal, μ, ν, mu, nu

struct PureHemi{T<:Real} <: Number
    m::T
    n::T
end
PureHemi(m::Real, n::Real) = PureHemi(promote(m,n)...)

struct Hemireal{T<:Real} <: Number
    r::T
    h::PureHemi{T}
end

Hemireal{T}(r::Real, m::Real, n::Real) where {T<:Real} = Hemireal{T}(r, PureHemi{T}(m, n))
Hemireal(r::T, m::T, n::T) where {T<:Real} = Hemireal{T}(r, PureHemi{T}(m, n))
Hemireal(r::Real, m::Real, n::Real) = Hemireal(promote(r, m, n)...)

function Hemireal(r::R, h::PureHemi{H}) where {R<:Real,H<:Real}
    T = promote_type(R,H)
    Hemireal{T}(r, PureHemi{T}(h.m, h.n))
end
Hemireal(h::PureHemi{R}) where {R<:Real} = Hemireal(zero(R), h)
Hemireal(r::R) where {R<:Real} = Hemireal(r, zero(PureHemi{R}))

const μ = PureHemi(true,false)
const ν = PureHemi(false,true)

## PureHemi implementation
convert(::Type{PureHemi{R}}, x::PureHemi) where R = PureHemi{R}(x.m, x.n)
convert(::Type{PureHemi{R}}, x::Real) where R = iszero(x) ? zero(PureHemi{R}) : throw(DomainError(x, "Non-zero reals cannot be converted to pure-hemi numbers"))

(-)(x::PureHemi) = PureHemi(-x.m, -x.n)

(+)(x::PureHemi, y::PureHemi) = PureHemi(x.m+y.m, x.n+y.n)
(-)(x::PureHemi, y::PureHemi) = PureHemi(x.m-y.m, x.n-y.n)
(*)(x::PureHemi, y::PureHemi) = x.m*y.n + x.n*y.m
(*)(c::Bool, x::PureHemi)     = PureHemi(c*x.m, c*x.n)
(*)(c::Real, x::PureHemi)     = PureHemi(c*x.m, c*x.n)
(*)(x::PureHemi, c::Real)     = c*x
(/)(x::PureHemi, c::Real)     = PureHemi(x.m/c, x.n/c)
(\)(c::Real, x::PureHemi)     = x/c
function (/)(x::PureHemi, y::PureHemi)
    @assert x.m*y.n == x.n*y.m
    y.n == 0 ? x.m/y.m : x.n/y.n
end
(\)(x::PureHemi, y::PureHemi) = y/x
(^)(x::PureHemi, p::Integer)  = (x*x)^(p/2)  # note x*x is much faster than x^2
(^)(x::PureHemi, p::Rational) = (x*x)^(p/2)
(^)(x::PureHemi, p::Real)     = (x*x)^(p/2)

# Symmetric division (there are other solutions y to y*x = c)
(/)(c::Real, x::PureHemi) = PureHemi(c/(2*x.n), c/(2*x.m))
(\)(x::PureHemi, c::Real) = c/x

real(::Type{PureHemi{R}}) where {R<:Real} = R
real(x::PureHemi{R}) where {R<:Real} = zero(R)
conj(x::PureHemi) = x
mu(x::PureHemi) = x.m
nu(x::PureHemi) = x.n

zero(::Type{PureHemi{R}}) where {R<:Real} = PureHemi{R}(0, 0)
zero(::PureHemi{R}) where {R<:Real} = zero(PureHemi{R})

isfinite(x::PureHemi) = isfinite(x.m) && isfinite(x.n)
# I'm not entirely certain whether the next two should even be defined.
# (note that abs2(x) != x*x)
abs2(x::PureHemi) = x.m*x.m + x.n*x.n
abs(x::PureHemi) = sqrt(abs2(x))

## Hemireal implementation
convert(::Type{Hemireal{R}}, x::Hemireal) where {R<:Real} = Hemireal{R}(x.r, x.h)
convert(::Type{Hemireal{R}}, x::PureHemi) where {R<:Real} = Hemireal{R}(0, x)
convert(::Type{Hemireal{R}}, x::Real) where {R<:Real} = Hemireal{R}(x, zero(PureHemi{R}))
convert(::Type{Hemireal{R}}, r::Real, m::Real, n::Real) where {R<:Real} = Hemireal{R}(r, PureHemi{R}(m, n))

(-)(x::Hemireal) = Hemireal(-x.r, -x.h)

(+)(x::Hemireal, y::Hemireal) = Hemireal(x.r+y.r, x.h+y.h)
(-)(x::Hemireal, y::Hemireal) = Hemireal(x.r-y.r, x.h-y.h)
(*)(x::Hemireal, y::Hemireal) = Hemireal(x.r*y.r + x.h*y.h, x.r*y.h + x.h*y.r)

(+)(x::Hemireal, y::PureHemi) = x+Hemireal(y)
(+)(x::PureHemi, y::Hemireal) = y+x
(-)(x::Hemireal, y::PureHemi) = x-Hemireal(y)
(-)(x::PureHemi, y::Hemireal) = Hemireal(x)-y
(*)(x::Hemireal, y::PureHemi) = x*Hemireal(y)
(*)(x::PureHemi, y::Hemireal) = y*x

(+)(x::Hemireal, y::Real) = x+Hemireal(y)
(+)(x::Real, y::Hemireal) = y+x
(-)(x::Hemireal, y::Real) = x-Hemireal(y)
(-)(x::Real, y::Hemireal) = Hemireal(x)-y
(*)(c::Bool, x::Hemireal) = Hemireal(c*x.r, c*x.h)
(*)(c::Real, x::Hemireal) = Hemireal(c*x.r, c*x.h)
(*)(x::Hemireal, c::Real) = c*x
(/)(x::Hemireal, c::Real) = Hemireal(x.r/c, x.h/c)
(\)(c::Real, x::Hemireal) = x/c

real(::Type{Hemireal{R}}) where {R<:Real} = R
real(x::Hemireal{R}) where {R<:Real} = x.r
conj(x::Hemireal) = x
mu(x::Hemireal) = mu(x.h)
nu(x::Hemireal) = nu(x.h)

zero(::Type{Hemireal{R}}) where {R<:Real} = Hemireal{R}(0, PureHemi{R}(0, 0))
zero(::Hemireal{R}) where {R<:Real} = zero(Hemireal{R})

isfinite(x::Hemireal) = isfinite(x.r) && isfinite(x.h)
# ?
abs2(x::Hemireal) = x.r*x.r + abs2(x.h)
abs(x::Hemireal) = sqrt(abs2(x))

for f in (:mu, :nu)
    @eval begin
        function $f(X::AbstractArray)
            out = Array{real(eltype(X))}(undef, size(X))
            for I in eachindex(X)
                out[I] = $f(X[I])
            end
            out
        end
    end
end

# (*){H1<:PureHemi,H2<:PureHemi}(A::StridedVecOrMat{H1}, B::StridedVector{H2}) = A_mul_B!(Array(promote_type(real(H1),real(H2)), size(A,1)), A, B)
# (*){H1<:PureHemi,H2<:PureHemi}(A::StridedVecOrMat{H1}, B::StridedMatrix{H2}) = A_mul_B!(Array(promote_type(real(H1),real(H2)), size(A,1), size(B,2)), A, B)
# A_mul_Bt{H1<:PureHemi,H2<:PureHemi}(A::StridedVecOrMat{H1}, B::StridedVecOrMat{H2}) = A_mul_Bt!(Array(promote_type(real(H1),real(H2)), size(A,1), size(B,1)), A, B)
# A_mul_Bc{H1<:PureHemi,H2<:PureHemi}(A::StridedVecOrMat{H1}, B::StridedVecOrMat{H2}) = A_mul_Bt!(Array(promote_type(real(H1),real(H2)), size(A,1), size(B,1)), A, B)
# At_mul_B{H1<:PureHemi,H2<:PureHemi}(A::StridedVecOrMat{H1}, B::StridedVector{H2}) = At_mul_B!(Array(promote_type(real(H1),real(H2)), size(A,2)), A, B)
# At_mul_B{H1<:PureHemi,H2<:PureHemi}(A::StridedVecOrMat{H1}, B::StridedMatrix{H2}) = At_mul_B!(Array(promote_type(real(H1),real(H2)), size(A,2), size(B,2)), A, B)
# Ac_mul_B{H1<:PureHemi,H2<:PureHemi}(A::StridedVecOrMat{H1}, B::StridedVector{H2}) = At_mul_B!(Array(promote_type(real(H1),real(H2)), size(A,2)), A, B)
# Ac_mul_B{H1<:PureHemi,H2<:PureHemi}(A::StridedVecOrMat{H1}, B::StridedMatrix{H2}) = At_mul_B!(Array(promote_type(real(H1),real(H2)), size(A,2), size(B,2)), A, B)
# At_mul_Bt{H1<:PureHemi,H2<:PureHemi}(A::StridedVecOrMat{H1}, B::StridedVecOrMat{H2}) = At_mul_Bt!(Array(promote_type(real(H1),real(H2)), size(A,2), size(B,1)), A, B)
# Ac_mul_Bc{H1<:PureHemi,H2<:PureHemi}(A::StridedVecOrMat{H1}, B::StridedVecOrMat{H2}) = At_mul_Bt!(Array(promote_type(real(H1),real(H2)), size(A,2), size(B,1)), A, B)

# promote_op(::Base.MulFun, ::Type{PureHemi{H1}}, ::Type{PureHemi{H2}})  where {H1<:Real,H2<:Real} = promote_type(H1,H2)
# promote_op(::Base.DotMulFun, ::Type{PureHemi{H1}}, ::Type{PureHemi{H2}}) where {H1<:Real,H2<:Real}  = promote_type(H1,H2)
# promote_op(::Base.MulFun, ::Type{R}, ::Type{PureHemi{H}}) where {R<:Real,H} = PureHemi{promote_type(R,H)}
# promote_op(::Base.DotMulFun, ::Type{R}, ::Type{PureHemi{H}}) where {R<:Real,H} = PureHemi{promote_type(R,H)}
# promote_op(::Base.MulFun, ::Type{PureHemi{H}}, ::Type{R}) where {R<:Real,H} = promote_op(Base.MulFun(), R, PureHemi{H})
# promote_op(::Base.DotMulFun, ::Type{PureHemi{H}}, ::Type{R}) where {R<:Real,H} = promote_op(Base.MulFun(), R, PureHemi{H})

promote_rule(::Type{Bool}, ::Type{PureHemi{H}}) where {H<:Real} = Hemireal{promote_type(Bool,H)}
promote_rule(::Type{R}, ::Type{PureHemi{H}}) where {R<:AbstractIrrational,H<:Real} = Hemireal{promote_type(R,H)}
promote_rule(::Type{R}, ::Type{PureHemi{H}}) where {R<:Real,H<:Real} = Hemireal{promote_type(R,H)}
promote_rule(::Type{PureHemi{H1}}, ::Type{PureHemi{H2}}) where {H1<:Real,H2<:Real} = PureHemi{promote_type(H1,H2)}
promote_rule(::Type{Bool}, ::Type{Hemireal{H}}) where {H<:Real} = Hemireal{promote_type(Bool,H)}
promote_rule(::Type{R}, ::Type{Hemireal{H}}) where {R<:AbstractIrrational,H<:Real} = Hemireal{promote_type(R,H)}
promote_rule(::Type{R}, ::Type{Hemireal{H}}) where {R<:Real,H<:Real} = Hemireal{promote_type(R,H)}
promote_rule(::Type{PureHemi{H1}}, ::Type{Hemireal{H2}}) where {H1<:Real,H2<:Real} = Hemireal{promote_type(H1,H2)}
promote_rule(::Type{Hemireal{H1}}, ::Type{Hemireal{H2}}) where {H1<:Real,H2<:Real} = Hemireal{promote_type(H1,H2)}

show(io::IO, x::PureHemi) = x.n >= 0 ? print(io, x.m, "μ + ", x.n, 'ν') : print(io, x.m, "μ - ", abs(x.n), 'ν')
show(io::IO, x::Hemireal) = x.h.m >= 0 ? print(io, x.r, " + ", x.h) : print(io, x.r, " - ", PureHemi(-x.h.m, x.h.n))

end # module
