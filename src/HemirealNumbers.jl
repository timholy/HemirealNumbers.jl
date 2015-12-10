__precompile__()

module HemirealNumbers

import Base: +, -, *, /, \, ^, abs, abs2, conj, convert, isfinite, promote_op, promote_rule, real, show, zero
import Base.LinAlg: A_mul_Bt, A_mul_Bc, At_mul_B, Ac_mul_B, At_mul_Bt, Ac_mul_Bc

export PureHemi, Hemireal, μ, ν, mu, nu

immutable PureHemi{T<:Real} <: Number
    m::T
    n::T
end
PureHemi(m::Real, n::Real) = PureHemi(promote(m,n)...)

immutable Hemireal{T<:Real} <: Number
    r::T
    h::PureHemi{T}
end

Hemireal(r::Real, m::Real, n::Real) = Hemireal(r, PureHemi(m, n))
function Hemireal{R,H}(r::R, h::PureHemi{H})
    T = promote_type(R,H)
    Hemireal{T}(r, PureHemi{T}(h.m, h.n))
end
Hemireal{R}(h::PureHemi{R}) = Hemireal(zero(R), h)
Hemireal{R}(r::R) = Hemireal(r, zero(PureHemi{R}))

const μ = PureHemi(true,false)
const ν = PureHemi(false,true)

## PureHemi implementation
convert{R}(::Type{PureHemi{R}}, x::PureHemi) = PureHemi{R}(x.m, x.n)
convert{R}(::Type{PureHemi{R}}, x::Real) = x == 0 ? zero(PureHemi{R}) : throw(DomainError())  # error("Non-zero reals cannot be converted to pure-hemi numbers")

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

real{R}(::Type{PureHemi{R}}) = R
real{R}(x::PureHemi{R}) = zero(R)
conj(x::PureHemi) = x
mu(x::PureHemi) = x.m
nu(x::PureHemi) = x.n

zero{R}(::Type{PureHemi{R}}) = PureHemi{R}(0, 0)
zero{R}(::PureHemi{R}) = zero(PureHemi{R})

isfinite(x::PureHemi) = isfinite(x.m) && isfinite(x.n)
# I'm not entirely certain whether the next two should even be defined.
# (note that abs2(x) != x*x)
abs2(x::PureHemi) = x.m*x.m + x.n*x.n
abs(x::PureHemi) = sqrt(abs2(x))

## Hemireal implementation
convert{R}(::Type{Hemireal{R}}, x::Hemireal) = Hemireal{R}(x.r, x.h)
convert{R}(::Type{Hemireal{R}}, x::PureHemi) = Hemireal{R}(0, x)
convert{R}(::Type{Hemireal{R}}, x::Real) = Hemireal{R}(x, zero(PureHemi{R}))
convert{R}(::Type{Hemireal{R}}, r::Real, m::Real, n::Real) = Hemireal{R}(r, PureHemi{R}(m, n))

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

real{R}(::Type{Hemireal{R}}) = R
real{R}(x::Hemireal{R}) = x.r
conj(x::Hemireal) = x
mu(x::Hemireal) = mu(x.h)
nu(x::Hemireal) = nu(x.h)

zero{R}(::Type{Hemireal{R}}) = Hemireal{R}(0, PureHemi{R}(0, 0))
zero{R}(::Hemireal{R}) = zero(Hemireal{R})

isfinite(x::Hemireal) = isfinite(x.r) && isfinite(x.h)
# ?
abs2(x::Hemireal) = x.r*x.r + abs2(x.h)
abs(x::Hemireal) = sqrt(abs2(x))

for f in (:mu, :nu)
    @eval begin
        function $f(X::AbstractArray)
            out = Array(real(eltype(X)), size(X))
            for I in eachindex(X)
                out[I] = $f(X[I])
            end
            out
        end
    end
end

(*){H1<:PureHemi,H2<:PureHemi}(A::StridedVecOrMat{H1}, B::StridedVector{H2}) = A_mul_B!(Array(promote_type(real(H1),real(H2)), size(A,1)), A, B)
(*){H1<:PureHemi,H2<:PureHemi}(A::StridedVecOrMat{H1}, B::StridedMatrix{H2}) = A_mul_B!(Array(promote_type(real(H1),real(H2)), size(A,1), size(B,2)), A, B)
A_mul_Bt{H1<:PureHemi,H2<:PureHemi}(A::StridedVecOrMat{H1}, B::StridedVecOrMat{H2}) = A_mul_Bt!(Array(promote_type(real(H1),real(H2)), size(A,1), size(B,1)), A, B)
A_mul_Bc{H1<:PureHemi,H2<:PureHemi}(A::StridedVecOrMat{H1}, B::StridedVecOrMat{H2}) = A_mul_Bt!(Array(promote_type(real(H1),real(H2)), size(A,1), size(B,1)), A, B)
At_mul_B{H1<:PureHemi,H2<:PureHemi}(A::StridedVecOrMat{H1}, B::StridedVector{H2}) = At_mul_B!(Array(promote_type(real(H1),real(H2)), size(A,2)), A, B)
At_mul_B{H1<:PureHemi,H2<:PureHemi}(A::StridedVecOrMat{H1}, B::StridedMatrix{H2}) = At_mul_B!(Array(promote_type(real(H1),real(H2)), size(A,2), size(B,2)), A, B)
Ac_mul_B{H1<:PureHemi,H2<:PureHemi}(A::StridedVecOrMat{H1}, B::StridedVector{H2}) = At_mul_B!(Array(promote_type(real(H1),real(H2)), size(A,2)), A, B)
Ac_mul_B{H1<:PureHemi,H2<:PureHemi}(A::StridedVecOrMat{H1}, B::StridedMatrix{H2}) = At_mul_B!(Array(promote_type(real(H1),real(H2)), size(A,2), size(B,2)), A, B)
At_mul_Bt{H1<:PureHemi,H2<:PureHemi}(A::StridedVecOrMat{H1}, B::StridedVecOrMat{H2}) = At_mul_Bt!(Array(promote_type(real(H1),real(H2)), size(A,2), size(B,1)), A, B)
Ac_mul_Bc{H1<:PureHemi,H2<:PureHemi}(A::StridedVecOrMat{H1}, B::StridedVecOrMat{H2}) = At_mul_Bt!(Array(promote_type(real(H1),real(H2)), size(A,2), size(B,1)), A, B)

promote_op{H1,H2}(::Base.MulFun, ::Type{PureHemi{H1}}, ::Type{PureHemi{H2}}) = promote_type(H1,H2)
promote_op{H1,H2}(::Base.DotMulFun, ::Type{PureHemi{H1}}, ::Type{PureHemi{H2}}) = promote_type(H1,H2)
promote_op{R<:Real,H}(::Base.MulFun, ::Type{R}, ::Type{PureHemi{H}}) = PureHemi{promote_type(R,H)}
promote_op{R<:Real,H}(::Base.DotMulFun, ::Type{R}, ::Type{PureHemi{H}}) = PureHemi{promote_type(R,H)}
promote_op{R<:Real,H}(::Base.MulFun, ::Type{PureHemi{H}}, ::Type{R}) = promote_op(Base.MulFun(), R, PureHemi{H})
promote_op{R<:Real,H}(::Base.DotMulFun, ::Type{PureHemi{H}}, ::Type{R}) = promote_op(Base.MulFun(), R, PureHemi{H})

promote_rule{H}(::Type{Bool}, ::Type{PureHemi{H}}) = Hemireal{promote_type(Bool,H)}
promote_rule{R,H}(::Type{Irrational{R}}, ::Type{PureHemi{H}}) = Hemireal{promote_type(Irrational{R},H)}
promote_rule{R<:Real,H}(::Type{R}, ::Type{PureHemi{H}}) = Hemireal{promote_type(R,H)}
promote_rule{H1,H2}(::Type{PureHemi{H1}}, ::Type{PureHemi{H2}}) = PureHemi{promote_type(H1,H2)}
promote_rule{H}(::Type{Bool}, ::Type{Hemireal{H}}) = Hemireal{promote_type(Bool,H)}
promote_rule{R,H}(::Type{Irrational{R}}, ::Type{Hemireal{H}}) = Hemireal{promote_type(Irrational{R},H)}
promote_rule{R<:Real,H}(::Type{R}, ::Type{Hemireal{H}}) = Hemireal{promote_type(R,H)}
promote_rule{H1,H2}(::Type{PureHemi{H1}}, ::Type{Hemireal{H2}}) = Hemireal{promote_type(H1,H2)}
promote_rule{H1,H2}(::Type{Hemireal{H1}}, ::Type{Hemireal{H2}}) = Hemireal{promote_type(H1,H2)}

show(io::IO, x::PureHemi) = x.n >= 0 ? print(io, x.m, "μ + ", x.n, 'ν') : print(io, x.m, "μ - ", abs(x.n), 'ν')
show(io::IO, x::Hemireal) = x.h.m >= 0 ? print(io, x.r, " + ", x.h) : print(io, x.r, " - ", PureHemi(-x.h.m, x.h.n))

end # module
