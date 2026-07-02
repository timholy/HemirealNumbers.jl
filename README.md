# HemiplexNumbers

[![CI](https://github.com/timholy/HemiplexNumbers.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/timholy/HemiplexNumbers.jl/actions/workflows/CI.yml)
[![codecov](https://codecov.io/gh/timholy/HemiplexNumbers.jl/graph/badge.svg?token=AcT0zNxN71)](https://codecov.io/gh/timholy/HemiplexNumbers.jl)

A hemiplex number `z` can be written
```jl
z = r + mμ + nν
```
where `r`, `m`, and `n` are real, and the special numbers `μ`, `ν` satisfy
```jl
μ*μ = ν*ν = 0, μ*ν = ν*μ = -1/2.
```
Addition, subtraction, and any operation involving real numbers are
defined "the obvious way," and the conjugate of `z = r + mμ + nν` is `r - mμ - nν`.
Multiplication of general hemiplex numbers is commutative but not
associative.  Hemiplex numbers with `ν=0` are the same as dual numbers.

The motivation for inventing/rediscovering the hemiplex numbers was to
solve, using finite numbers, what would otherwise be singular
equations.
