# HemirealNumbers

[![Build Status](https://travis-ci.org/timholy/HemirealNumbers.jl.svg?branch=master)](https://travis-ci.org/timholy/HemirealNumbers.jl)
[![codecov.io](https://codecov.io/github/timholy/HemirealNumbers.jl/coverage.svg?branch=master)](https://codecov.io/github/timholy/HemirealNumbers.jl?branch=master)

A hemireal number `z` can be written
```jl
z = r + mμ + nν
```
where `r`, `m`, and `n` are real, and the special numbers `μ`, `ν` satisfy
```jl
μ*μ = ν*ν = 0, μ*ν = ν*μ = 1.
```
Addition, subtraction, and any operation involving real numbers are
defined "the obvious way," and the conjugate of `z` is just `z`.
Multiplication of general hemireals is commutative but not
associative.  Hemireals with `ν=0` are the same as dual numbers.

The motivation for inventing/rediscovering (?) the hemireals was to
solve, using finite numbers, what would otherwise be singular
equations.
