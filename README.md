# HemiplexNumbers

[![CI](https://github.com/timholy/HemiplexNumbers.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/timholy/HemiplexNumbers.jl/actions/workflows/CI.yml)
[![codecov](https://codecov.io/gh/timholy/HemiplexNumbers.jl/graph/badge.svg?token=AcT0zNxN71)](https://codecov.io/gh/timholy/HemiplexNumbers.jl)

A hemiplex number `z` can be written
```jl
z = r + mμ + nν
```
where `r`, `m`, and `n` are real or complex and the special numbers `μ`, `ν` satisfy
```jl
μ*μ = ν*ν = 0, μ*ν = ν*μ = -1/2.
```
Addition, subtraction, and scalar multiplication act coefficientwise.
Multiplication is commutative but, in general, not associative.  Hemiplex
numbers with `ν=0` are dual numbers.

Conjugation combines hemi-conjugation with conjugation of the coefficients:
`conj(z) = conj(r) - conj(m)μ - conj(n)ν`.  The function `coefftype(z)`
returns the coefficient type, while `real(typeof(z))` returns its underlying
real type.  For complex-coefficient `Hemiplex` numbers, `z*conj(z)` generally
retains a hemi-part, so `abs2` and `abs` are supported only for real
coefficients.  (`abs2` of a `PureHemi` is real-valued with either coefficient
type.)  Inversion works with both: `inv` uses the coefficient-linear
hemi-conjugate `r - mμ - nν`, whose product with `z` is the scalar `r² + mn`.

The motivation for inventing/rediscovering the hemiplex numbers was to
solve, using finite numbers, what would otherwise be singular
equations.  They are described in [Edelman & Holy,
arXiv:2607.21383](https://arxiv.org/abs/2607.21383);
[HemiplexFactorizations.jl](https://github.com/timholy/HemiplexFactorizations.jl)
applies them to the Cholesky decomposition of arbitrary symmetric
matrices.

## Citation

If you use this package in published work, please cite:

> A. Edelman and T. E. Holy, "Jordan algebras, hemiplex numbers, and the
> Cholesky decomposition of arbitrary symmetric matrices," arXiv:2607.21383
> (2026). https://doi.org/10.48550/arXiv.2607.21383

```bibtex
@article{EdelmanHoly2026,
  author  = {Edelman, Alan and Holy, Timothy E.},
  title   = {Jordan algebras, hemiplex numbers, and the {Cholesky}
             decomposition of arbitrary symmetric matrices},
  journal = {arXiv},
  year    = {2026},
  doi     = {10.48550/arXiv.2607.21383},
  url     = {https://arxiv.org/abs/2607.21383},
}
```

Machine-readable metadata is in [`CITATION.cff`](CITATION.cff).
