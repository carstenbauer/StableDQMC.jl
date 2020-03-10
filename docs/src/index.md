# Introduction

This is a library of numerically stable linear algebra routines for performing inversions as they appear in the calculation of Green's functions in [determinant Quantum Monte Carlo](https://en.wikipedia.org/wiki/Quantum_Monte_Carlo).

For more details, check out the accompanyig [paper](https://github.com/crstnbr/StableDQMC.jl/raw/master/paper/stabledqmc.pdf), in which we describe and benchmark a few specific algorithms. The plots in this paper have been generated with the notebooks in [the notebook directory](https://github.com/crstnbr/StableDQMC.jl/tree/master/notebooks) of this repository.

Feel free to give feedback, open issues, or contribute useful algorithms yourself! 🙂

## Installation

```
] add StableDQMC
```

The package has only one dependency, Requires.jl.

## Getting started

```julia
julia> using LinearAlgebra, StableDQMC

julia> B = rand(ComplexF64, 100,100);

julia> Bfact = udt(B);

julia> G = inv_one_plus_loh(Bfact);
```

Since the matrix `B` is well-conditioned in this case, we have

```julia
julia> G ≈ inv(I + B)
true
```

## Why should I care?

```@raw html
<table>
  <tr>
    <td><img src="https://github.com/crstnbr/StableDQMC.jl/raw/master/paper/figures/naive_vs_stable.png" width=500px></td>
    <td><img src="https://github.com/crstnbr/StableDQMC.jl/raw/master/paper/figures/decomp_comparison_simple.png" width=500px></td>
  </tr>
</table>
```

**Left:** Slice matrix product. **Right:** Equal-times Green's function.

Note that "SVD (D&C)" is the algorithm used by Julia's `svd` function.
