<div align="center">
  <a href="https://crstnbr.github.io/StableDQMC.jl/dev">
    <img src="https://raw.githubusercontent.com/crstnbr/StableDQMC.jl/master/docs/src/assets/logo_large.png" alt="StableDQMC.jl" width="400">
  </a>
</div>

<h2 align="center">StableDQMC.jl</h2>
<p align="center">
  <span>Numerical stabilization routines for determinant quantum Monte Carlo</span>
</p>
<br>

> Nothing brings fear to my heart more than a floating point number. — Gerald Jay Sussman
| **Documentation**                                                               | **Build Status**                                                                                |  **Community**                                                                                |
|:-------------------------------------------------------------------------------:|:-----------------------------------------------------------------------------------------------:|:-----------------------------------------------------------------------------------------------:|
| [![][docs-dev-img]][docs-dev-url] | ![][lifecycle-img] [![][github-ci-img]][github-ci-url] [![][codecov-img]][codecov-url] [![][pkgeval-img]][pkgeval-url] | [![][slack-img]][slack-url] [![][license-img]][license-url] |

[docs-dev-img]: https://img.shields.io/badge/docs-dev-blue.svg
[docs-dev-url]: https://crstnbr.github.io/StableDQMC.jl/dev
[github-ci-img]: https://github.com/crstnbr/StableDQMC.jl/workflows/Run%20tests/badge.svg
[github-ci-url]: https://github.com/crstnbr/StableDQMC.jl/actions?query=workflow%3A%22Run+tests%22
[codecov-img]: https://img.shields.io/codecov/c/github/crstnbr/StableDQMC.jl/master.svg?label=codecov
[codecov-url]: http://codecov.io/github/crstnbr/StableDQMC.jl?branch=master

[pkgeval-img]: https://juliaci.github.io/NanosoldierReports/pkgeval_badges/S/StableDQMC.svg
[pkgeval-url]: https://juliaci.github.io/NanosoldierReports/pkgeval_badges/report.html

[slack-url]: https://slackinvite.julialang.org/
[slack-img]: https://img.shields.io/badge/chat-on%20slack-yellow.svg
[license-img]: https://img.shields.io/badge/License-MIT-red.svg
[license-url]: https://opensource.org/licenses/MIT

[lifecycle-img]: https://img.shields.io/badge/lifecycle-stable-blue.svg

This is a library of numerically stable linear algebra routines for performing inversions as they appear in the calculation of Green's functions in [determinant Quantum Monte Carlo](https://en.wikipedia.org/wiki/Quantum_Monte_Carlo).

For more details, check out the (unfinished!) accompanyig [notes](https://github.com/crstnbr/StableDQMC.jl/raw/master/paper/stabilization.pdf), in which we describe and benchmark a few specific algorithms. The plots in the notes have been generated with the notebooks in [the notebook directory](https://github.com/crstnbr/StableDQMC.jl/tree/master/notebooks) of this repository.

Feel free to give feedback, open issues, or contribute useful algorithms yourself! 🙂

### Installation

```
] add StableDQMC
```

The package has only one dependency, Requires.jl.

### Why should I care?

<table>
  <tr>
    <td><img src="paper/figures/naive_vs_stable.png" width=500px></td>
    <td><img src="paper/figures/decomp_comparison_simple.png" width=500px></td>
  </tr>
</table>

**Left:** Slice matrix product. **Right:** Equal-times Green's function.

Note that "SVD (D&C)" is the algorithm used by Julia's `svd` function.

### Inversions

In DQMC, we commonly perform inversions like `G = [1 + B]^-1` to obtain the equal-times Green's function and `G = [A + B]^-1` for the time-displaced pendant. The following methods are exported to facilitate these tasks.

- `inv_one_plus`, `inv_one_plus!`
- `inv_sum`, `inv_sum!`

When function names are suffixed with `_loh`, i.e. `inv_one_plus_loh`, a more sophisticated method is used for numerical stabilization (see the paper linked above for more details).

### Short example

```julia
julia> using LinearAlgebra

julia> B = rand(ComplexF64, 100,100);

julia> Bfact = udt(B);

julia> G = inv_one_plus_loh(Bfact);
```

Since the matrix `B` is well-conditioned in this case, we have

```julia
julia> G ≈ inv(I + B)
true
```

### UDT (QR) Decomposition

Based on the QR decomposition, we introduce a `UDT` factorization, where `U` is unitary, `D` is real-valued and diagonal, and `T` is upper-triangular. To decompose a given matrix `M` the `udt` function is exported.

```julia
julia> M = rand(10,10);

julia> udt(M)
UDT{Float64,Float64,Array{Float64,2}}([-0.246588 0.12668 … 0.582208 0.206435; -0.373953 -0.300804 … 0.152994 0.0523203; … ; -0.214686 -0.403362 … -0.124248 -0.390502; -0.40412 -0.147009 … 0.1839 0.197964], [2.15087, 1.47129, 1.14085, 0.911765, 0.850504, 0.620149, 0.545588, 0.412213, 0.305983, 0.148787], [-0.597235 -1.0 … -0.678767 -0.59054; -0.385741 0.0 … -1.0 -0.361263; … ; 0.0 0.0 … 0.0 0.0; 0.0 0.0 … 0.0 0.0])
```

In our tests (see `paper/`), this decomposition turns out to be superior to `SVD` for DQMC.

### SVD Decompositions

The package provides convenient access to several LAPACK algorithms for calculating singular value decompositions (SVDs):

* `gesdd`, `gesdd!`: Divide and conquer
* `gesvd`, `gesvd!`: Regular
* `gesvj`, `gesvj!`: Jacobi (based on [JacobiSVD.jl](https://github.com/RalphAS/JacobiSVD.jl))

Furthermore, you can access a type-generic, pure Julia implementation,

* `genericsvd`, `genericsvd!` (based on [GenericSVD.jl](https://github.com/JuliaLinearAlgebra/GenericSVD.jl))

However, to keep the dependencies of the package minimal, only the first two are available by default and loading of the Jacobi or type-generic SVD is opt-in. We provide convenience functions `StableDQMC.addJacobiSVD()` and `StableDQMC.addGenericSVD()` to facilitate this process. See below for a quick demonstration.

####  Automagically opt-in to JacobiSVD/GenericSVD

```julia
julia> using StableDQMC

julia> gesvj
ERROR: UndefVarError: gesvj not defined

julia> StableDQMC.addJacobiSVD()
  Updating registry at C:\Users\carsten\.julia\registries\General
  Updating git-repo https://github.com/JuliaRegistries/General.git
  Updating git-repo https://github.com/RalphAS/JacobiSVD.jl
 Resolving package versions...
  Updating C:\Users\carsten\Desktop\stabledqmctest\Project.toml
  [2ca068c6] + JacobiSVD v0.0.0 #master (https://github.com/RalphAS/JacobiSVD.jl)
  Updating C:\Users\carsten\Desktop\stabledqmctest\Manifest.toml
  [2ca068c6] + JacobiSVD v0.0.0 #master (https://github.com/RalphAS/JacobiSVD.jl)
┌ Warning: Package StableDQMC does not have JacobiSVD in its dependencies:
│ - If you have StableDQMC checked out for development and have
│   added JacobiSVD as a dependency but haven't updated your primary
│   environment's manifest file, try Pkg.resolve().
│ - Otherwise you may need to report an issue with StableDQMC
└ Loading JacobiSVD into StableDQMC from project dependency, future warnings for StableDQMC are suppressed.
[ Info: Recompiling stale cache file C:\Users\carsten.julia\compiled\v1.1\JacobiSVD\Frhox.ji for JacobiSVD [2ca068c6-2156-5cf0-8317-c67edf277a2c]

julia> gesvj # gesvj and gesvj! are available now
gesvj (generic function with 1 method)

```

#### Manual opt-in

Provided that you have [JacobiSVD.jl](https://github.com/RalphAS/JacobiSVD.jl) (or [GenericSVD.jl](https://github.com/JuliaLinearAlgebra/GenericSVD.jl)) installed, you can get the LAPACK access functions `gesvj`, `gesvj!` (or `genericsvd`, `genericsvd!`) simply by `import JacobiSVD` (`import GenericSVD`).

```julia
julia> using StableDQMC

julia> gesvj
ERROR: UndefVarError: gesvj not defined

julia> import JacobiSVD # using might lead to name conflicts

julia> gesvj # gesvj and gesvj! are available now
gesvj (generic function with 1 method)
```
