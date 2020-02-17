# Matrix Decompositions

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
