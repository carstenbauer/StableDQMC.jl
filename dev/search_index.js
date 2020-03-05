var documenterSearchIndex = {"docs":
[{"location":"decompositions/#Matrix-Decompositions-1","page":"Matrix Decompositions","title":"Matrix Decompositions","text":"","category":"section"},{"location":"decompositions/#UDT-(QR)-factorization-1","page":"Matrix Decompositions","title":"UDT (QR) factorization","text":"","category":"section"},{"location":"decompositions/#","page":"Matrix Decompositions","title":"Matrix Decompositions","text":"Based on the (pivoted) QR decomposition, we introduce a UDT factorization,","category":"page"},{"location":"decompositions/#","page":"Matrix Decompositions","title":"Matrix Decompositions","text":"beginaligned\n\tB = QR = UDT\nendaligned","category":"page"},{"location":"decompositions/#","page":"Matrix Decompositions","title":"Matrix Decompositions","text":"where we have split R into a diagonal piece D and upper triangular piece T. Hence, U = Q is unitary, D = textrmdiag(R) is a real diagonal matrix, and T is upper triangular.","category":"page"},{"location":"decompositions/#","page":"Matrix Decompositions","title":"Matrix Decompositions","text":"To decompose a given matrix M the udt function is exported.","category":"page"},{"location":"decompositions/#","page":"Matrix Decompositions","title":"Matrix Decompositions","text":"julia> M = rand(10,10);\n\njulia> udt(M)\nUDT{Float64,Float64,Array{Float64,2}}([-0.246588 0.12668 … 0.582208 0.206435; -0.373953 -0.300804 … 0.152994 0.0523203; … ; -0.214686 -0.403362 … -0.124248 -0.390502; -0.40412 -0.147009 … 0.1839 0.197964], [2.15087, 1.47129, 1.14085, 0.911765, 0.850504, 0.620149, 0.545588, 0.412213, 0.305983, 0.148787], [-0.597235 -1.0 … -0.678767 -0.59054; -0.385741 0.0 … -1.0 -0.361263; … ; 0.0 0.0 … 0.0 0.0; 0.0 0.0 … 0.0 0.0])","category":"page"},{"location":"decompositions/#","page":"Matrix Decompositions","title":"Matrix Decompositions","text":"In our tests (see paper/), this decomposition turns out to be superior to SVD for DQMC.","category":"page"},{"location":"decompositions/#SVD-factorization-1","page":"Matrix Decompositions","title":"SVD factorization","text":"","category":"section"},{"location":"decompositions/#","page":"Matrix Decompositions","title":"Matrix Decompositions","text":"A singular value decomposition (SVD) is given by","category":"page"},{"location":"decompositions/#","page":"Matrix Decompositions","title":"Matrix Decompositions","text":"beginaligned\n\tB = USV^dagger\nendaligned","category":"page"},{"location":"decompositions/#","page":"Matrix Decompositions","title":"Matrix Decompositions","text":"where U is unitary, S is a real diagonal matrix, and V^dagger is unitary.","category":"page"},{"location":"decompositions/#","page":"Matrix Decompositions","title":"Matrix Decompositions","text":"The package provides convenient access to several LAPACK algorithms for calculating singular value decompositions (SVDs):","category":"page"},{"location":"decompositions/#","page":"Matrix Decompositions","title":"Matrix Decompositions","text":"gesdd, gesdd!: Divide and conquer\ngesvd, gesvd!: Regular\ngesvj, gesvj!: Jacobi (based on JacobiSVD.jl)","category":"page"},{"location":"decompositions/#","page":"Matrix Decompositions","title":"Matrix Decompositions","text":"Furthermore, you can access a type-generic, pure Julia implementation,","category":"page"},{"location":"decompositions/#","page":"Matrix Decompositions","title":"Matrix Decompositions","text":"genericsvd, genericsvd! (based on GenericSVD.jl)","category":"page"},{"location":"decompositions/#","page":"Matrix Decompositions","title":"Matrix Decompositions","text":"However, to keep the dependencies of the package minimal, only the first two are available by default and loading of the Jacobi or type-generic SVD is opt-in. We provide convenience functions StableDQMC.addJacobiSVD() and StableDQMC.addGenericSVD() to facilitate this process. See below for a quick demonstration.","category":"page"},{"location":"decompositions/#Automagically-opt-in-to-JacobiSVD/GenericSVD-1","page":"Matrix Decompositions","title":"Automagically opt-in to JacobiSVD/GenericSVD","text":"","category":"section"},{"location":"decompositions/#","page":"Matrix Decompositions","title":"Matrix Decompositions","text":"julia> using StableDQMC\n\njulia> gesvj\nERROR: UndefVarError: gesvj not defined\n\njulia> StableDQMC.addJacobiSVD()\n  Updating registry at C:\\Users\\carsten\\.julia\\registries\\General\n  Updating git-repo https://github.com/JuliaRegistries/General.git\n  Updating git-repo https://github.com/RalphAS/JacobiSVD.jl\n Resolving package versions...\n  Updating C:\\Users\\carsten\\Desktop\\stabledqmctest\\Project.toml\n  [2ca068c6] + JacobiSVD v0.0.0 #master (https://github.com/RalphAS/JacobiSVD.jl)\n  Updating C:\\Users\\carsten\\Desktop\\stabledqmctest\\Manifest.toml\n  [2ca068c6] + JacobiSVD v0.0.0 #master (https://github.com/RalphAS/JacobiSVD.jl)\n┌ Warning: Package StableDQMC does not have JacobiSVD in its dependencies:\n│ - If you have StableDQMC checked out for development and have\n│   added JacobiSVD as a dependency but haven't updated your primary\n│   environment's manifest file, try Pkg.resolve().\n│ - Otherwise you may need to report an issue with StableDQMC\n└ Loading JacobiSVD into StableDQMC from project dependency, future warnings for StableDQMC are suppressed.\n[ Info: Recompiling stale cache file C:\\Users\\carsten.julia\\compiled\\v1.1\\JacobiSVD\\Frhox.ji for JacobiSVD [2ca068c6-2156-5cf0-8317-c67edf277a2c]\n\njulia> gesvj # gesvj and gesvj! are available now\ngesvj (generic function with 1 method)\n","category":"page"},{"location":"decompositions/#Manual-opt-in-1","page":"Matrix Decompositions","title":"Manual opt-in","text":"","category":"section"},{"location":"decompositions/#","page":"Matrix Decompositions","title":"Matrix Decompositions","text":"Provided that you have JacobiSVD.jl (or GenericSVD.jl) installed, you can get the LAPACK access functions gesvj, gesvj! (or genericsvd, genericsvd!) simply by import JacobiSVD (import GenericSVD).","category":"page"},{"location":"decompositions/#","page":"Matrix Decompositions","title":"Matrix Decompositions","text":"julia> using StableDQMC\n\njulia> gesvj\nERROR: UndefVarError: gesvj not defined\n\njulia> import JacobiSVD # using might lead to name conflicts\n\njulia> gesvj # gesvj and gesvj! are available now\ngesvj (generic function with 1 method)","category":"page"},{"location":"exports/#Methods-1","page":"Exports","title":"Methods","text":"","category":"section"},{"location":"exports/#Index-1","page":"Exports","title":"Index","text":"","category":"section"},{"location":"exports/#","page":"Exports","title":"Exports","text":"Pages = [\"exports.md\"]","category":"page"},{"location":"exports/#Documentation-1","page":"Exports","title":"Documentation","text":"","category":"section"},{"location":"exports/#","page":"Exports","title":"Exports","text":"Modules = [StableDQMC]\nPrivate = false\nOrder   = [:function, :type]","category":"page"},{"location":"exports/#StableDQMC.calc_Bchain-Tuple{Any,Any}","page":"Exports","title":"StableDQMC.calc_Bchain","text":"calc_Bchain(B,N) -> (R, svs)\n\nCalculate B^N as BBBB...*B from right to left.\n\nReturns the result R and the singular values (column of svs) of all intermediate powers (rows of svs).\n\n\n\n\n\n","category":"method"},{"location":"exports/#StableDQMC.calc_Bchain_qr-Tuple{Any,Any}","page":"Exports","title":"StableDQMC.calc_Bchain_qr","text":"Same as calc_Bchain but stabilizes the matrix products by intermediate matrix decompositions.\n\n\n\n\n\n","category":"method"},{"location":"exports/#StableDQMC.calc_Bchain_svd-Tuple{Any,Any}","page":"Exports","title":"StableDQMC.calc_Bchain_svd","text":"Same as calc_Bchain but stabilizes the matrix products by intermediate matrix decompositions.\n\n\n\n\n\n","category":"method"},{"location":"exports/#StableDQMC.calc_tdgf_qr-Tuple{Any,Any,Any}","page":"Exports","title":"StableDQMC.calc_tdgf_qr","text":"Calculate fake \"G(τ, 0)\", i.e. [B^-N1 + B^N2]^-1\n\n\n\n\n\n","category":"method"},{"location":"exports/#StableDQMC.calc_tdgf_svd-Tuple{Any,Any,Any}","page":"Exports","title":"StableDQMC.calc_tdgf_svd","text":"Calculate fake \"G(τ, 0)\", i.e. [B^-N1 + B^N2]^-1\n\n\n\n\n\n","category":"method"},{"location":"exports/#StableDQMC.fact_mult-Tuple{LinearAlgebra.SVD,LinearAlgebra.SVD}","page":"Exports","title":"StableDQMC.fact_mult","text":"fact_mult(A::SVD, B::SVD) -> SVD\n\nStabilized multiplication of two SVD decompositions. Returns a SVD factorization object.\n\n\n\n\n\n","category":"method"},{"location":"exports/#StableDQMC.fact_mult-Tuple{UDT,UDT}","page":"Exports","title":"StableDQMC.fact_mult","text":"fact_mult(A::UDT, B::UDT) -> UDT\n\nStabilized multiplication of two UDT decompositions. Returns a UDT factorization object.\n\n\n\n\n\n","category":"method"},{"location":"exports/#StableDQMC.gesdd!-Tuple{Any}","page":"Exports","title":"StableDQMC.gesdd!","text":"Same as gesdd but saves space by overwriting the input matrix.\n\n\n\n\n\n","category":"method"},{"location":"exports/#StableDQMC.gesdd-Tuple{Any}","page":"Exports","title":"StableDQMC.gesdd","text":"Calculates the SVD deomposition of a matrix by employing a divide and conquer approach. Returns a SVD factorization object.\n\n\n\n\n\n","category":"method"},{"location":"exports/#StableDQMC.gesvd!-Tuple{Any}","page":"Exports","title":"StableDQMC.gesvd!","text":"Same as gesvd but saves space by overwriting the input matrix.\n\n\n\n\n\n","category":"method"},{"location":"exports/#StableDQMC.gesvd-Tuple{Any}","page":"Exports","title":"StableDQMC.gesvd","text":"Calculates the SVD deomposition of a matrix. Returns a SVD factorization object.\n\n\n\n\n\n","category":"method"},{"location":"exports/#StableDQMC.inv!-Union{Tuple{M}, Tuple{Er}, Tuple{E}, Tuple{M,UDT{E,Er,M}}} where M where Er where E","page":"Exports","title":"StableDQMC.inv!","text":"inv!(res, F::UDT) -> res\n\nSame as inv but writes result into preallocated res.\n\n\n\n\n\n","category":"method"},{"location":"exports/#StableDQMC.inv_one_plus!-Tuple{Any,LinearAlgebra.SVD}","page":"Exports","title":"StableDQMC.inv_one_plus!","text":"inv_one_plus!(res, F::SVD) -> res\n\nStabilized calculation of [1 + USVt]^(-1):\n\nUse one intermediate SVD decomposition.\n\nWrites the result into res.\n\nSee svd_inv_one_plus for preallocation options.\n\n\n\n\n\n","category":"method"},{"location":"exports/#StableDQMC.inv_one_plus!-Tuple{Any,UDT,UDT}","page":"Exports","title":"StableDQMC.inv_one_plus!","text":"inv_one_plus!(res, A::UDT, Bdagger::UDT) -> res\n\nStabilized calculation of [1 + UlDlTl(UrDrTr)^†]^(-1). Writes the result into res.\n\nSee udt_inv_one_plus for preallocation options.\n\n\n\n\n\n","category":"method"},{"location":"exports/#StableDQMC.inv_one_plus!-Tuple{Any,UDT}","page":"Exports","title":"StableDQMC.inv_one_plus!","text":"inv_one_plus!(res, F::UDT) -> res\n\nSame as inv_one_plus but stores the result in preallocated res.\n\n\n\n\n\n","category":"method"},{"location":"exports/#StableDQMC.inv_one_plus-Tuple{LinearAlgebra.SVD}","page":"Exports","title":"StableDQMC.inv_one_plus","text":"inv_one_plus(F::SVD) -> AbstractMatrix\n\nStabilized calculation of [1 + USVt]^(-1):\n\nUse one intermediate SVD decomposition.\n\nSee svd_inv_one_plus for preallocation options.\n\n\n\n\n\n","category":"method"},{"location":"exports/#StableDQMC.inv_one_plus-Tuple{UDT,UDT}","page":"Exports","title":"StableDQMC.inv_one_plus","text":"inv_one_plus(A::UDT, Bdagger::UDT) -> AbstractMatrix\n\nStabilized calculation of [1 + UlDlTl(UrDrTr)^†]^(-1).\n\nSee udt_inv_one_plus for preallocation options.\n\n\n\n\n\n","category":"method"},{"location":"exports/#StableDQMC.inv_one_plus-Tuple{UDT}","page":"Exports","title":"StableDQMC.inv_one_plus","text":"inv_one_plus(F::UDT) -> AbstractMatrix\n\nStabilized calculation of [1 + UDT]^(-1):\n\nUse one intermediate UDT decomposition.\n\nFaster but potentially less accurate than inv_one_plus_loh.\n\nSee udt_inv_one_plus for preallocation options.\n\n\n\n\n\n","category":"method"},{"location":"exports/#StableDQMC.inv_one_plus_loh!-Tuple{Any,LinearAlgebra.SVD}","page":"Exports","title":"StableDQMC.inv_one_plus_loh!","text":"inv_one_plus_loh!(res, F::SVD) -> res\n\nStabilized calculation of [1 + USVt]^(-1):\n\nSeparate scales larger and smaller than unity\nUse two intermediate SVD decompositions.\n\nWrites the result into res.\n\nSee svd_inv_one_plus_loh for preallocation options.\n\n\n\n\n\n","category":"method"},{"location":"exports/#StableDQMC.inv_one_plus_loh!-Tuple{Any,UDT}","page":"Exports","title":"StableDQMC.inv_one_plus_loh!","text":"inv_one_plus_loh!(res, F::UDT) -> res\n\nStabilized calculation of [1 + UDT]^(-1):\n\nSeparate scales larger and smaller than unity\nUse two intermediate UDT decompositions.\n\nWrites the result into res.\n\nSee udt_inv_one_plus_loh for preallocation options.\n\n\n\n\n\n","category":"method"},{"location":"exports/#StableDQMC.inv_one_plus_loh-Tuple{LinearAlgebra.SVD}","page":"Exports","title":"StableDQMC.inv_one_plus_loh","text":"inv_one_plus_loh(F::SVD) -> AbstractMatrix\n\nStabilized calculation of [1 + USVt]^(-1):\n\nSeparate scales larger and smaller than unity\nUse two intermediate SVD decompositions.\n\nSee svd_inv_one_plus_loh for preallocation options.\n\n\n\n\n\n","category":"method"},{"location":"exports/#StableDQMC.inv_one_plus_loh-Tuple{UDT}","page":"Exports","title":"StableDQMC.inv_one_plus_loh","text":"inv_one_plus_loh(F::UDT) -> AbstractMatrix\n\nStabilized calculation of [1 + UDT]^(-1):\n\nSeparate scales larger and smaller than unity\nUse two intermediate UDT decompositions.\n\nSee udt_inv_one_plus_loh for preallocation options.\n\n\n\n\n\n","category":"method"},{"location":"exports/#StableDQMC.inv_sum!-Tuple{Any,LinearAlgebra.SVD,LinearAlgebra.SVD}","page":"Exports","title":"StableDQMC.inv_sum!","text":"inv_sum!(res, A::SVD, B::SVD) -> res\n\nStabilized calculation of [UaSaVta + UbSbVtb]^(-1):\n\nUse one intermediate SVD decompositions.\n\nWrites the result into res.\n\nSee svd_inv_sum for preallocation options.\n\n\n\n\n\n","category":"method"},{"location":"exports/#StableDQMC.inv_sum!-Tuple{Any,UDT,UDT}","page":"Exports","title":"StableDQMC.inv_sum!","text":"inv_sum!(res, A::UDT, B::UDT) -> res\n\nStabilized calculation of [UaDaTa + UbDbTb]^(-1):\n\nUse one intermediate UDT decompositions.\n\nSee udt_inv_sum for preallocation options.\n\n\n\n\n\n","category":"method"},{"location":"exports/#StableDQMC.inv_sum-Tuple{LinearAlgebra.SVD,LinearAlgebra.SVD}","page":"Exports","title":"StableDQMC.inv_sum","text":"inv_sum(A::SVD, B::SVD) -> res\n\nStabilized calculation of [UaSaVta + UbSbVtb]^(-1):\n\nUse one intermediate SVD decompositions.\n\nSee svd_inv_sum for preallocation options.\n\n\n\n\n\n","category":"method"},{"location":"exports/#StableDQMC.inv_sum-Tuple{UDT,UDT}","page":"Exports","title":"StableDQMC.inv_sum","text":"inv_sum(A::UDT, B::UDT) -> AbstractMatrix\n\nStabilized calculation of [UaDaTa + UbDbTb]^(-1):\n\nUse one intermediate UDT decompositions.\n\nSee udt_inv_sum for preallocation options.\n\n\n\n\n\n","category":"method"},{"location":"exports/#StableDQMC.inv_sum_loh!-Tuple{Any,LinearAlgebra.SVD,LinearAlgebra.SVD}","page":"Exports","title":"StableDQMC.inv_sum_loh!","text":"inv_sum_loh!(res, A::SVD, B::SVD) -> res\n\nStabilized calculation of [UaSaVta + UbSbVtb]^(-1):\n\nSeparate scales larger and smaller than unity\nUse two intermediate SVD decompositions.\n\nWrites the result into res.\n\nSee svd_inv_sum_loh for preallocation options.\n\n\n\n\n\n","category":"method"},{"location":"exports/#StableDQMC.inv_sum_loh!-Tuple{Any,UDT,UDT}","page":"Exports","title":"StableDQMC.inv_sum_loh!","text":"inv_sum_loh!(res, A::UDT, B::UDT) -> res\n\nStabilized calculation of [UaDaTa + UbDbTb]^(-1):\n\nSeparate scales larger and smaller than unity\nUse two intermediate UDT decompositions.\n\nWrites the result into res.\n\nSee udt_inv_sum_loh for preallocation options.\n\n\n\n\n\n","category":"method"},{"location":"exports/#StableDQMC.inv_sum_loh-Tuple{LinearAlgebra.SVD,LinearAlgebra.SVD}","page":"Exports","title":"StableDQMC.inv_sum_loh","text":"inv_sum_loh(A::SVD, B::SVD) -> res\n\nStabilized calculation of [UaSaVta + UbSbVtb]^(-1):\n\nSeparate scales larger and smaller than unity\nUse two intermediate SVD decompositions.\n\nSee svd_inv_sum_loh for preallocation options.\n\n\n\n\n\n","category":"method"},{"location":"exports/#StableDQMC.inv_sum_loh-Tuple{UDT,UDT}","page":"Exports","title":"StableDQMC.inv_sum_loh","text":"inv_sum_loh(A::UDT, B::UDT) -> AbstractMatrix\n\nStabilized calculation of [UaDaTa + UbDbTb]^(-1):\n\nSeparate scales larger and smaller than unity\nUse two intermediate UDT decompositions.\n\nSee udt_inv_sum_loh for preallocation options.\n\n\n\n\n\n","category":"method"},{"location":"exports/#StableDQMC.svd_inv_one_plus-Tuple{LinearAlgebra.SVD}","page":"Exports","title":"StableDQMC.svd_inv_one_plus","text":"svd_inv_one_plus(F::SVD) -> SVD\n\nStabilized calculation of [1 + USVt]^(-1):\n\nUse one intermediate SVD decomposition.\n\nOptions for preallocation via keyword arguments:\n\nt = similar(F.U)\na = similar(F.U)\nc = similar(F.U)\n\n\n\n\n\n","category":"method"},{"location":"exports/#StableDQMC.svd_inv_one_plus_loh-Tuple{LinearAlgebra.SVD}","page":"Exports","title":"StableDQMC.svd_inv_one_plus_loh","text":"svd_inv_one_plus_loh(F::SVD) -> SVD\n\nStabilized calculation of [1 + USVt]^(-1):\n\nSeparate scales larger and smaller than unity\nUse two intermediate SVD decompositions.\n\nOptions for preallocation via keyword arguments:\n\nl = similar(F.U)\nr = similar(F.U)\nDp = similar(F.D)\nDm = similar(F.D)\n\n\n\n\n\n","category":"method"},{"location":"exports/#StableDQMC.svd_inv_sum-Tuple{LinearAlgebra.SVD,LinearAlgebra.SVD}","page":"Exports","title":"StableDQMC.svd_inv_sum","text":"svd_inv_sum(A::SVD, B::SVD) -> SVD\n\nStabilized calculation of [UaSaVta + UbSbVtb]^(-1):\n\nUse one intermediate SVD decompositions.\n\nOptions for preallocations via keyword arguments:\n\nm1 = similar(A.U)\nm2 = similar(A.U)\n\n\n\n\n\n","category":"method"},{"location":"exports/#StableDQMC.svd_inv_sum_loh-Tuple{LinearAlgebra.SVD,LinearAlgebra.SVD}","page":"Exports","title":"StableDQMC.svd_inv_sum_loh","text":"svd_inv_sum_loh(A::SVD, B::SVD) -> SVD\n\nStabilized calculation of [UaSaVta + UbSbVtb]^(-1):\n\nSeparate scales larger and smaller than unity\nUse two intermediate SVD decompositions.\n\nOptions for preallocations via keyword arguments:\n\nNone (yet)\n\n\n\n\n\n","category":"method"},{"location":"exports/#StableDQMC.udt!-Union{Tuple{AbstractArray{C,2}}, Tuple{C}} where C<:Number","page":"Exports","title":"StableDQMC.udt!","text":"udv! is the same as svd, but saves space by overwriting the input A, instead of creating a copy.\n\n\n\n\n\n","category":"method"},{"location":"exports/#StableDQMC.udt-Union{Tuple{AbstractArray{C,2}}, Tuple{C}} where C<:Number","page":"Exports","title":"StableDQMC.udt","text":"Compute the UDT decomposition of A and return an UDT object.\n\nU, D, and T, can be obtained from the factorization F with F.U, F.D, and F.T such that A = U * Diagonal(D) * T.\n\nIterating the decomposition produces the components U, D, and V.\n\nNote that T is upper triangular only up to permutations of columns of T.\n\n\n\n\n\n","category":"method"},{"location":"exports/#StableDQMC.udt_inv_one_plus-Tuple{UDT,UDT}","page":"Exports","title":"StableDQMC.udt_inv_one_plus","text":"udt_inv_one_plus(A::UDT, Bdagger::UDT) -> UDT\n\nStabilized calculation of [1 + UlDlTl(UrDrTr)^†]^(-1). Returns and UDT factorization object.\n\nOptional preallocations via keyword arguments:\n\ntmp = similar(A.U)\ntmp2 = similar(A.U)\ntmp3 = similar(A.U)\n\n\n\n\n\n","category":"method"},{"location":"exports/#StableDQMC.udt_inv_one_plus-Tuple{UDT}","page":"Exports","title":"StableDQMC.udt_inv_one_plus","text":"udt_inv_one_plus(F::UDT) -> UDT\n\nStabilized calculation of [1 + UDT]^(-1). Returns and UDT factorization object.\n\nOptional preallocations via keyword arguments:\n\nu = similar(F.U)\nt = similar(F.T)\n\n\n\n\n\n","category":"method"},{"location":"exports/#StableDQMC.udt_inv_one_plus_loh-Tuple{UDT}","page":"Exports","title":"StableDQMC.udt_inv_one_plus_loh","text":"udt_inv_one_plus_loh(F::UDT) -> UDT\n\nStabilized calculation of [1 + UDT]^(-1):\n\nSeparate scales larger and smaller than unity\nUse two intermediate UDT decompositions.\n\nOptions for preallocation via keyword arguments:\n\nl = similar(F.U)\nr = similar(F.U)\nDp = similar(F.D)\nDm = similar(F.D)\n\n\n\n\n\n","category":"method"},{"location":"exports/#StableDQMC.udt_inv_sum-Tuple{UDT,UDT}","page":"Exports","title":"StableDQMC.udt_inv_sum","text":"udt_inv_sum(A::UDT, B::UDT) -> UDT\n\nStabilized calculation of [UaDaTa + UbDbTb]^(-1):\n\nUse one intermediate UDT decompositions.\n\nOptional preallocations via keyword arguments:\n\nm2 = similar(A.U)\n\n\n\n\n\n","category":"method"},{"location":"exports/#StableDQMC.udt_inv_sum_loh-Tuple{UDT,UDT}","page":"Exports","title":"StableDQMC.udt_inv_sum_loh","text":"udt_inv_sum_loh(A::UDT, B::UDT) -> UDT\n\nStabilized calculation of [UaDaTa + UbDbTb]^(-1):\n\nSeparate scales larger and smaller than unity\nUse two intermediate UDT decompositions.\n\nOptions for preallocations via keyword arguments:\n\nmat2 = similar(A.U)\nDap = similar(A.D)\nDam = similar(A.D)\nDbp = similar(A.D)\nDbm = similar(A.D)\n\n\n\n\n\n","category":"method"},{"location":"#Introduction-1","page":"Home","title":"Introduction","text":"","category":"section"},{"location":"#","page":"Home","title":"Home","text":"This is a library of numerically stable linear algebra routines for performing inversions as they appear in the calculation of Green's functions in determinant Quantum Monte Carlo.","category":"page"},{"location":"#","page":"Home","title":"Home","text":"For more details, check out the (unfinished!) accompanyig notes, in which we describe and benchmark a few specific algorithms. The plots in the notes have been generated with the notebooks in the notebook directory of this repository.","category":"page"},{"location":"#","page":"Home","title":"Home","text":"Feel free to give feedback, open issues, or contribute useful algorithms yourself! 🙂","category":"page"},{"location":"#Installation-1","page":"Home","title":"Installation","text":"","category":"section"},{"location":"#","page":"Home","title":"Home","text":"] add StableDQMC","category":"page"},{"location":"#","page":"Home","title":"Home","text":"The package has only one dependency, Requires.jl.","category":"page"},{"location":"#Getting-started-1","page":"Home","title":"Getting started","text":"","category":"section"},{"location":"#","page":"Home","title":"Home","text":"julia> using LinearAlgebra, StableDQMC\n\njulia> B = rand(ComplexF64, 100,100);\n\njulia> Bfact = udt(B);\n\njulia> G = inv_one_plus_loh(Bfact);","category":"page"},{"location":"#","page":"Home","title":"Home","text":"Since the matrix B is well-conditioned in this case, we have","category":"page"},{"location":"#","page":"Home","title":"Home","text":"julia> G ≈ inv(I + B)\ntrue","category":"page"},{"location":"#Why-should-I-care?-1","page":"Home","title":"Why should I care?","text":"","category":"section"},{"location":"#","page":"Home","title":"Home","text":"<table>\n  <tr>\n    <td><img src=\"https://github.com/crstnbr/StableDQMC.jl/raw/master/paper/figures/naive_vs_stable.png\" width=500px></td>\n    <td><img src=\"https://github.com/crstnbr/StableDQMC.jl/raw/master/paper/figures/decomp_comparison_simple.png\" width=500px></td>\n  </tr>\n</table>","category":"page"},{"location":"#","page":"Home","title":"Home","text":"Left: Slice matrix product. Right: Equal-times Green's function.","category":"page"},{"location":"#","page":"Home","title":"Home","text":"Note that \"SVD (D&C)\" is the algorithm used by Julia's svd function.","category":"page"},{"location":"inversions/#Stable-Inversions-1","page":"Stable Inversion","title":"Stable Inversions","text":"","category":"section"},{"location":"inversions/#Overview-1","page":"Stable Inversion","title":"Overview","text":"","category":"section"},{"location":"inversions/#","page":"Stable Inversion","title":"Stable Inversion","text":"In DQMC, to obtain the equal-times Green's function we commonly perform the inversion","category":"page"},{"location":"inversions/#","page":"Stable Inversion","title":"Stable Inversion","text":"G = left1 + Bright^-1","category":"page"},{"location":"inversions/#","page":"Stable Inversion","title":"Stable Inversion","text":"Similarily, we compute","category":"page"},{"location":"inversions/#","page":"Stable Inversion","title":"Stable Inversion","text":"G = leftA + Bright^-1","category":"page"},{"location":"inversions/#","page":"Stable Inversion","title":"Stable Inversion","text":"for the time-displaced Green's function. The following methods are exported to facilitate these tasks.","category":"page"},{"location":"inversions/#","page":"Stable Inversion","title":"Stable Inversion","text":"inv_one_plus, inv_one_plus!\ninv_sum, inv_sum!","category":"page"},{"location":"inversions/#","page":"Stable Inversion","title":"Stable Inversion","text":"When function names are suffixed with _loh, i.e. inv_one_plus_loh, a more sophisticated method is used for numerical stabilization (see the paper linked above for more details).","category":"page"},{"location":"inversions/#Details-1","page":"Stable Inversion","title":"Details","text":"","category":"section"},{"location":"inversions/#inv_one_plus-1","page":"Stable Inversion","title":"inv_one_plus","text":"","category":"section"},{"location":"inversions/#","page":"Stable Inversion","title":"Stable Inversion","text":"beginaligned\n\tG = mathbb1 + UDX^-1 \n\t= Uunderbrace(U^dagger X^-1 + D)_udxX^-1\n\t= (Uu)d(xX)^-1\n\t= U_r D_r X_r\nendaligned","category":"page"},{"location":"inversions/#","page":"Stable Inversion","title":"Stable Inversion","text":"with U_r = (xX)^-1, D_r = d^-1, X_r = (Uu)^-1.","category":"page"},{"location":"inversions/#inv_one_plus_loh-1","page":"Stable Inversion","title":"inv_one_plus_loh","text":"","category":"section"},{"location":"inversions/#","page":"Stable Inversion","title":"Stable Inversion","text":"beginaligned\n\tG = mathbb1 + UDX^-1\n\t= mathbb1 + UD_mD_pX^-1\n\t= (X^-1 D_p^-1 + U D_m) D_p X^-1\n\t= X^-1 underbraceD_p^-1 (underbraceX^-1 D_p^-1 + UD_m_udx)^-1_udx \n\t= U_r D_r X_r\nendaligned","category":"page"},{"location":"inversions/#","page":"Stable Inversion","title":"Stable Inversion","text":"with D_m = min(D 1), D_p = max(D 1), U_r = X^-1u, D_r = d, and X_r = x. [Loh2005, Loh1989]","category":"page"},{"location":"inversions/#inv_sum-1","page":"Stable Inversion","title":"inv_sum","text":"","category":"section"},{"location":"inversions/#","page":"Stable Inversion","title":"Stable Inversion","text":"beginaligned\n\tG(tau_1 tau_2) = U_L D_L X_L + U_R D_R X_R^-1\n\t= U_L underbrace(D_L X_L X_R^-1 + U_L^dagger U_R D_R)_udx X_R ^-1\n\t= (U_L u) d^-1 (x X_R)^-1\n\t= U_r D_r X_r\nendaligned","category":"page"},{"location":"inversions/#","page":"Stable Inversion","title":"Stable Inversion","text":"where U_r = (x X_R)^-1, D_r = d^-1, and X_r = (U_L u)^-1.","category":"page"},{"location":"inversions/#inv_sum_loh-1","page":"Stable Inversion","title":"inv_sum_loh","text":"","category":"section"},{"location":"inversions/#","page":"Stable Inversion","title":"Stable Inversion","text":"beginaligned\n\tG(tau_1 tau_2) = U_L D_L X_L + U_R D_R X_R^-1 \n\t= U_L D_Lm D_Lp X_L + U_R D_Rm D_Rp X_R^-1 \n\t= leftU_L D_Lp underbraceleft( dfracD_LmD_Rp X_L X_R^-1 + U_L^dagger U_R dfracD_RmD_Lp right)_udx X_R D_Rp right^-1 \n\t= X_R^-1 underbracedfrac1D_Rp udx^-1 dfrac1D_Lp_udx U_L^dagger\n\t= U_r D_r X_r\nendaligned","category":"page"},{"location":"inversions/#","page":"Stable Inversion","title":"Stable Inversion","text":"with D_Rm = min(D_R 1), D_Rp = max(D_R 1), D_Lm = min(D_L 1), D_Lp = max(D_L 1), U_r = X_R^-1 u, D_r = d, and X_r = x U_L^dagger. [Loh2005]","category":"page"},{"location":"inversions/#Resources-1","page":"Stable Inversion","title":"Resources","text":"","category":"section"},{"location":"inversions/#Loh2005-1","page":"Stable Inversion","title":"[Loh2005]","text":"","category":"section"},{"location":"inversions/#","page":"Stable Inversion","title":"Stable Inversion","text":"@article{Loh2005,\n\tauthor = {Loh, E. Y. and Gubernatis, J. E. and Scalettar, R. T. and White, S. R. and Scalapino, D. J. and Sugar, R. L.},\n\ttitle = {{Numerical Stability and the Sign Problem in the Determinant Quantum Monte Carlo Method}},\n\tjournal = {International Journal of Modern Physics C},\n\tvolume = {16},\n\tnumber = {08},\n\tpages = {1319--1327},\n\tmonth = {aug},\n\tyear = {2005}\n\tissn = {0129-1831},\n\tdoi = {10.1142/S0129183105007911},\n\turl = {http://www.worldscientific.com/doi/abs/10.1142/S0129183105007911},\n}","category":"page"},{"location":"inversions/#Loh1989-1","page":"Stable Inversion","title":"[Loh1989]","text":"","category":"section"},{"location":"inversions/#","page":"Stable Inversion","title":"Stable Inversion","text":"@incollection{Loh1989,\n\tauthor = {Loh, E. Y. and Gubernatis, J. E. and Scalettar, R. T. and Sugar, R. L. and White, S. R.},\n\ttitle = {{Stable Matrix-Multiplication Algorithms for Low-Temperature Numerical Simulations of Fermions}},\n\tyear = {1989}\n\tpages = {55--60},\n\tdoi = {10.1007/978-1-4613-0565-1_8},\n\turl = {http://link.springer.com/10.1007/978-1-4613-0565-1{\\_}8},\n}","category":"page"}]
}
