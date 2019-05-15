module StableDQMC

using SparseArrays, LinearAlgebra
using GenericSVD, JacobiSVD

include("linalg.jl")
include("qr.jl")
include("svd.jl")

export decompose_udv!,
    inv_one_plus_udv_alt,
    inv_sum_udvs!
    UDT_to_mat,
    inv_one_plus_udv_scalettar,
    inv_udt,
    UDT_to_mat!,
    inv_one_plus_udv_scalettar!,
    inv_udt!
    UDV_to_mat,
    inv_one_plus_two_udts!,
    inv_sum_udts,
    inv_udv,
    UDV_to_mat!,
    inv_one_plus_udt,
    inv_sum_udts!,
    inv_udv!,
    decompose_udt,
    inv_one_plus_udt!,
    inv_sum_udts_scalettar,
    multiply_safely,
    decompose_udt!,
    inv_one_plus_udt_scalettar!,
    inv_sum_udts_scalettar!,
    multiply_safely_udv,
    decompose_udv,
    inv_one_plus_udv,
    inv_sum_udvs,
    gesdd!,
    gesdd,
    gesvd!,
    gesvd,
    genericsvd!,
    genericsvd,
    gesvj!,
    gesvj

end # module
