module StableDQMC

using SparseArrays, LinearAlgebra
using GenericSVD, JacobiSVD

include("linalg.jl")
include("qr.jl")
include("svd.jl")
include("Bchain.jl")


# Slice matrix chain B_M .... B_1
export calc_product_chain
export calc_product_chain_stabilized


# QR / UDT
export decompose_udt!
export decompose_udt
export multiply_safely
export UDT_to_mat!
export UDT_to_mat


export inv_udt!
export inv_udt
export inv_one_plus_udt!
export inv_one_plus_udt
export inv_one_plus_two_udts!
export inv_sum_udts!
export inv_sum_udts
export inv_sum_udts_loh!
export inv_sum_udts_loh
export inv_one_plus_udt_loh!



# SVD / UDV
export gesdd!
export gesdd
export gesvd!
export gesvd
export genericsvd!
export genericsvd
export gesvj!
export gesvj

export decompose_udv!
export decompose_udv
export UDV_to_mat!
export UDV_to_mat
export multiply_safely_udv

export inv_udv!
export inv_udv
export inv_one_plus_udv
export inv_one_plus_udv_alt
export inv_one_plus_udv_loh!
export inv_one_plus_udv_loh
export inv_sum_udvs!
export inv_sum_udvs





end # module