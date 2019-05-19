module StableDQMC

using SparseArrays, LinearAlgebra
using GenericSVD, JacobiSVD

include("linalg.jl")
include("udt_type.jl")
include("qr.jl")
include("svd.jl")
include("Bchain.jl")


# Slice matrix chain B_M .... B_1
export calc_product_chain
export calc_product_chain_stabilized


# QR / UDT
export UDT, inv!, udt_mult, Matrix!
export udt!, udt

export inv_one_plus!
export inv_one_plus
export udt_inv_one_plus
export udt_inv_sum
export inv_sum!
export inv_sum

export udt_inv_one_plus_loh
export inv_one_plus_loh!
export inv_one_plus_loh
export udt_inv_sum_loh
export inv_sum_loh!
export inv_sum_loh



# SVD / UDV
export gesdd!
export gesdd
export gesvd!
export gesvd
export genericsvd!
export genericsvd
export gesvj!
export gesvj

export udv!
export udv
export UDV_to_mat!
export UDV_to_mat
export mult_stable_udv

export inv_udv!
export inv_udv
export inv_one_plus_udv
export inv_one_plus_udv_alt
export inv_one_plus_udv_loh!
export inv_one_plus_udv_loh
export inv_sum_udvs!
export inv_sum_udvs





end # module