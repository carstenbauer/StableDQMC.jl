module StableDQMC

using LinearAlgebra, Pkg
using Requires

include("helpers.jl")
include("qr_udt.jl")
include("svd.jl")


# Slice matrix chain B_M .... B_1
export calc_Bchain, calc_Bchain_qr, calc_Bchain_svd
export calc_tdgf_qr, calc_tdgf_svd


# QR / UDT
export UDT, inv!, fact_mult, Matrix!
export udt!, udt

export udt_inv_one_plus
export inv_one_plus!
export inv_one_plus

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

export svd_inv_one_plus
export inv_one_plus!
export inv_one_plus

export svd_inv_one_plus_loh
export inv_one_plus_loh!
export inv_one_plus_loh

export svd_inv_sum
export inv_sum!
export inv_sum

export svd_inv_sum_loh
export inv_sum_loh!
export inv_sum_loh


# Optional JacobiSVD + GenericSVD
function addJacobiSVD()
    pkg"add https://github.com/RalphAS/JacobiSVD.jl"
    @eval import JacobiSVD
end
rmJacobiSVD() = pkg"rm JacobiSVD"

function addGenericSVD()
    pkg"add GenericSVD"
    @eval import GenericSVD
end
rmGenericSVD() = pkg"rm GenericSVD"

function __init__()
    @require JacobiSVD="2ca068c6-2156-5cf0-8317-c67edf277a2c" include("svd_jacobi.jl")
    @require GenericSVD="01680d73-4ee2-5a08-a1aa-533608c188bb" include("svd_generic.jl")
end


end # module
