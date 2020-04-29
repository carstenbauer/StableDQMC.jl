#
#
#   These methods are for testing purposes only.
#   They are not optimized at all.
#
#


"""
Calculates `B_M B_M-1 ... B_1` and stabilizes the
matrix products by intermediate matrix decompositions.
Assumes that the input `Bs` are ordered as `[B_1, B_2, ..., B_M]`.
"""
function calc_Bchain_svd(Bs; svdalg = svd)
    M = length(Bs)
    svs = zeros(Float64, M,size(Bs[1],1))
    F = svdalg(Bs[1])
    U, S, Vt = F.U, F.S, F.Vt
    svs[1,:] = log.(S)
    svc = 2

    for k in 2:M
        # multiply B from the left to USVt
        # and update USVt
        U = Bs[k] * U
        U *= Diagonal(S)
        M = svdalg(U)
        U, S, Vtnew = M.U, M.S, M.Vt
        Vt = Vtnew * Vt

        # keep singular values
        svs[svc,:] = log.(S)
        svc += 1
    end

    return SVD(U, S, Vt), svs
end


function calc_Bchain_svd_faster(Bs; svdalg = svd)
    M = length(Bs)
    flv = 4
    N = 100
    G = eltype(eltype(Bs))

    start = 1
    stop = M

    U = Matrix{G}(I, flv*N, flv*N)
    S = ones(Float64, flv*N)
    Vt = Matrix{G}(I, flv*N, flv*N)

    tmp = copy(U)
    tmp2 = copy(Vt)

    svs = zeros(length(start:stop),flv*N)
    svc = 1
    for k in start:stop
        if mod(k,1) == 0 || k == stop # always decompose in the end
            mul!(tmp, Bs[k], U) # U = Bs[k] * U
            rmul!(tmp, Diagonal(S))
            M = svdalg(tmp)
            U, S, Vtnew = M.U, M.S, M.Vt
            mul!(tmp2, Vtnew, Vt)
            Vt .= tmp2

            svs[svc,:] = log.(S)
            svc += 1
        else
            U = Bs[k] * U
        end
    end
    return SVD(U, S, Vt), svs
end

"""
Calculates `B_M B_M-1 ... B_1` and stabilizes the
matrix products by intermediate matrix decompositions.
Assumes that the input `Bs` are ordered as `[B_1, B_2, ..., B_M]`.
"""
function calc_Bchain_qr(Bs)
    M = length(Bs)
    svs = zeros(Float64, M,size(Bs[1],1))
    U, D, T = udt(Bs[1])
    svs[1,:] = log.(D)
    svc = 2

    for k in 2:M
        # multiply B from the left to UDT
        # and update UDT
        U = Bs[k] * U
        U *= Diagonal(D)
        U, D, Tnew = udt(U)
        T = Tnew * T

        # keep singular values
        svs[svc,:] = log.(D)
        svc += 1
    end

    return UDT(U, D, T), svs
end



# """
# Calculate fake "G(τ, 0)", i.e. [B^-N1 + B^N2]^-1
# """
# function calc_tdgf_qr(B, N1, N2; inv_sum_method = inv_sum)
#   if N1 != 0
#     udtl = calc_Bchain_qr(inv(B), N1)[1]
#   else
#     udtl = udt!(Matrix(I*one(eltype(B)), size(B)))
#   end
#
#   if N2 != 0
#     udtr = calc_Bchain_qr(B, N2)[1]
#   else
#     udtr = udt!(Matrix(I*one(eltype(B)), size(B)))
#   end
#
#   inv_sum_method(udtl, udtr)
# end
#
#
# """
# Calculate fake "G(τ, 0)", i.e. [B^-N1 + B^N2]^-1
# """
# function calc_tdgf_svd(B, N1, N2; inv_sum_method = inv_sum, svdalg_inv = svd!, svdalg_chain = svd)
#   # svdalg_nonmutating = Symbol(replace(string(nameof(svdalg)), "!" => ""))
#   if N1 != 0
#     svdl = calc_Bchain_svd(inv(B), N1; svdalg = svdalg_chain)[1]
#   else
#     svdl = svdalg(Matrix(I*one(eltype(B)), size(B)))
#   end
#
#   if N2 != 0
#     svdr = calc_Bchain_svd(B, N2; svdalg = svdalg_chain)[1]
#   else
#     svdr = svdalg(Matrix(I*one(eltype(B)), size(B)))
#   end
#
#   inv_sum_method(svdl, svdr; svdalg = svdalg_inv)
# end
