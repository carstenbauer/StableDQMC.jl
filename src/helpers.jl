#
#
#   These methods are for testing purposes only.
#   They are not optimized at all.
#
#

"""
Calculate the condition number from the given singular values
"""
function LinearAlgebra.cond(svs::AbstractVector)
    mi, ma = extrema(svs)
    ma / mi
end


"""

calc_Bchain(B,N) -> (R, svs)

Calculate B^N as B*B*B*B*...*B from right to left.

Returns the result `R` and the singular values (column of `svs`) of all intermediate powers (rows of `svs`).
"""
function calc_Bchain(B, N)
    bs = size(B,1)
    R = Matrix{eltype(B)}(I, bs, bs)
    svs = zeros(eltype(B), N,bs)
    svc = 1
    for k in 1:N
        R = B * R
        svs[svc,:] = log.(svdvals(R))
        svc += 1
    end
    return (R, svs)
end


"""
Same as `calc_Bchain` but stabilizes the
matrix products by intermediate matrix decompositions.
"""
function calc_Bchain_svd(B, N; svdalg = svd)
    svs = zeros(eltype(B), N,size(B,1))
    F = svdalg(B)
    U, S, Vt = F.U, F.S, F.Vt
    svs[1,:] = log.(S)
    svc = 2

    for k in 2:N
        # multiply B from the left to USVt
        # and update USVt
        U = B * U
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



"""
Same as `calc_Bchain` but stabilizes the
matrix products by intermediate matrix decompositions.
"""
function calc_Bchain_qr(B, N)
    svs = zeros(eltype(B), N,size(B,1))
    U, D, T = udt(B)
    svs[1,:] = log.(D)
    svc = 2

    for k in 2:N
        # multiply B from the left to UDT
        # and update UDT
        U = B * U
        U *= Diagonal(D)
        U, D, Tnew = udt(U)
        T = Tnew * T
        
        # keep singular values
        svs[svc,:] = log.(D)
        svc += 1
    end
    
    return UDT(U, D, T), svs
end



"""
Calculate fake "G(τ, 0)", i.e. [B^-N1 + B^N2]^-1
"""
function calc_tdgf_qr(B, N1, N2; inv_sum_method = inv_sum)
  if N1 != 0
    udtl = calc_Bchain_qr(inv(B), N1) # TODO
  else
    udtl = udt!(Matrix(I*one(eltype(B)), size(B)))
  end

  if N2 != 0
    udtr = calc_Bchain_qr(B, N2)[1]
  else
    udtr = udt!(Matrix(I*one(eltype(B)), size(B)))
  end

  inv_sum_method(udtl, udtr)
end


"""
Calculate fake "G(τ, 0)", i.e. [B^-N1 + B^N2]^-1
"""
function calc_tdgf_svd(B, N1, N2; inv_sum_method = inv_sum)
  if N1 != 0
    svdl = calc_Bchain_svd(inv(B), N1) # TODO
  else
    svdl = svd!(Matrix(I*one(eltype(B)), size(B)))
  end

  if N2 != 0
    svdr = calc_Bchain_svd(B, N2)[1]
  else
    svdr = svd!(Matrix(I*one(eltype(B)), size(B)))
  end

  inv_sum_method(svdl, svdr)
end