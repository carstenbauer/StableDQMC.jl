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