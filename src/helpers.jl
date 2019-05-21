"""
Calculate the condition number from the given singular values
"""
function LinearAlgebra.cond(svs::AbstractVector)
    mi, ma = extrema(svs)
    ma / mi
end


"""

calc_product_chain(B,N) -> (R, svs)

Calculate B^N as B*B*B*B*...*B from right to left.

Returns the result `R` and the singular values (column of `svs`) of all intermediate powers (rows of `svs`).
"""
function calc_product_chain(B, N)
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
Does the same as `calc_product_chain` but stabilizes the
matrix products by intermediate matrix decompositions.
"""
function calc_product_chain_stabilized(B, N, decomposition_method)
    svs = zeros(eltype(B), N,size(B,1))
    U, D, X = decomposition_method(B)
    svs[1,:] = log.(D)
    svc = 2

    for k in 2:N
        # multiply B from the left to UDX
        # and update UDX
        U = B * U
        U *= Diagonal(D)
        U, D, Xnew = decomposition_method(U)
        X = Xnew * X
        
        # keep singular values
        svs[svc,:] = log.(D)
        svc += 1
    end
    
    return U, D, X, svs
end