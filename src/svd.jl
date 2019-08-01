"""
Calculates the SVD deomposition of a matrix by employing
a divide and conquer approach.
Returns a `SVD` factorization object.
"""
gesdd(x) = gesdd!(copy(x))
"""
Same as `gesdd` but saves space by overwriting the input matrix.
"""
function gesdd!(x)
  u,s,vt = LinearAlgebra.LAPACK.gesdd!('A', x)
  SVD(u,s,vt)
end


"""
Calculates the SVD deomposition of a matrix.
Returns a `SVD` factorization object.
"""
gesvd(x) = gesvd!(copy(x))
"""
Same as `gesvd` but saves space by overwriting the input matrix.
"""
function gesvd!(x)
  u,s,vt = LinearAlgebra.LAPACK.gesvd!('A', 'A', x)
  SVD(u,s,vt)
end




# extra operations for SVD factorizations
"""
    fact_mult(A::SVD, B::SVD) -> SVD

Stabilized multiplication of two SVD decompositions.
Returns a `SVD` factorization object.
"""
function fact_mult(A::SVD, B::SVD)
    mat = A.Vt * B.U
    lmul!(Diagonal(A.S), mat)
    rmul!(mat, Diagonal(B.S))
    F = svd!(mat)
    SVD(A.U * F.U, F.S, F.Vt * B.Vt)
end


# """
#     *(A::SVD, B::SVD)

# Stabilized multiplication of two SVD decompositions.

# (Type piracy if `LinearAlgebra` defines `*` for `SVD`s in the future!)
# """
# function Base.:*(A::SVD, B::SVD)
#     mat = A.Vt * B.U
#     lmul!(Diagonal(A.S), mat)
#     rmul!(mat, Diagonal(B.S))
#     F = svd!(mat)
#     (A.U * F.U) * Diagonal(F.S) * (F.Vt * B.Vt)
# end










##############################################################
#
#                   SVD / USVt
#
##############################################################

"""
    svd_inv_one_plus(F::SVD) -> SVD

Stabilized calculation of [1 + USVt]^(-1):

  * Use one intermediate SVD decomposition.

Options for preallocation via keyword arguments:

  * `t = similar(F.U)`
  * `a = similar(F.U)`
  * `c = similar(F.U)`
"""
function svd_inv_one_plus(F::SVD;
                          svdalg = svd!,
                          t = similar(F.U),
                          a = similar(F.U),
                          c = similar(F.U),
                          internaluse=false)
  U, S, V = F
  mul!(t, U', V)
  t[diagind(t)] .+= S
  u, s, v = svdalg(t)
  mul!(a, V, v)
  s .= 1 ./ s
  mul!(c, u', U')
  if internaluse
    SVD(a, s, c)
  else
    SVD(copy(a), s, copy(c))
  end
end


"""
    inv_one_plus!(res, F::SVD) -> res

Stabilized calculation of [1 + USVt]^(-1):

  * Use one intermediate SVD decomposition.

Writes the result into `res`.

See `svd_inv_one_plus` for preallocation options.
"""
function inv_one_plus!(res, F::SVD; kwargs...)
  M = svd_inv_one_plus(F; internaluse=true, kwargs...)
  rmul!(M.U, Diagonal(M.S))
  mul!(res, M.U, M.Vt)
  res
end


"""
    inv_one_plus(F::SVD) -> AbstractMatrix

Stabilized calculation of [1 + USVt]^(-1):

  * Use one intermediate SVD decomposition.

See `svd_inv_one_plus` for preallocation options.
"""
function inv_one_plus(F::SVD; kwargs...)
  res = similar(F.U)
  inv_one_plus!(res, F; kwargs...)
  res
end



"""
    svd_inv_one_plus_loh(F::SVD) -> SVD

Stabilized calculation of [1 + USVt]^(-1):

  * Separate scales larger and smaller than unity
  * Use two intermediate SVD decompositions.

Options for preallocation via keyword arguments:

  * `l = similar(F.U)`
  * `r = similar(F.U)`
  * `Dp = similar(F.D)`
  * `Dm = similar(F.D)`
"""
function svd_inv_one_plus_loh(F::SVD;
                              svdalg = svd!,
                              Sp = similar(F.S),
                              Sm = similar(F.S),
                              l = similar(F.V),
                              r = similar(F.U))
  U, S, V = F

  copyto!(l, V)
  copyto!(r, U)

  Sp .= max.(S, 1)
  Sm .= min.(S, 1)

  Sp .\= 1 # Sp now Spinv!!!

  rmul!(l, Diagonal(Sp))
  rmul!(r, Diagonal(Sm))

  M = svdalg(l + r)
  m = inv(M)
  lmul!(Diagonal(Sp), m)
  u, s, v = svdalg(m)
  mul!(m, V, u)

  SVD(m, s, v')
end


"""
  inv_one_plus_loh!(res, F::SVD) -> res

Stabilized calculation of [1 + USVt]^(-1):

  * Separate scales larger and smaller than unity
  * Use two intermediate SVD decompositions.

Writes the result into `res`.

See `svd_inv_one_plus_loh` for preallocation options.
"""
function inv_one_plus_loh!(res, F::SVD; kwargs...)
  M = svd_inv_one_plus_loh(F; kwargs...)
  rmul!(M.U, Diagonal(M.S))
  mul!(res, M.U, M.Vt)
  res
end


"""
  inv_one_plus_loh(F::SVD) -> AbstractMatrix

Stabilized calculation of [1 + USVt]^(-1):

  * Separate scales larger and smaller than unity
  * Use two intermediate SVD decompositions.

See `svd_inv_one_plus_loh` for preallocation options.
"""
function inv_one_plus_loh(F::SVD; kwargs...)
  res = similar(F.U)
  inv_one_plus_loh!(res, F; kwargs...)
end








"""
    svd_inv_sum(A::SVD, B::SVD) -> SVD

Stabilized calculation of [UaSaVta + UbSbVtb]^(-1):

  * Use one intermediate SVD decompositions.

Options for preallocations via keyword arguments:

  * `m1 = similar(A.U)`
  * `m2 = similar(A.U)`

"""
function svd_inv_sum(A::SVD, B::SVD; m1 = similar(A.U), m2 = similar(A.U), svdalg = svd!, internaluse = false)
  Ua, Sa, Va = A
  Ub, Sb, Vb = B

  mul!(m1, Va', Vb)
  lmul!(Diagonal(Sa), m1)

  mul!(m2, Ua', Ub)
  rmul!(m2, Diagonal(Sb))

  u,s,v = svdalg(m1 + m2)

  mul!(m1, Ua, u)
  mul!(m2, Vb, v)

  if internaluse
    SVD(m2, 1 ./ s, m1')
  else
    SVD(copy(m2), 1 ./ s, copy(m1'))
  end
end


"""
    inv_sum!(res, A::SVD, B::SVD) -> res

Stabilized calculation of [UaSaVta + UbSbVtb]^(-1):

  * Use one intermediate SVD decompositions.

Writes the result into `res`.

See `svd_inv_sum` for preallocation options.
"""
function inv_sum!(res, A::SVD, B::SVD; kwargs...)
  M = svd_inv_sum(A, B; internaluse = true, kwargs...)
  rmul!(M.U, Diagonal(M.S))
  mul!(res, M.U, M.Vt)
  res
end

"""
    inv_sum(A::SVD, B::SVD) -> res

Stabilized calculation of [UaSaVta + UbSbVtb]^(-1):

  * Use one intermediate SVD decompositions.

See `svd_inv_sum` for preallocation options.
"""
function inv_sum(A::SVD, B::SVD; kwargs...)
  res = similar(A.U)
  inv_sum!(res, A, B; kwargs...)
  res
end





"""
    svd_inv_sum_loh(A::SVD, B::SVD) -> SVD

Stabilized calculation of [UaSaVta + UbSbVtb]^(-1):

  * Separate scales larger and smaller than unity
  * Use two intermediate SVD decompositions.

Options for preallocations via keyword arguments:

  * None (yet)

"""
function svd_inv_sum_loh(A::SVD, B::SVD; svdalg = svd!)
    Ua, Sa, Va = A
    Ub, Sb, Vb = B

    d = length(Sa)

    Sap = max.(Sa,1.)
    Sam = min.(Sa,1.)
    Sbp = max.(Sb,1.)
    Sbm = min.(Sb,1.)

    mat1 = Va' * Vb
    for j in 1:d, k in 1:d
        mat1[j,k]=mat1[j,k] * Sam[j]/Sbp[k]
    end

    mat2 = adjoint(Ua) * Ub
    for j in 1:d, k in 1:d
        mat2[j,k]=mat2[j,k] * Sbm[k]/Sap[j]
    end

    mat1 = mat1 + mat2

    M = svdalg(mat1)
    mat1 = inv(M)

    for j in 1:d, k in 1:d
        mat1[j,k]=mat1[j,k] / Sbp[j] / Sap[k]
    end

    u, s, v = svdalg(mat1)
    mul!(mat1, Vb, u)
    mul!(mat2, v', adjoint(Ua))

    SVD(mat1, s, mat2)
end


"""
    inv_sum_loh!(res, A::SVD, B::SVD) -> res

Stabilized calculation of [UaSaVta + UbSbVtb]^(-1):

  * Separate scales larger and smaller than unity
  * Use two intermediate SVD decompositions.

Writes the result into `res`.

See `svd_inv_sum_loh` for preallocation options.
"""
function inv_sum_loh!(res, A::SVD, B::SVD; kwargs...)
  M = svd_inv_sum_loh(A, B; kwargs...)
  rmul!(M.U, Diagonal(M.S))
  mul!(res, M.U, M.Vt)
  res

end

"""
    inv_sum_loh(A::SVD, B::SVD) -> res

Stabilized calculation of [UaSaVta + UbSbVtb]^(-1):

  * Separate scales larger and smaller than unity
  * Use two intermediate SVD decompositions.

See `svd_inv_sum_loh` for preallocation options.
"""
function inv_sum_loh(A::SVD, B::SVD; kwargs...)
  res = similar(A.U)
  inv_sum_loh!(res, A, B; kwargs...)
  res
end
