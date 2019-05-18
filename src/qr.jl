##############################################################
#
#                   QR / UDT
#
##############################################################
"""

    inv_one_plus(F::UDT) -> AbstractMatrix

Stable calculation of `[1 + UDT]^(-1)`:

  * Use one intermediate UDT decomposition.

Faster but potentially less accurate than `inv_one_plus_loh`.

Optional preallocations possible via keyword arguments:

  * `d = similar(F.D)`
  * `u = similar(F.U)`
  * `t = similar(F.T)`
"""
function inv_one_plus(F::UDT; kwargs...)
  res = similar(F.U)
  inv_one_plus!(res, F; kwargs...)
  return res
end




"""

  inv_one_plus!(res, F::UDT) -> res

Same as `inv_one_plus` but stores the result in preallocated `res`.
"""
function inv_one_plus!(res, F::UDT; u = similar(F.U), d = similar(F.D), t = similar(F.T))
  # @warn "Calling potentially inaccurate `inv_one_plus_udt!`"
  # d = mc.s.d
  # u = mc.s.tmp
  # t = mc.s.tmp2
  U, D, T = F

  m = U' / T
  m[diagind(m)] .+= D

  utmp, d, ttmp = udt!(m)
  mul!(u, U, utmp)
  mul!(t, ttmp, T)
    
  ldiv!(m, lu!(t), Diagonal(1 ./ d))
    
  mul!(res, m, u')
  res
end










"""

    udt_inv_one_plus(A::UDT, Bdagger::UDT) -> UDT

Stable calculation of [1 + UlDlTl(UrDrTr)^†]^(-1). Returns and
`UDT` factorization object.

Optional preallocations via keyword arguments:

  * tmp = similar(A.U)
  * tmp2 = similar(A.U)
  * tmp3 = similar(A.U)
"""
function udt_inv_one_plus(A::UDT, Bdagger::UDT; tmp = similar(A.U), tmp2 = similar(A.U), tmp3 = similar(A.U))
  # s = mc.s
  # tmp = mc.s.tmp
  # tmp2 = mc.s.tmp2
  # tmp3 = mc.s.curr_U
  Ul,Dl,Tl = A
  Ur,Dr,Tr = Bdagger

  mul!(tmp, Tl, adjoint(Tr))
  rmul!(tmp, Diagonal(Dr))
  lmul!(Diagonal(Dl), tmp)
  U1, D1, T1 = udt!(tmp)

  mul!(tmp3, Ul, U1)
  mul!(tmp2, T1, adjoint(Ur))
  mul!(tmp, adjoint(tmp3), inv(tmp2))

  tmp .+= Diagonal(D1)

  u, d, t = udt!(tmp)
  mul!(tmp, t, tmp2)
  mul!(tmp2, tmp3, u)

  UDT(inv(tmp), 1 ./ d, tmp2')
end





"""

  inv_one_plus!(res, A::UDT, Bdagger::UDT) -> res

Stable calculation of [1 + UlDlTl(UrDrTr)^†]^(-1).
Writes the result into `res`.

See `udt_inv_one_plus` for preallocation options.
"""
function inv_one_plus!(res, A::UDT, Bdagger::UDT; kwargs...)
  F = udt_inv_one_plus(A, Bdagger; kwargs...)
  rmul!(F.U, Diagonal(F.D))
  mul!(res, F.U, F.T)
  res
end


"""

  inv_one_plus(A::UDT, Bdagger::UDT) -> AbstractMatrix

Stable calculation of [1 + UlDlTl(UrDrTr)^†]^(-1).

See `udt_inv_one_plus` for preallocation options.
"""
function inv_one_plus(A::UDT, Bdagger::UDT; kwargs...)
  res = similar(A.U)
  inv_one_plus!(res, A, Bdagger; kwargs...)
  return res
end






################################################# Stopped here





"""

  inv_one_plus_udt_loh!(mc, res, Ua, Da, Ta, Ub, Db, Tb) -> nothing

Stable calculation of [1 + UDT]^(-1):

  * Separate scales larger and smaller than unity
  * Use two intermediate UDT decompositions.

Uses preallocated memory in `mc`. Writes the result into `res`.
"""
function inv_one_plus_udt_loh!(mc, res, U,D,T)
  d = mc.s.d
  r = mc.s.tmp
  l = mc.s.tmp2

  Dp = max.(D,1.)
  Dm = min.(D,1.)
  
  Dp .\= 1
  Dpinv = Dp # renaming

  ldiv!(l, lu(T), Diagonal(Dpinv)) # Don't use lu!(T) because T is input

  mul!(r, U, Diagonal(Dm))
  r .+= l

  u, t = udt!(r, d)

  ldiv!(r, lu!(t), Diagonal(1 ./ d))
  mul!(l, r, u')

  lmul!(Diagonal(Dpinv), l)
  u, t = udt!(l, d)

  ldiv!(l, lu(T), u)

  rmul!(l, Diagonal(d))
  mul!(res, l, t)
  nothing
end









"""

  inv_sum_udts(Ua, Da, Ta, Ub, Db, Tb) -> U, D, T

Stable calculation of [UaDaTa + UbDbTb]^(-1):

  * Use one intermediate UDT decompositions.

Faster but less accurate than `inv_sum_udts_loh`.

Consider `inv_sum_udts!` as an efficient (not one-to-one) replacement.
"""
function inv_sum_udts(Ua,Da,Ta,Ub,Db,Tb)
  @warn "Calling somewhat inefficient and potentially inaccurate `inv_sum_udts`"

  m1 = Ta / Tb
  lmul!(Diagonal(Da), m1)

  m2 = Ua' * Ub
  rmul!(m2, Diagonal(Db))

  u,d,t = udt(m1 + m2)

  mul!(m1, Ua, u)
  mul!(m2, t, Tb)

  return inv(m2), 1 ./ d, m1'
end









"""

  inv_sum_udts!(mc, res, Ua, Da, Ta, Ub, Db, Tb) -> nothing

Stable calculation of [UaDaTa + UbDbTb]^(-1):

  * Use one intermediate UDT decompositions.

Uses preallocated memory in `mc`. Writes the result into `res`.

Much faster (~40%) than `inv_sum_udts_loh!` but less accurate.
"""
function inv_sum_udts!(mc, res, Ua,Da,Ta,Ub,Db,Tb)
  @warn "Calling potentially inaccurate `inv_sum_udts!`"

  d = mc.s.d
  m1 = mc.s.tmp
  m2 = mc.s.tmp2
  tmp = mc.s.U

  m1 = Ta / Tb
  lmul!(Diagonal(Da), m1)

  mul!(m2, Ua', Ub)
  rmul!(m2, Diagonal(Db))

  u,t = udt!(m1 + m2, d)

  mul!(m1, Ua, u)
  mul!(m2, t, Tb)

  ldiv!(tmp, lu!(m2), Diagonal(1 ./ d))
  mul!(res, tmp, m1')

  nothing
end









"""

  inv_sum_udts_loh(Ua, Da, Ta, Ub, Db, Tb) -> U, D, T

Stable calculation of [UaDaTa + UbDbTb]^(-1):

  * Separate scales larger and smaller than unity
  * Use two intermediate UDT decompositions.

Consider `inv_sum_udts_loh!` as an efficient (not one-to-one) replacement.
"""
function inv_sum_udts_loh(Ua, Da, Ta, Ub, Db, Tb)
    @warn "Calling somewhat inefficient `inv_sum_udts_loh`"

    d=length(Da)
    
    # separating scales larger and smaller than unity
    Dap = max.(Da,1.)
    Dam = min.(Da,1.)
    Dbp = max.(Db,1.)
    Dbm = min.(Db,1.)

    # mat1 = Dam * Vda * Vdb' / Dbp
    mat1 = Ta * inv(Tb)
    @inbounds for j in 1:d, k in 1:d
        mat1[j,k]=mat1[j,k] * Dam[j]/Dbp[k]
    end

    # mat2 = 1/(Dap) * Ua' * Ub * Dbm
    mat2 = Ua' * Ub
    @inbounds for j in 1:d, k in 1:d
        mat2[j,k]=mat2[j,k] * Dbm[k]/Dap[j]
    end
    
    #mat1 = mat1 + mat2
    mat1 = mat1 + mat2
    
    # decompose mat1: U, D, T
    U, D, T = udt(mat1)

    # invert inner part: mat1 = (U D T)^(-1) = mat1^(-1)
    # was UDT_to_mat!(mat1, U, D, T, invert=true)
    lmul!(Diagonal(D), T)
    ldiv!(mat1, lu!(T), U') # mat1 = T \ (U')

    # mat1 = 1/Dbp * mat1 /Dap
    @inbounds for j in 1:d, k in 1:d
        mat1[j,k]=mat1[j,k] / Dbp[j] / Dap[k]
    end

    #mat1 = U D T
    U, D, T = udt(mat1)

    # U = Tb^(-1) * U , T = T * Ua'
    mul!(mat1, inv(Tb), U)
    mul!(mat2, T, Ua')
    U = mat1
    T = mat2

    return U, D, T
end









"""

  inv_sum_udts_loh!(mc, res, Ua, Da, Ta, Ub, Db, Tb) -> nothing

Stable calculation of [UaDaTa + UbDbTb]^(-1):

  * Separate scales larger and smaller than unity
  * Use two intermediate UDT decompositions.

Uses preallocated memory in `mc`. Writes the result into `res`.
"""
function inv_sum_udts_loh!(mc, res, Ua, Da, Ta, Ub, Db, Tb)
    # optimization rounds: 1
    mat1 = mc.s.tmp
    mat2 = mc.s.U
    D = mc.s.d
    to = mc.a.to
    
    d=length(Da)
        
    # separating scales larger and smaller than unity
    Dap = max.(Da,1.)
    Dam = min.(Da,1.)
    Dbp = max.(Db,1.)
    Dbm = min.(Db,1.)

    # mat1 = Dam * Vda * Vdb' / Dbp
    mat1 = Ta / Tb
    @inbounds for j in 1:d, k in 1:d
        mat1[j,k]=mat1[j,k] * Dam[j]/Dbp[k]
    end

    # mat2 = 1/(Dap) * Ua' * Ub * Dbm
    mul!(mat2, Ua', Ub)    
    @inbounds for j in 1:d, k in 1:d
        mat2[j,k]=mat2[j,k] * Dbm[k]/Dap[j]
    end
    
    # mat1 = mat1 + mat2
    mat1 .+= mat2
    
    # decompose mat1: U, D, T
    U, T = udt!(mat1, D)

    # invert and combine inner part: mat1 = (U D T)^(-1)
    lmul!(Diagonal(D), T)
    ldiv!(mat1, lu!(T), U') # mat1 = T \ (U')

    # mat1 = 1/Dbp * mat1 /Dap
    @inbounds for j in 1:d, k in 1:d
        mat1[j,k]=mat1[j,k] / Dbp[j] / Dap[k]
    end

    #mat1 = U D T
    U, T = udt!(mat1, D)

    # U = Tb^(-1) * U , T = T * Ua'
    ldiv!(mat1, lu!(Tb), U) # mat1 = Tb \ U
    mul!(mat2, T, Ua')
    U = mat1
    T = mat2

    # combine UDT into res
    rmul!(U, Diagonal(D))
    mul!(res, U, T)
    nothing
end