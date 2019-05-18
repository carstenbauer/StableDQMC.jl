function udt(A::AbstractMatrix{C}) where C<:Number
  F = qr(A, Val(true))
  p = F.p
  @views p[p] = collect(1:length(p))
  D = abs.(real(diag(F.R)))
  T = (Diagonal(1 ./ D) * F.R)[:, p]
  return Matrix(F.Q), D, T
end

function udt!(A::AbstractMatrix{C}, D) where C<:Number
  n = length(D)
  F = qr!(A, Val(true))
  R = F.R # F.R is of regular matrix type

  @views F.p[F.p] = 1:n

  @inbounds for i in 1:n
    D[i] = abs(real(R[i,i]))
  end

  lmul!(Diagonal(1 ./ D), R)

  return Matrix(F.Q), R[:, F.p] # Q, (D is modified in-place), T
end




##############################################################
#
#                   QR / UDT
#
##############################################################
"""
Safely multiply two UDT decompositions
"""
function mult_stable(Ul,Dl,Tl, Ur,Dr,Tr)
  mat = Tl * Ur
  lmul!(Diagonal(Dl), mat)
  rmul!(mat, Diagonal(Dr))
  U, D, T = udt(mat)
  return Ul*U, D, T*Tr
end








"""

  inv_udt(U, D, T) -> result

Calculate (UDT)^(-1).
"""
function inv_udt(U,D,T)
  res = similar(U)
  inv_udt!(res, U, D, T)
  res
end

"""

  inv_udt!(res, U, D, T) -> nothing

Calculate (UDT)^(-1) and store result in `res`
"""
function inv_udt!(res, U,D,T)
  tmp = similar(U)
  ldiv!(tmp, lu(T), Diagonal(1 ./ D))
  mul!(res, tmp, U')
  nothing
end









"""

  inv_one_plus_udt(U, D, T) -> result

Stable calculation of [1 + UDT]^(-1):

  * Use one intermediate UDT decomposition.

Faster but less accurate than the loh approach.

Consider `inv_one_plus_udt!` as an efficient (not one-to-one) replacement.
"""
function inv_one_plus_udt(U,D,T)
  @warn "Calling somewhat inefficient and potentially inaccurate `inv_one_plus_udt`"

  m = U' / T
  m[diagind(m)] .+= D
  u,d,t = udt(m)
  u = U * u
  t = t * T
  ldiv!(m, lu!(t), Diagonal(1 ./ d))
  m * u'
end











"""

  inv_one_plus_udt!(mc, res, U, D, T) -> nothing

Stable calculation of [1 + UDT]^(-1):

  * Use one intermediate UDT decomposition.

Uses preallocated memory in `mc`. Writes the result into `res`.

Much faster (~50%) than `inv_one_plus_udt_loh!` but less accurate.
"""
function inv_one_plus_udt!(mc, res, U,D,T)
  @warn "Calling potentially inaccurate `inv_one_plus_udt!`"

  d = mc.s.d
  u = mc.s.tmp
  t = mc.s.tmp2

  m = U' / T
  m[diagind(m)] .+= D

  utmp,ttmp = udt!(m, d)
  mul!(u, U, utmp)
  mul!(t, ttmp, T)
    
  ldiv!(m, lu!(t), Diagonal(1 ./ d))
    
  mul!(res, m, u')
  nothing
end















"""

  inv_one_plus_two_udts!(mc, U, D, T, Ul, Dl, Tl, Ur, Dr, Tr) -> nothing

Stable calculation of [1 + UlDlTl(UrDrTr)^†]^(-1).

Uses preallocated memory in `mc`. Writes the result into `U`, `D`, and `T`.
"""
function inv_one_plus_two_udts!(mc, U,D,T, Ul,Dl,Tl, Ur,Dr,Tr)
  s = mc.s
  tmp = mc.s.tmp
  tmp2 = mc.s.tmp2
  tmp3 = mc.s.curr_U

  mul!(tmp, Tl, adjoint(Tr))
  rmul!(tmp, Diagonal(Dr))
  lmul!(Diagonal(Dl), tmp)
  U1, T1 = udt!(tmp, s.D)

  mul!(tmp3, Ul, U1)
  mul!(tmp2, T1, adjoint(Ur))
  mul!(tmp, adjoint(tmp3), inv(tmp2))

  tmp .+= Diagonal(s.D)

  u, t = udt!(tmp, D)

  mul!(tmp, t, tmp2)
  copyto!(U, inv(tmp))

  mul!(tmp, tmp3, u)
  copyto!(T, adjoint(tmp))

  D .= 1 ./ D

  nothing
end




"""

  inv_one_plus_two_udts!(mc, res, Ul, Dl, Tl, Ur, Dr, Tr) -> nothing

Stable calculation of [1 + UlDlTl(UrDrTr)^†]^(-1).

Uses preallocated memory in `mc`. Writes the result into `res`.
"""
function inv_one_plus_two_udts!(mc, res, Ul,Dl,Tl, Ur,Dr,Tr)
  inv_one_plus_two_udts!(mc, s.U, s.d, s.T, Ul, Dl, Tl, Ur, Dr, Tr)
  rmul!(s.U, Diagonal(s.d))
  mul!(res, s.U, s.T)
  nothing
end















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
















"""
  UDT_to_mat!(mat, U, D, T[; invert=false]) -> nothing

Combine given UDT to a matrix. Result will be stored in `mat`.

If `invert = true` (UDT)^(-1) is calculated.

Note that if parts of the input UDT may be overwritten it's
more efficient to inline this function and remove unnecessary copies!
"""
function UDT_to_mat!(mat, U, D, T; invert=false) 
    # TODO: Check where this is called and inline manually for better speed. (remove copy())
    @warn "UDT_to_mat! probably shouldn't be called here" # TODO remove once checked
    if !invert
        mat1 = copy(U)
        rmul!(mat1, Diagonal(D))
        mul!(mat,mat1,T)
    else # (DT)^-1 * U^dagger
        mat1 = copy(T)
        lmul!(Diagonal(D), mat1)
        mat .= mat1 \ U'
    end
    nothing
end

function UDT_to_mat(U, D, T; kw...)
  res = copy(U)
  UDT_to_mat!(res, U, D, T; kw...)
  res
end