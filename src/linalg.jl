#### SVD, i.e. UDV decomposition
decompose_udv!(A::AbstractMatrix{<:Number}) = LinearAlgebra.LAPACK.gesvd!('A','A',A)
decompose_udv(A::AbstractMatrix{T}) where T<:Number = decompose_udv!(copy(A))


#### QR, i.e. UDT decomposition
function decompose_udt(A::AbstractMatrix{C}) where C<:Number
  F = qr(A, Val(true))
  p = F.p
  @views p[p] = collect(1:length(p))
  D = abs.(real(diag(F.R)))
  T = (Diagonal(1 ./ D) * F.R)[:, p]
  return Matrix(F.Q), D, T
end

function decompose_udt!(A::AbstractMatrix{C}, D) where C<:Number
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



# See https://discourse.julialang.org/t/asymmetric-speed-of-in-place-sparse-dense-matrix-product/10256/3
import LinearAlgebra.mul!
function mul!(C::StridedMatrix, X::StridedMatrix, A::SparseMatrixCSC)
    mX, nX = size(X)
    nX == A.m || throw(DimensionMismatch())
    fill!(C, zero(eltype(C)))
    rowval = A.rowval
    nzval = A.nzval
    @inbounds for  col = 1:A.n, k=A.colptr[col]:(A.colptr[col+1]-1)
        ki=rowval[k]
        kv=nzval[k]
        for multivec_row=1:mX
            C[multivec_row, col] += X[multivec_row, ki] * kv
        end
    end
    C
end

# LinearAlgebra.Diagonal(I::UniformScaling{T}) where T <: Number = I.λ









##############################################################
#
#                   QR / UDT
#
##############################################################
"""
Safely multiply two UDT decompositions
"""
function multiply_safely(Ul,Dl,Tl, Ur,Dr,Tr)
  mat = Tl * Ur
  lmul!(Diagonal(Dl), mat)
  rmul!(mat, Diagonal(Dr))
  U, D, T = decompose_udt(mat)
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

Faster but less accurate than the scalettar approach.

Consider `inv_one_plus_udt!` as an efficient (not one-to-one) replacement.
"""
function inv_one_plus_udt(U,D,T)
  @warn "Calling somewhat inefficient and potentially inaccurate `inv_one_plus_udt`"

  m = U' / T
  m[diagind(m)] .+= D
  u,d,t = decompose_udt(m)
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

Much faster (~50%) than `inv_one_plus_udt_scalettar!` but less accurate.
"""
function inv_one_plus_udt!(mc, res, U,D,T)
  @warn "Calling potentially inaccurate `inv_one_plus_udt!`"

  d = mc.s.d
  u = mc.s.tmp
  t = mc.s.tmp2

  m = U' / T
  m[diagind(m)] .+= D

  utmp,ttmp = decompose_udt!(m, d)
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
  U1, T1 = decompose_udt!(tmp, s.D)

  mul!(tmp3, Ul, U1)
  mul!(tmp2, T1, adjoint(Ur))
  mul!(tmp, adjoint(tmp3), inv(tmp2))

  tmp .+= Diagonal(s.D)

  u, t = decompose_udt!(tmp, D)

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

  inv_one_plus_udt_scalettar!(mc, res, Ua, Da, Ta, Ub, Db, Tb) -> nothing

Stable calculation of [1 + UDT]^(-1):

  * Separate scales larger and smaller than unity
  * Use two intermediate UDT decompositions.

Uses preallocated memory in `mc`. Writes the result into `res`.
"""
function inv_one_plus_udt_scalettar!(mc, res, U,D,T)
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

  u, t = decompose_udt!(r, d)

  ldiv!(r, lu!(t), Diagonal(1 ./ d))
  mul!(l, r, u')

  lmul!(Diagonal(Dpinv), l)
  u, t = decompose_udt!(l, d)

  ldiv!(l, lu(T), u)

  rmul!(l, Diagonal(d))
  mul!(res, l, t)
  nothing
end









"""

  inv_sum_udts(Ua, Da, Ta, Ub, Db, Tb) -> U, D, T

Stable calculation of [UaDaTa + UbDbTb]^(-1):

  * Use one intermediate UDT decompositions.

Faster but less accurate than `inv_sum_udts_scalettar`.

Consider `inv_sum_udts!` as an efficient (not one-to-one) replacement.
"""
function inv_sum_udts(Ua,Da,Ta,Ub,Db,Tb)
  @warn "Calling somewhat inefficient and potentially inaccurate `inv_sum_udts`"

  m1 = Ta / Tb
  lmul!(Diagonal(Da), m1)

  m2 = Ua' * Ub
  rmul!(m2, Diagonal(Db))

  u,d,t = decompose_udt(m1 + m2)

  mul!(m1, Ua, u)
  mul!(m2, t, Tb)

  return inv(m2), 1 ./ d, m1'
end









"""

  inv_sum_udts!(mc, res, Ua, Da, Ta, Ub, Db, Tb) -> nothing

Stable calculation of [UaDaTa + UbDbTb]^(-1):

  * Use one intermediate UDT decompositions.

Uses preallocated memory in `mc`. Writes the result into `res`.

Much faster (~40%) than `inv_sum_udts_scalettar!` but less accurate.
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

  u,t = decompose_udt!(m1 + m2, d)

  mul!(m1, Ua, u)
  mul!(m2, t, Tb)

  ldiv!(tmp, lu!(m2), Diagonal(1 ./ d))
  mul!(res, tmp, m1')

  nothing
end









"""

  inv_sum_udts_scalettar(Ua, Da, Ta, Ub, Db, Tb) -> U, D, T

Stable calculation of [UaDaTa + UbDbTb]^(-1):

  * Separate scales larger and smaller than unity
  * Use two intermediate UDT decompositions.

Consider `inv_sum_udts_scalettar!` as an efficient (not one-to-one) replacement.
"""
function inv_sum_udts_scalettar(Ua, Da, Ta, Ub, Db, Tb)
    @warn "Calling somewhat inefficient `inv_sum_udts_scalettar`"

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
    U, D, T = decompose_udt(mat1)

    # invert inner part: mat1 = (U D T)^(-1) = mat1^(-1)
    # was UDT_to_mat!(mat1, U, D, T, invert=true)
    lmul!(Diagonal(D), T)
    ldiv!(mat1, lu!(T), U') # mat1 = T \ (U')

    # mat1 = 1/Dbp * mat1 /Dap
    @inbounds for j in 1:d, k in 1:d
        mat1[j,k]=mat1[j,k] / Dbp[j] / Dap[k]
    end

    #mat1 = U D T
    U, D, T = decompose_udt(mat1)

    # U = Tb^(-1) * U , T = T * Ua'
    mul!(mat1, inv(Tb), U)
    mul!(mat2, T, Ua')
    U = mat1
    T = mat2

    return U, D, T
end









"""

  inv_sum_udts_scalettar!(mc, res, Ua, Da, Ta, Ub, Db, Tb) -> nothing

Stable calculation of [UaDaTa + UbDbTb]^(-1):

  * Separate scales larger and smaller than unity
  * Use two intermediate UDT decompositions.

Uses preallocated memory in `mc`. Writes the result into `res`.
"""
function inv_sum_udts_scalettar!(mc, res, Ua, Da, Ta, Ub, Db, Tb)
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
    U, T = decompose_udt!(mat1, D)

    # invert and combine inner part: mat1 = (U D T)^(-1)
    lmul!(Diagonal(D), T)
    ldiv!(mat1, lu!(T), U') # mat1 = T \ (U')

    # mat1 = 1/Dbp * mat1 /Dap
    @inbounds for j in 1:d, k in 1:d
        mat1[j,k]=mat1[j,k] / Dbp[j] / Dap[k]
    end

    #mat1 = U D T
    U, T = decompose_udt!(mat1, D)

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
































##############################################################
#
#                   SVD / UDV(t)
#
##############################################################
# multiplies two UDVds -> UDVd
function multiply_safely_udv(Ul,Dl,Vdl,Ur,Dr,Vdr)
  tmp = adjoint(Vdl) * Ur
  rmul!(tmp, Diagonal(Dr))
  lmul!(Diagonal(Dl), tmp)
  U, D, Vd = decompose_udv!(tmp)
  U = Ul * U
  Vd = Vd * Vdr
  return U, D, Vd
end


# Calculates (UDVd)^-1, where U, D, Vd come from SVD decomp.
function inv_udv(U,D,Vd)
  m = copy(Vd')
  rmul!(m, Diagonal(1 ./ D))
  res = similar(m)
  mul!(res, m, U')
  res
end


# Calculates (UDVd)^-1, where U, D, Vd come from SVD decomp.
function inv_udv!(res, U,D,Vd)

  # copy here isn't necessary but Vd would be overwritten if we drop it
  # also, m could be preallocated in stack if necessary
  m = copy(Vd')
  rmul!(m, Diagonal(1 ./ D))
  mul!(res, m, U')
  nothing
end



# Calculates (1 + UDVd)^-1, where U, D, Vd come from SVD decomp.
# !! Breaks down for large spread in D (i.e. low temperatures).
function inv_one_plus_udv(U,D,Vd)
  inner = copy(Vd')
  inner .+= U * Diagonal(D)
  I = decompose_udv!(inner)
  u = copy(adjoint(I[3] * Vd))
  d = 1 ./ I[2]
  vd = adjoint(I[1])

  rmul!(u,Diagonal(d))
  u * vd
end

# same as inv_one_plus_udt but separating both U AND Vd from D
# !! Breaks down for large spread in D (i.e. low temperatures). Slightly better than normal version.
function inv_one_plus_udv_alt(U,D,Vd)
  inner = copy((Vd*U)')
  inner[diagind(inner)] .+= D
  u, d, vd = decompose_udv!(inner)

  t1 = adjoint(vd*Vd)
  t2 = adjoint(U*u)
  rmul!(t1, Diagonal(1 ./ d))
  t1*t2
end

# Calculates (1 + UDVd)^-1, where U, D, Vd come from SVD decomp.
# More controlled handling of scales, however also slower.
function inv_one_plus_udv_scalettar(U,D,Vd)
  Dp = max.(D,1.)
  Dm = min.(D,1.)
  Dpinv = 1 ./ Dp

  l = copy(Vd')
  rmul!(l, Diagonal(Dpinv))

  r = copy(U)
  rmul!(r, Diagonal(Dm))

  u, d, vd = decompose_udv!(l+r)

  m = inv_udv(u,d,vd)
  lmul!(Diagonal(Dpinv), m)
  u, d, vd = decompose_udv!(m)

  mul!(m, Vd', u)
  # return m, d, vd
  rmul!(m, Diagonal(d))
  m*vd
end


# TODO: Optimize!
# I only made the function overwrite res. Otherwise it's unchanged compared to the one above.
function inv_one_plus_udv_scalettar!(mc, res, U,D,Vd)
  # all similars here could go into mc.s
  Dp = similar(D)
  Dm = similar(D)
  l = similar(Vd)
  l .= Vd'
  r = similar(U)
  r .= U
  # u, d, vd below could be preallocated, but does this matter for speed?
  # m could be preallocated, I guess

  Dp .= max.(D, 1)
  Dm .= min.(D, 1)

  Dp .\= 1 # Dp now Dpinv!!!

  rmul!(l, Diagonal(Dp))
  rmul!(r, Diagonal(Dm))

  u, d, vd = decompose_udv!(l + r)

  m = inv_udv(u,d,vd) # TODO: optimize
  lmul!(Diagonal(Dp), m)
  u, d, vd = decompose_udv!(m)

  mul!(m, Vd', u)
  # return m, d, vd
  rmul!(m, Diagonal(d))
  # res .= m*vd
  mul!(res, m, vd)
  nothing
end







function UDV_to_mat!(mat, U, D, Vd; invert=false) 
    if !invert
        mat1 = copy(U)
        rmul!(mat1, Diagonal(D))
        mul!(mat, mat1, Vd)
    else #V D^(-1) Ud = (D^-1 *Vd)^(dagger) *Ud
        mat1 = copy(Vd)
        lmul!(Diagonal(1 ./ D), mat1)
        mul!(mat, mat1', U')
    end
    nothing
end


function UDV_to_mat(U, D, Vd; kw...)
  res = copy(U)
  UDV_to_mat!(res, U, D, Vd; kw...)
  res
end





# Calculates (UaDaVda + UbDbVdb)^-1
function inv_sum_udvs(Ua, Da, Vda, Ub, Db, Vdb)
    
    d=length(Da)
    
    # separating scales larger and smaller than unity
    Dap = max.(Da,1.)
    Dam = min.(Da,1.)
    Dbp = max.(Db,1.)
    Dbm = min.(Db,1.)

    # mat1 = Dam * Vda * Vdb' / Dbp

    mat1 = Vda * adjoint(Vdb)
    for j in 1:d, k in 1:d
        mat1[j,k]=mat1[j,k] * Dam[j]/Dbp[k]
    end

    # mat2 = 1/(Dap) * Ua' * Ub * Dbm
    mat2 = adjoint(Ua) * Ub
    for j in 1:d, k in 1:d
        mat2[j,k]=mat2[j,k] * Dbm[k]/Dap[j]
    end
    
    #mat1 = mat1 + mat2
    mat1 = mat1 + mat2
    
    # decompose mat1: U, D, Vd
    U, D, Vd = decompose_udv!(mat1)

    # invert inner part: mat1 = (U D Vd)^(-1) = mat1^(-1)
    UDV_to_mat!(mat1, U, D, Vd, invert=true)

    # mat1 = 1/Dbp * mat1 /Dap
    for j in 1:d, k in 1:d
        mat1[j,k]=mat1[j,k] / Dbp[j] / Dap[k]
    end

    #mat1 = U D Vd
    U, D, Vd = decompose_udv!(mat1)

    # U = Vdb' * U , Vd = Vd * Ua'
    mul!(mat1, Vdb', U)
    mul!(mat2, Vd, Ua')
    U = mat1
    Vd = mat2

    return U, D, Vd
end


# Calculates (UaDaVda + UbDbVdb)^-1
# TODO: Optimize!
# I only made the function overwrite res. Otherwise it's unchanged compared to the one above.
function inv_sum_udvs!(mc, res, Ua, Da, Vda, Ub, Db, Vdb)
    
    d=length(Da)
    
    #DXp = max (X%D, 1) ,DXm = min (X%D, 1) and similarly for Y.
    Dap = max.(Da,1.)
    Dam = min.(Da,1.)
    Dbp = max.(Db,1.)
    Dbm = min.(Db,1.)

    #mat1= DXm * X%Vd * (Y%Vd)^dagger /DYp
    mat1 = Vda * adjoint(Vdb)
    for j in 1:d, k in 1:d
        mat1[j,k]=mat1[j,k] * Dam[j]/Dbp[k]
    end

    #mat2 = 1/(DXp) * (X%U)^dagger * Y%U * DYm
    mat2 = adjoint(Ua) * Ub
    for j in 1:d, k in 1:d
        mat2[j,k]=mat2[j,k] * Dbm[k]/Dap[j]
    end
    
    #mat1 = mat1+mat2
    mat1 = mat1 + mat2
    
    #invert mat1: mat1=mat1^(-1)
    U, D, Vd = decompose_udv!(mat1)
    UDV_to_mat!(mat1, U, D, Vd, invert=true)

    #mat1 = 1/DYp * mat1 /DXp
    for j in 1:d, k in 1:d
        mat1[j,k]=mat1[j,k] / Dbp[j] / Dap[k]
    end

    #mat1 = U D Vd
    U, D, Vd = decompose_udv!(mat1)

    # U = (Y%Vd)^dagger * U , Vd = Vd * (X%U)^dagger
    mul!(mat1,adjoint(Vdb),U)
    mul!(mat2,Vd,adjoint(Ua))
    U=mat1
    Vd=mat2

    # return U, D, Vd
    UDV_to_mat!(res, U, D, Vd, invert=false)
    nothing
end