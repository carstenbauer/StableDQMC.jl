gesdd = x -> LinearAlgebra.LAPACK.gesdd!('A', copy(x))
gesvd = x -> LinearAlgebra.LAPACK.gesvd!('A', 'A', copy(x))
genericsvd = x -> (F = svd(x); return (F.U, F.S, F.Vt))
gesvj = x -> JacobiSVD.gesvj!('G','U','V', copy(x))


decompose_udv!(A::AbstractMatrix{<:Number}) = LinearAlgebra.LAPACK.gesvd!('A','A',A)
decompose_udv(A::AbstractMatrix{T}) where T<:Number = decompose_udv!(copy(A))



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