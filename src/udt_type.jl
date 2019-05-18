# UDT decomposition (basically QR)
struct UDT{E,Er<:Real,M<:AbstractArray{E}} <: Factorization{E}
    U::M
    D::Vector{Er}
    T::M
    function UDT{E,Er,M}(U, D, T) where {E,Er,M<:AbstractArray{E}}
        new{E,Er,M}(U, D, T)
    end
end

UDT(U::AbstractArray{E}, D::Vector{Er}, T::AbstractArray{E}) where {E,Er<:Real} = UDT{E,Er,typeof(U)}(U, D, T)
function UDT{E}(U::AbstractArray, D::AbstractVector{Er}, T::AbstractArray) where {E,Er<:Real}
    UDT(convert(AbstractArray{E}, U),
        convert(Vector{Er}, D),
        convert(AbstractArray{E}, T))
end


# Conversion
Base.AbstractMatrix(F::UDT) = (F.U * Diagonal(F.D)) * F.T
Base.AbstractArray(F::UDT) = AbstractMatrix(F)
Base.Matrix(F::UDT) = Array(AbstractArray(F))
Base.Array(F::UDT) = Matrix(F)


# iteration for destructuring into components
Base.iterate(S::UDT) = (S.U, Val(:D))
Base.iterate(S::UDT, ::Val{:D}) = (S.D, Val(:T))
Base.iterate(S::UDT, ::Val{:T}) = (S.T, Val(:done))
Base.iterate(S::UDT, ::Val{:done}) = nothing


Base.size(A::UDT, dim::Integer) = dim == 1 ? size(A.U, dim) : size(A.T, dim)
Base.size(A::UDT) = (size(A, 1), size(A, 2))

Base.similar(A::UDT) = UDT(similar(A.U), similar(A.D), similar(A.T))



# decomposition functions
"""
Compute the UDT decomposition of `A` and return an `UDT` object.

`U`, `D`, and `T`, can be obtained from the factorization `F`
with `F.U`, `F.D`, and `F.T` such that `A = U * Diagonal(D) * T`.

Iterating the decomposition produces the components `U`, `D`, and `V`.

Note that `T` is upper triangular only up to permutations of columns of `T`.
"""
function udt(A::AbstractMatrix{C}) where {C<:Number}
  F = qr(A, Val(true))
  _qr_to_udt(A, F)
end


"""
`udv!` is the same as `svd`, but saves space by overwriting the input `A`, instead of creating a
copy.
"""
function udt!(A::AbstractMatrix{C}) where {C<:Number}
  F = qr!(A, Val(true))
  _qr_to_udt(A, F)
end


@inline function _qr_to_udt(A::AbstractMatrix{C}, F::QRPivoted) where {C<:Number}
    n = size(A, 1)
    D = Vector{real(C)}(undef, n)
    R = F.R # F.R has regular matrix type
    @views F.p[F.p] = 1:n

    @inbounds for i in 1:n
    D[i] = abs(real(R[i,i]))
    end
    lmul!(Diagonal(1 ./ D), R)
    UDT(Matrix(F.Q), D, R[:, F.p])
end


function udt(x::Number)
    UDT(x == 0 ? fill(one(x), 1, 1) : fill(x/abs(x), 1, 1), [abs(x)], fill(one(x), 1, 1))
end
function udt(x::Integer)
    udt(float(x))
end


# operations
"""
    inv(F::UDT) -> AbstractMatrix

Computes the inverse matrix of the `UDT` decomposition of a matrix.
"""
function Base.inv(F::UDT)
    inv!(similar(F.U), F)
end


"""
    inv!(res, F::UDT) -> res

Same as `inv` but writes result into preallocated `res`.
"""
function inv!(res::M, F::UDT{E, Er, M}) where {E,Er,M}
    tmp = similar(F.U)
    ldiv!(tmp, lu(F.T), Diagonal(1 ./ F.D))
    mul!(res, tmp, F.U')
    return res
end

"""
    udt_mult(A::UDT, B::UDT) -> UDT

Stabilized multiplication of two `UDT` decompositions.
Returns a `UDT` factorization object.
"""
function udt_mult(A::UDT, B::UDT)
    mat = A.T * B.U
    lmul!(Diagonal(A.D), mat)
    rmul!(mat, Diagonal(B.D))
    F = udt!(mat)
    UDT(A.U * F.U, F.D, F.T * B.T)
end


"""
    *(A::UDT, B::UDT)

Stabilized multiplication of two `UDT` decompositions.
"""
function Base.:*(A::UDT, B::UDT)
    mat = A.T * B.U
    lmul!(Diagonal(A.D), mat)
    rmul!(mat, Diagonal(B.D))
    F = udt!(mat)
    (A.U * F.U) * Diagonal(F.D) * (F.T * B.T)
end