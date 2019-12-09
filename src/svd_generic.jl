"""
Calculates the SVD deomposition of a matrix in a type generic manner.
Can for example be used with a `Matrix{BigFloat}`.
Returns a `SVD` factorization object.
"""
genericsvd(x) = genericsvd!(copy(x))
export genericsvd
"""
Same as `genericsvd` but saves space by overwriting the input matrix.
"""
function genericsvd!(x::AbstractArray{T}) where {T<:Union{BigFloat, BigInt}}
  GenericSVD.svd!(x; alg=nothing)
end
function genericsvd!(x::AbstractArray)
  GenericSVD.svd!(x)
end
export genericsvd!
