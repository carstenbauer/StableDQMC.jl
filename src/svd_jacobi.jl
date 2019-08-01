"""
Calculates the SVD deomposition of a matrix by using the
Jacobi method. Returns a `SVD` factorization object.
"""
gesvj(x) = gesvj!(copy(x))
export gesvj
"""
Same as `gesvj` but saves space by overwriting the input matrix.
"""
function gesvj!(x)
  JacobiSVD.jsvd!(x)
end
export gesvj!
