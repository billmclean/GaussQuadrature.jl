"""
Finds the eigenvalues and first components of the normalised
eigenvectors of a symmetric tridiagonal matrix by the implicit
QR method.

a[i]   On entry, holds the ith diagonal entry of the matrix.
       On exit, holds the ith eigenvalue.

b[i]   On entry, holds the [i,i-1] entry of the matrix for
       i = 2, 3, ..., n.  (The value of e[1] is not used.)
       On exit, e is overwritten.

z[i]   On exit, holds the first component of the ith normalised
       eigenvector associated with d[i].

maxits The maximum number of QR iterations per eigenvalue.

Martin and Wilkinson, Numer. Math. 12: 377-383 (1968).
Dubrulle, Numer. Math. 15: 450 (1970).
Handbook for Automatic Computation, Vol ii, Linear Algebra, pp. 241-248, 1971.
"""
function special_eigenproblem!(a::Vector{T}, b::Vector{T}, z::Vector{T},
	                       maxits::Integer) where T <: AbstractFloat
    n = lastindex(a)
    fill!(z, zero(T))
    z[1] = one(T)
    if n == 1 # Nothing to do for 1x1 matrix
	return
    end
    for k = n:-1:3
	converged = false
	for p = 1:maxits
	    if abs(b[k]) < eps(T) * ( abs(a[k]) + abs(a[k-1]) )
		converged = true
		break
	    end
	    # Wilkinson shift
	    λ₊, λ₋ = eig2x2(a, b, k-1) 
	    μ = ( abs(λ₊ - a[k]) < abs(λ₋ - a[k]) ) ? λ₊ : λ₋
	    implicit_Q_step!(a, b, z, μ, k)
	end
	if !(converged)
	    error("""Eigenvalue λ[$k] = $(a[k]) failed to converge after
		     $maxits iterations.  Try increasing maxits""")
	end
    end
    λ₊, λ₋ = eig2x2!(a, b, z)
    a[1] = λ₊
    a[2] = λ₋
    b[2] = zero(T)
end

function implicit_Q_step!(a::Vector{T}, b::Vector{T}, z::Vector{T}, μ::T, 
	                  k::Integer) where T <: AbstractFloat
    c, s, d = initial_Givens_step!(a, b, μ)
    Givens_update!(z, c, s, 1)
    for j = 2:k-2
	c, s, d = general_Givens_step!(a, b, d, j)
        Givens_update!(z, c, s, j)
    end
    c, s = final_Givens_step!(a, b, d, k-1)
    Givens_update!(z, c, s, k-1)
end

function Givens_update!(z::Vector{T}, c::T, s::T, 
	                j::Integer) where T <: AbstractFloat
    t1, t2 = z[j], z[j+1]
    z[j]   =  c * t1 + s * t2
    z[j+1] = -s * t1 + c * t2
end

function initial_Givens_step!(a::Vector{T}, b::Vector{T}, μ::T
	                     ) where T <: AbstractFloat
    h = hypot(a[1] - μ, b[2])
    c = (a[1] - μ) / h
    s = b[2] / h
    t1, t2 = a[1], a[2]
    a[1] = c^2 * t1 + 2c * s * b[2] + s^2 * t2
    a[2] = s^2 * t1 - 2c * s * b[2] + c^2 * t2
    b[2] = (c+s) * (c-s) * b[2] + c * s * (t2 - t1)
    d = s * b[3]
    b[3] *= c
    return c, s, d
end

function general_Givens_step!(a::Vector{T}, b::Vector{T}, d::T, 
	                      j::Integer) where T <: AbstractFloat
    h = hypot(b[j], d)
    c = b[j] / h
    s = d / h
    t1, t2 = a[j], a[j+1]
    a[j]   = c^2 * t1 + 2c * s * b[j+1] + s^2 * t2
    a[j+1] = s^2 * t1 - 2c * s * b[j+1] + c^2 * t2
    b[j]   = c * b[j] + s * d
    b[j+1] = (c+s) * (c-s) * b[j+1] + c * s * (t2 - t1)
    d = s * b[j+2]
    b[j+2] *= c
    return c, s, d
end

function final_Givens_step!(a::Vector{T}, b::Vector{T}, d::T, 
	                    j::Integer) where T <: AbstractFloat
    h = hypot(b[j], d)
    c = b[j] / h
    s = d / h
    t1, t2 = a[j], a[j+1]
    a[j]   = c^2 * t1 + 2c * s * b[j+1] + s^2 * t2
    a[j+1] = s^2 * t1 - 2c * s * b[j+1] + c^2 * t2
    b[j]   = c * b[j] + s * d
    b[j+1] = (c+s) * (c-s) * b[j+1] + c * s * (t2 - t1)
    return c, s, d
end

"""
    eig2x2(a, b, j) -> λ₊, λ₋

Compute the eigenvalues of the symmetric 2x2 matrix

    a[j]     b[j+1]
    b[j+1]   a[j+1].

If the optional keyword argument `vals_only` is set to `false`, then
the function returns λ₊, λ₋, v₊, v₋ where v₊ and v₋ are the normalised
eigenvectors.
"""
function eig2x2(a::Vector{T}, b::Vector{T}, j::Int64;
                vals_only=true) where T <: AbstractFloat
    a1, a2, b1 = a[j], a[j+1], b[j+1]
    S = a1 + a2
    D = hypot(a1-a2, 2b1)
    λ₊ = ( S + D ) / 2
    λ₋ = ( S - D ) / 2
    if vals_only
        return λ₊, λ₋
    else
        v₊ = [2b1, a2-a1+D]
        v₋ = [-v₊[2], 2b1]
        normv = norm(v₊)
        v₊ /= normv
        v₋ /= normv
        return λ₊, λ₋, v₊, v₋
    end
end

function eig2x2!(a::Vector{T}, b::Vector{T}, z::Vector{T}
                 ) where T <: AbstractFloat
    λ₊, λ₋, v₊, v₋ = eig2x2(a, b, 1, vals_only=false)
    t1, t2 = z[1], z[2]
    z[1] = v₊[1] * t1 + v₊[2] * t2
    z[2] = v₋[1] * t1 + v₋[2] * t2
    return λ₊, λ₋
end


