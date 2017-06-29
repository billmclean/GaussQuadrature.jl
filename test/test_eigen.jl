
const n = 8
a, b = legendre_coefs(Float64, n)

Tn = SymTridiagonal(copy(a), b[2:n])
F = eigfact(Tn)

w = zeros(n)
special_eigenproblem!(a, b, w, 30)

idx = sortperm(a)
@test a[idx] ≈ F[:values]
@test abs.(w[idx]) ≈ abs.(F[:vectors][1,:])
