include("modified_Chebyshev.jl")
n = 5
a, b = shifted_legendre_coefs(Float64, n)
ρ = 3.5
ν = modified_moments(n, ρ)

import GaussQuadrature as GQ
ν_old = GQ.modified_moments(n, ρ)

using Polynomials
using QuadGK

P = OffsetVector{Polynomial}(undef, 0:6)
P[0] = Polynomial([1])
P[1] = Polynomial([0,1])
for j = 2:5
    P[j] = Polynomial( (2j-1) * P[1] * P[j-1] - (j-1)*P[j-2] ) / j
end
ν_bf = OffsetVector{Float64}(undef,0:5)
err = OffsetVector{Float64}(undef,0:5)
for j = 0:5
    ν_bf[j], err[j] = quadgk(x -> P[j](2x-1) * x^ρ * log(1/x) * sqrt(2j+1), 0, 1)
end

