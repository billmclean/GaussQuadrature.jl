using Polynomials
using QuadGK

let
    n = 5
    P = OffsetVector{Polynomial}(undef, 0:n)
    P[0] = Polynomial([1])
    P[1] = Polynomial([0, 1])
    for j = 2:n
        P[j] = Polynomial( (2j-1) * P[1] * P[j-1] - (j-1)*P[j-2] ) / j
    end

    ρ = 2.5
    ν = modified_moments(n, ρ)
    for j = 0:n
        correct, err = quadgk(0, 1) do x
            return P[j](2x-1) * x^ρ * log(1/x) * sqrt(2j+1)
        end
        @ test ν[j] ≈ correct
    end
end
    
