import OffsetArrays: OffsetVector, OffsetMatrix

function shifted_legendre_coefs(::Type{T}, n) where T <: AbstractFloat
    a = fill(1/convert(T, 2), n)
    b = OffsetVector{T}(undef, 0:n)
    b[0] = one(T)
    for k = 1:n
	rk = convert(T, k)
	b[k] = rk / ( 2 * sqrt((2rk-1)*(2rk+1)) )
    end
    return a, b
end

function modified_moments(n::Integer, ρ::T) where T <: AbstractFloat
    r = (ρ<0) ? 0 : round(Integer, ρ)
    m = min(2n+1, r)
    S = 1 / (1 + ρ)
    B = S
    ν = OffsetVector{T}(undef, 0:2n+1)
    ν[0] = B * S
    for j = 1:m
	B *= (ρ + 1 - j) / (ρ + 1 + j)
	S -= 2j / ( (ρ + 1 - j)*(ρ + 1 + j) )
	tjp1 = convert(T, 2j+1)
	ν[j] = B * S * sqrt(tjp1)
        println("j = $j, B = $B, S = $S")
    end
    if 2n+1 > r
	j = r + 1
	X = 2(r+1) / (ρ + r + 2)^2
	tjp1 = convert(T, 2j+1)
	ν[j] = B * ( S * (ρ - r)/(ρ + r + 2) - X ) * sqrt(tjp1)
	for j = r+2:2n+1
	    B *= (ρ + 1 - j) / (ρ + 1 + j)
	    S -= 2j / ( (ρ + 1 - j)*(ρ + 1 + j) )
	    tjp1 = convert(T, 2j+1)
	    ν[j] = B * ( S * (ρ - r)/(ρ + r + 2) - X ) * sqrt(tjp1)
	end
    end
    return ν
end

function modified_moments(::Type{T}, n::Integer, 
                          r::Integer) where T <: AbstractFloat
    m = min(2n, r+1)
    rr = convert(T, r)
    S = 1 / (1 + rr)
    B = S
    ν = OffsetVector{T}(undef, 0:2n+1)
    ν[0] = B * S
    for j = 1:m
        B *= (rr + 1 - j) / (rr + 1 + j)
        S -= 2j / ( (rr + 1 - j) * (rr + 1 + j) )
        tjp1 = convert(T, 2j+1)
        ν[j] = B * S * sqrt(tjp1)
    end
    if 2n+1 > r
        product = one(T)
        for k = 1:r
            product *= convert(T, k) / (2 * (2k+1))
        end
        sq = sqrt(2rr + 3)
        ν[r+1] = - ( product / (2*(r+1)) ) * sq
        for j = r+2:2n-1
            tjp1 = convert(T, 2j+1) 
            new_sq = sqrt(tjp1)
            ν[j] = -ν[j-1] * (neq_sq/sq) * ( (j - rr - 1)/(j + rr + 1) )
            sq = new_sq
        end
    end
    return ν
end
    
end

function modified_Chebyshev(a::Vector{T}, b::OffsetVector{T}, 
	                    ν::OffsetVector{T}) where T <: AbstractFloat
    n = lastindex(a)
    σ = OffsetMatrix{T}(undef, 0:2n+1, 0:n)
    fill!(σ, zero(T))
    α = Vector{T}(undef, n)
    β = OffsetVector{T}(undef, 0:n)
    σ[1:2n-1,1] .= ν[1:2n-1]
    α[1] = a[1] + b[1] * σ[2,1] / σ[1,1]
    β[0] = sqrt(b[0] * ν[0])
    t = b[2]*σ[3,1] + (a[2] - α[1]) * σ[2,1] - β[0] * σ[2,0]
    β[1] = b[1] * sqrt( 1 + t / ( b[1] * σ[1,1] ) )
    for k = 2:n
	for j = k:2n+2-k
	    σ[j,k] = ( b[j] * σ[j+1,k-1] + (a[j] - α[k-1]) * σ[j,k-1]
                     + b[j-1] * σ[j-1,k-1] - β[k-2] * σ[j,k-2] ) / β[k-1]
	end
	α[k] = a[k] + ( b[k] * σ[k+1,k] - β[k-1] * σ[k,k-1] ) / σ[k,k]
	β[k] = b[k] * sqrt( 1 + t / ( b[k] * σ[k,k] ) )
    end
    return α, β, σ
end
