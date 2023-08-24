
let
    a = Float64[3,  -2, 6,  0, 4]
    b = Float64[NaN, 1, 8, -5, 3]
    A = SymTridiagonal(a, b[2:end])
    F = eigen(A)
    correct = F.values
    z = Vector{Float64}(undef, 5)
    maxits = 10
    special_eigenproblem!(a, b, z, maxits)
    idx = sortperm(a)
    @test a[idx] ≈ correct
    @test abs.(z[idx]) ≈ abs.(F.vectors[1,:])
end
