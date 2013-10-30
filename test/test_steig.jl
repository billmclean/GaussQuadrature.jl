
using GaussQuadrature

d = [3.0, -2.0, 7.0, 6.0]
e = [-1.0, 0.0, 4.0, -2.0]
z = zero(d)

A = SymTridiagonal(d, e[1:end-1])
D, V = eig(A)

steig!(d, e, z, 30)
idx = sortperm(d)
d = d[idx]
z = z[idx]

@printf("\nEigenvalues:\n\n")
@printf("%20s  %20s  %12s\n\n", "eig", "steig!", "difference")
for i=1:length(d)
    @printf("%20.15f  %20.15f  %12.2e\n", d[i], D[i], d[i]-D[i])
end

@printf("\nFirst component of eigenvectors:\n\n")
@printf("%20s  %20s  %12s\n\n", "eig", "steig!", "difference")
for i=1:length(d)
    @printf("%20.15f  %20.15f  %12.2e\n", z[i], V[1,i], z[i]-V[1,i])
end
