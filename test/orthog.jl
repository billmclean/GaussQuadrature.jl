# Orthogonality tests

using GaussQuadrature

T = BigFloat
xlegendre(n, alpha, beta, endpt) = legendre(T, n, endpt)
xchebyshev_first(n, alpha, beta, endpt) = chebyshev(T, n, 1, endpt)
xchebyshev_second(n, alpha, beta, endpt) = chebyshev(T, n, 2, endpt)
xhermite(n, alpha, beta, endpt) = hermite(T, n)
xlaguerre(n, alpha, beta, endpt) = laguerre(n, alpha, endpt)

xlegendre_coeff(n, alpha, beta, endpt) = legendre_coeff(T, n, endpt)
xchebyshev_first_coeff(n, alpha, beta, endpt) = chebyshev_coeff(T, 
                                                n, 1, endpt)
xchebyshev_second_coeff(n, alpha, beta, endpt) = chebyshev_coeff(T, 
                                                 n, 2, endpt)
xhermite_coeff(n, alpha, beta, endpt) = hermite_coeff(T, n)
xlaguerre_coeff(n, alpha, beta, endpt) = laguerre_coeff(n, alpha, endpt)

rule = [xlegendre, xchebyshev_first, xchebyshev_second, 
        xhermite, jacobi, xlaguerre]

coeff = [ xlegendre_coeff, xchebyshev_first_coeff, 
          xchebyshev_second_coeff, xhermite_coeff, jacobi_coeff, 
          xlaguerre_coeff ]

name = ["Legendre", "Chebyshev (first)", "Chebyshev (second)",
         "Hermite", "Jacobi", "Laguerre" ]

alpha = one(T)
beta  = one(T)/2
xargs = [(),(1),(2), (), (alpha, beta), (alpha)]

function discrepancy(n, dop, w, p)
    d = 0.0
    for k = 0:dop
        for j = 0:min(k-1, dop-k)
            s = 0.0
            for l = 1:n
                s += w[l] * p[l,j+1] * p[l,k+1]
            end
            d = max(d, abs(s))
        end 
        if 2*k <= dop
            s = 0.0
            for l = 1:n
                s += w[l] * p[l,k+1]^2
            end
            d = max(d, abs(s-1.0))
        end
    end 
    return d
end

function dop(n, endpt)
    if endpt == neither
        return 2*n - 1
    elseif endpt == left || endpt == right
        return 2*n - 2
    else
        return 2*n - 3
    end
end

function test_rule(descr, nmin, nmax, rule, coeff, name, endpt)
    println("\nTesting ", descr, " rules with ", nmin, " to ", 
            nmax, " points")
    for case in zip(rule, coeff, name)
        rulefunc, coeffunc, rulename = case
        maxd = 0.0
        for n = nmin:nmax
            x, w = rulefunc(n, alpha, beta, endpt)
            a, b, muzero = coeffunc(2n, alpha, beta, endpt)
            p = orthonormal_poly(x, a, b, muzero)
            d = discrepancy(n, dop(n, endpt), w, p)
            maxd = max(d, maxd)
        end
        @printf("\t%s: %12.4e\n", rpad(rulename, 20, ' '), maxd)
    end
end

print("\nFloating point type data type is ", T, '\n')

test_rule("plain Gauss", 1, 5, rule, coeff, name, neither)

filter!(x->(x!=xhermite), rule)
filter!(x->(x!=xhermite_coeff), coeff)
filter!(x->(x!="Hermite"), name)

test_rule("left Radau", 1, 5, rule, coeff, name, left)

filter!(x->(x!=xlaguerre), rule)
filter!(x->(x!=xlaguerre_coeff), coeff)
filter!(x->(x!="Laguerre"), name)

test_rule("right Radau", 1, 5, rule, coeff, name, right)
test_rule("Lobatto", 2, 5, rule, coeff, name, both)
