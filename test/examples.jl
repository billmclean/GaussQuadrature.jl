using GaussQuadrature

function variant(endpt)
    if endpt == neither
        return "Default    "
    elseif endpt == left
        return "Left Radau "
    elseif endpt == right
        return "Right Radau"
    elseif endpt == both
        return "Lobatto    "
    else
        error("Unknown endpt")
    end
end

function table(f, rule, name, ans, n, endpts)
    println("\nTesting ", name, " rule with ", n, " points:")
    for endpt in endpts
        x, w = rule(n, endpt)
        integral = 0
        for j = 1:n
            integral += w[j] * f(x[j])
        end
        relerr = ( integral - ans ) / ans
        @printf("\t%s  %12.2e\n", variant(endpt), relerr)
    end
end

Beta(x,y) = gamma(x) * gamma(y) / gamma(x+y)

jacobiintegral(alpha, beta, c) = 2^(alpha+beta+1)*Beta(1+alpha,1+beta) / ( 
                                (c-1)^(1+alpha) * (c+1)^(1+beta) )
function laguerrfunc(x, alpha)
    if x > eps(typeof(x))
        r = -expm1(-x) / x
    else
        r = 1.0
    end
    return r^alpha
end

laguerreintegral(alpha) = Beta(1, 1+alpha)

hermiteintegral(a) = sqrt(pi) * exp(a^2)

endpts = [neither, left, right, both]
table(x -> 1/(1+x^2), legendre, "Legendre", pi/2, 18, endpts)

c = 3.0
table(x -> 1.0/(x+c), 
      (n, endpt) -> chebyshev(n, 1, endpt), "Chebyshev (kind=1)",
      jacobiintegral(-0.5, -0.5, c), 18, endpts)

table(x -> 1.0/(x+c)^3.0, 
      (n, endpt) -> chebyshev(n, 2, endpt), "Chebyshev (kind=2)",
      jacobiintegral(0.5, 0.5, c), 18, endpts)

alpha = 0.4
beta = -0.2
table(x -> 1.0/(x+c)^(alpha+beta+2), 
      (n, endpt) -> jacobi(n, alpha, beta, endpt), "Jacobi", 
      jacobiintegral(alpha, beta, c), 18, endpts)

table(x -> laguerrfunc(x, alpha), 
      (n, endpt) -> laguerre(n, alpha, endpt), "Laguerre",
      laguerreintegral(alpha), 20, [neither, left])

a = 1.2
table(x -> exp(2*a*x), (n, endpt) -> hermite(n), "Hermite",
      hermiteintegral(a), 18, [neither])

