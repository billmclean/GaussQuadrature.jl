using GaussQuadrature

T = BigFloat

print("\nFloating point data type is ", T, '\n')
@printf("\teps = %0.2e\n\n", eps(T))

npts = { Float32 => { "Legendre"  => 10, 
                      "Chebyshev" => 10,
                      "Jacobi"    => 10,
                      "Laguerre"  => 10,
                      "Hermite"   => 10 },
         Float64 => { "Legendre"  => 20, 
                      "Chebyshev" => 12,
                      "Jacobi"    => 12,
                      "Laguerre"  => 20,
                      "Hermite"   => 16 },
         BigFloat => { "Legendre" => 100, 
                      "Chebyshev" => 60,
                      "Jacobi"    => 60,
                      "Laguerre"  => 600,
                      "Hermite"   => 80 }
}

const half = one(T) / 2

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

function table{T}(f, rule, name, ans::T, n::Integer, endpts)
    println("\nTesting ", name, " rule with ", n, " points:")
    for endpt in endpts
        x, w = rule(n, endpt)
        integral = zero(typeof(x[1])) 
        for j = 1:n
            integral += w[j] * f(x[j])
        end
        relerr = ( integral - ans ) / ans
        @printf("\t%s  %12.2e\n", variant(endpt), relerr)
    end
end

Beta{T}(x::T, y::T) = gamma(x) * gamma(y) / gamma(x+y)

legendrefunc{T}(x::T) = one(T) / ( one(T) + x^2 )

jacobifunc{T}(x::T, c::T, alpha::T, beta::T) = (
             (x+c)^convert(T,-alpha-beta-2) )

jacobiintegral{T}(alpha::T, beta::T, c::T) = ( 2^(alpha+beta+1)
      * Beta(1+alpha,1+beta) / ( (c-1)^(1+alpha) * (c+1)^(1+beta) ) )

function laguerrefunc{T}(x::T, alpha::T)
    if x > eps(T)
        r = -expm1(-x) / x
    else
        r = one(T)
    end
    return r^alpha
end

laguerreintegral{T}(alpha::T) = Beta(one(T), 1+alpha)

hermiteintegral{T}(a::T) = sqrt(convert(T, pi)) * exp(a^2)

endpts = [neither, left, right, both]
table(legendrefunc, (n, endpt) -> legendre(T, n, endpt), 
      "Legendre", convert(T, pi)/2, npts[T]["Legendre"], endpts)

c = convert(T, 3.0)
table(x -> jacobifunc(x, c, -half, -half),
      (n, endpt) -> chebyshev(T, n, 1, endpt), "Chebyshev (kind=1)",
      jacobiintegral(-half, -half, c), npts[T]["Chebyshev"], endpts)

table(x -> jacobifunc(x, c, half, half),
      (n, endpt) -> chebyshev(T, n, 2, endpt), "Chebyshev (kind=2)",
      jacobiintegral(half, half, c), npts[T]["Chebyshev"], endpts)

alpha = convert(T, 2)  / 5
beta  = convert(T, -1) / 5 
table(x -> jacobifunc(x, c, alpha, beta),
      (n, endpt) -> jacobi(n, alpha, beta, endpt), "Jacobi", 
      jacobiintegral(alpha, beta, c), npts[T]["Jacobi"], endpts)

table(x -> laguerrefunc(x, alpha), 
      (n, endpt) -> laguerre(n, alpha, endpt), "Laguerre",
      laguerreintegral(alpha), npts[T]["Laguerre"], [neither, left])

a = convert(T, 6) / 5
table(x -> exp(2*a*x), (n, endpt) -> hermite(T, n), "Hermite",
      hermiteintegral(a), npts[T]["Hermite"], [neither])

