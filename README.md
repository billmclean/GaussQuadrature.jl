GaussQuadrature.jl
==================

Julia package to compute weights and points for the classical Gauss 
quadrature rules.

Handles the Legendre, Chebyshev, Jacobi, Laguerre and Hermite weights, 
and also provides the Lobatto and Radau variants of the Gauss rules.  
For example, to obtain a plain 5-point Gauss-Legendre rule do

    julia> using GaussQuadrature
    julia> x, w = legendre(5)

Read the initial comments in the src/GaussQuadrature.jl module
for full details, or read the help documentation for the individual
functions called legendre, chebyshev, jacobi, laguerre and hermite.  

Changelog: 

2013-10-12 
Initial release, version 0.1, tested with Julia 0.2 pre-release.

2013-11-01
Version 0.2, with support for arbitrary FloatingPoint arithmetic.
Works with Julia 0.3.

2015-09-25
Version 0.3, remove deprecated features for Julia 0.4, add 
documentation for individual functions.

2015-10-17
Version 0.3.1, add basic test for Travis, register package.
