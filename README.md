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
for full details.

Changelog: 

2013-10-12 
Initial release, version 0.1, tested with Julia 0.2 pre-release.

2013-11-01
Version 0.2, with support for arbitrary FloatingPoint arithmetic.
