using GaussQuadrature
using Test
import LinearAlgebra: SymTridiagonal, eigen
import SpecialFunctions: besseli

include("test_eigenproblem1.jl")
include("test_eigenproblem2.jl")
include("test_legendre.jl")
include("test_chebyshev.jl")
include("test_jacobi.jl")
include("test_laguerre.jl")
include("test_hermite.jl")
include("test_logweight.jl")
