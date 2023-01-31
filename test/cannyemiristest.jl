#include("CannyEmiris.jl")
#Main.CannyEmiris
using AlgebraicSolvers



using DynamicPolynomials
X = @polyvar u v 
F = [u*v-v+1.0,2u+4v-2,3*u*v+3*u-3.0]


function Coeff0(i::Int64, E::Vector{Int64} )
    coefficient(F[i], prod(X.^E))
end

D0 = [[1,1] [1,1] [1,1]]
N0 = [1,1]

CE0, PM0 = CannyEmiris.Multihomogeneous(D0, N0, Coeff0)

A = [[1,1] [1,1] [1,1]]
H = [[1,0] [0,1]]

CE, PM = CannyEmiris.Zonotopes(A,H)


