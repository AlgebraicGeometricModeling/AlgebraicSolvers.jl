#include("CannyEmiris.jl")
#Main.CannyEmiris
using AlgebraicSolvers

#using SparseArrays

# With numeric coefficients
using DynamicPolynomials
X = @polyvar u v 
F = [u*v-v, 2u-4v+2,3*u*v+3*u-6.0]


NumCoeff = (F,X) -> function (i::Int64, E::Vector{Int64} )
    coefficient(F[i], prod(X.^E))
end

D0 = [[1,1] [1,1] [1,1]]
N0 = [1,1]

# With numeric coefficients
CE0, PM0 = CannyEmiris.Multihomogeneous(D0, N0, NumCoeff(F,X))

# With symbolic coefficients
function Coeff1(i, E)
    PolyVar{true}("a(" * string(i) * ")" * string(E) * "")+0;
end
CE1, PM1 = CannyEmiris.Multihomogeneous(D0, N0, Coeff1)


A = [[1,1] [1,1] [1,1]]
H = [[1,0] [1,1]]

# With random coefficients

G = randn(3,4)*[1,u,v,u*v]
CE2, PM2 = CannyEmiris.Zonotopes(A,H,NumCoeff(G,X))


# With symbolic coefficients
CE, PM = CannyEmiris.Zonotopes(A,H,Coeff1)


