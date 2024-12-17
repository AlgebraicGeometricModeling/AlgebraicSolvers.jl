using DynamicPolynomials, AlgebraicSolvers

#using SparseArrays
# With numeric coefficients

X = @polyvar u v 
F = [u*v-v, 2u-4v+2,3*u*v+3*u-6.0]



D0 = [[1,1] [1,1] [1,1]]
N0 = [1,1]

# With numeric coefficients
CE0, PM0 = CannyEmiris.Multihomogeneous(D0, N0,CannyEmiris.NumCoeff(F,X))

# With symbolic coefficients
function Coeff1(i, E)
    DynamicPolynomials.Variable{DynamicPolynomials.Commutative{DynamicPolynomials.CreationOrder}, Graded{LexOrder}}("a(" * string(i) * ")" * string(E) * "")+0;
end

CE1, PM1 = CannyEmiris.Multihomogeneous(D0, N0, Coeff1)


A = [[1,1] [1,1] [1,1]]
H = [[1,0] [1,1]]

# With random coefficients

G = randn(3,4)*[1,u,v,u*v]
CE2, PM2 = CannyEmiris.Zonotopes(A,H,CannyEmiris.NumCoeff(G,X))


# With symbolic coefficients
CE, PM = CannyEmiris.Zonotopes(A,H,Coeff1)


