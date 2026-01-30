using DynamicPolynomials, AlgebraicSolvers, LinearAlgebra, SparseArrays

X = @polyvar x0 x1 x2 x3 x4 

P=[ x0^2 + 2*x1^2 + 2*x2^2 + 2*x3^2 + 2*x4^2 - x0,
 2*x0*x1 + 2*x1*x2 + 2*x2*x3 + 2*x3*x4 - x1,
 x1^2 + 2*x0*x2 + 2*x1*x3 + 2*x2*x4 - x2,
 2*x1*x2 + 2*x0*x3 + 2*x1*x4 - x3,
 x0 + 2*x1 + 2*x2 + 2*x3 + 2*x4 - 1
]

R, L = res_matrix(Mc,P)
N, L = tnf(Mc,P)

B = quot_basis(:macaulay,P)

N, _ = LinearAlgebra.nullspace(R)

F = qr!(N')
U = F.R 
IB = column_basis(U)

M = mult_matrices(Mc,P)
Xi = eigdiag(M)

Xi, ms = AlgebraicSolvers.solve(Mc, P; verbose=true)


Er = rel_error(P, Xi, X)
println("-- Rel error: ", norm(Er,Inf));
println("-- Mult sols: ", ms);

#println("-- Sol: ", Xi)
