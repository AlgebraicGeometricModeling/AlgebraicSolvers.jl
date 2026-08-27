using DynamicPolynomials, LinearAlgebra, AlgebraicSolvers

X = @polyvar x y
P = [x^2+x*y-2*x, y^2]

P = [2y + 2x + x*y + x^2, y^2]

Mc = Macaulay()

R, L = res_matrix(P, Mc)

N, IB = nullspace(R)

N, L, IB = tnf(P, Mc)

B = quot_basis(P, Mc)

M = mult_matrices(P, variables(P), Mc)

Xi, ms = AlgebraicSolvers.solve(P, Mc; verbose=true)

println("-- Sol mult: ",ms)



