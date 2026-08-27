using DynamicPolynomials, LinearAlgebra, AlgebraicSolvers, Groebner

DP = DynamicPolynomials

N=5
verbose = true

P = Groebner.Examples.katsuran(N)
P = convert_DP(P)

Mth = Macaulay()

X = DP.variables(P)

R, L = res_matrix(P, Mth)

N, Ib = LinearAlgebra.nullspace(R)

B = quot_basis(P, Mth)

M = mult_matrices(P, variables(P), Mth)

#ms = multiplicities(v)

Xi, ms  = AlgebraicSolvers.solve(P, Macaulay();verbose=true)

println("-- Sol mult: ", ms);

