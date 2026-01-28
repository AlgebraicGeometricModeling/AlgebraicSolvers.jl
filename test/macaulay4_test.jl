using DynamicPolynomials, AlgebraicSolvers, LinearAlgebra, Groebner

DP = DynamicPolynomials

N=5
verbose = true

P = Groebner.Examples.katsuran(N)
P = convert_DP(P)

Mth = Macaulay()

X = DP.variables(P)

R, L = res_matrix(Mth, P)

N, _ = LinearAlgebra.nullspace(R)

B = quot_basis(Mth,P)

M = mult_matrices(Mth,P)

Xi  = AlgebraicSolvers.solve(Macaulay(),P,verbose=true)

Er = rel_error(P, Xi, X)
println("-- Rel error: ", norm(Er,Inf));
