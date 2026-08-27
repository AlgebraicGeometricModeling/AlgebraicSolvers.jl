using DynamicPolynomials, AlgebraicSolvers, LinearAlgebra

X = @polyvar x y

P = [1-3*x+y+2x*y,x^2+x+y-1]
Tr = Toric()

R, L = res_matrix(P, Tr)
N, L = tnf(P, Tr)
B = quot_basis(P, Tr)
M = mult_matrices(P,X,Tr)

Xi, ms = AlgebraicSolvers.solve(P, Tr; verbose=true)

println("-- Mult sols: ", ms);
