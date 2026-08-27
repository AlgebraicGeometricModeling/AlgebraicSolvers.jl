using DynamicPolynomials, AlgebraicSolvers, LinearAlgebra

X = @polyvar x y

P = [x*y-x^2, y^2+x*y-x+y]

Mc = Macaulay()

B= quot_basis(P,Mc)


N, L = tnf(P,Mc)
B = quot_basis(P,Mc)
M = mult_matrices(P,variables(P), Mc)

Xi, ms = AlgebraicSolvers.solve(P, Mc; verbose=true)

println("-- Mult sols: ", ms);


