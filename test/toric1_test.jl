using DynamicPolynomials, AlgebraicSolvers, LinearAlgebra

X = @polyvar x y

P = [1-3*x+y+2x*y,x^2+x+y-1]
Tr = Toric()

R, L = res_matrix(Tr,P)
N, L = tnf(Tr,P)
B = quot_basis(Tr,P)
M = mult_matrices(Tr,P,X)

Xi, ms = AlgebraicSolvers.solve(Tr, P; verbose=true)
Er = rel_error(P,Xi)
println("-- Rel error: ", norm(Er,Inf));
println("-- Mult sols: ", ms);
