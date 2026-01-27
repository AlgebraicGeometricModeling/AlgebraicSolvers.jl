using DynamicPolynomials, AlgebraicSolvers, LinearAlgebra

X = @polyvar x y

P = [1-3y+x+2x*y,x*y-x+y-1]
Tr = Toric()

R, L = res_matrix(Tr,P)
N, L = tnf(Tr,P)
B = quot_basis(Tr,P)
M = mult_matrices(Tr,P,X)

Xi = solve(Tr, P; verbose=true)
Er = rel_error(P,Xi,X)
println("-- Rel error: ", norm(Er,Inf));
Xi
