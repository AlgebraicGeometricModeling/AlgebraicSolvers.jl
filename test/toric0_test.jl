using DynamicPolynomials, AlgebraicSolvers, LinearAlgebra

X = @polyvar x y

P = [1-3y+2x*y,x*y-x+y-1]
T = Toric()

R, L = res_matrix(T,P)
N, L = tnf(T,P)
B = quo_basis(T,P)
M = mult_matrices(T,P,X)

Xi = solve(T, P; verbose=true)

Er = rel_error(P,Xi,X)
println("-- Rel error: ", norm(Er,Inf));

Xi
