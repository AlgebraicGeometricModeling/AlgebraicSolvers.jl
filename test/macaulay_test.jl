using DynamicPolynomials, AlgebraicSolvers, LinearAlgebra

X = @polyvar x y

P = [2-x-x*y-x^2, y^2+ x*y+x^2-x+y-1]

Mc = Macaulay()

R, L = res_matrix(Mc,P)
N, L = tnf(Mc,P)
B = quotient_basis(:macaulay,P)
M = mult_matrices(Mc,P)

Xi = solve(:macaulay, P; verbose=false)

Er = rel_error(P,Xi,X)
println("-- Rel error: ", norm(Er,Inf));

Xi
