using DynamicPolynomials, AlgebraicSolvers, LinearAlgebra

X = @polyvar x y

P = [x*y-x^2, y^2+x*y-x+y]

Mc = Macaulay()

B= quot_basis(P,Mc)


N, L = tnf(P,Mc)
B = quot_basis(P,Mc)
M = mult_matrices(P,variables(P), Mc)


Xi, E, Info = diagonalization(M)

Er = rel_error(P,Xi,X)
println("-- Rel error: ", norm(Er,Inf));


