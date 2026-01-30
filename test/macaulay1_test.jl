using DynamicPolynomials, AlgebraicSolvers, LinearAlgebra

X = @polyvar x y

P = [x*y-x^2, y^2+x*y-x+y]

Mc = Macaulay()

B= quot_basis(Mc,P)


N, L = tnf(Mc,P)
B = quot_basis(Mc,P)
M = mult_matrices(Mc,P)

Xi, ms = AlgebraicSolvers.solve(:macaulay, P; verbose=true)
Er = rel_error(P,Xi,X)
println("-- Rel error: ", norm(Er,Inf));

Xi

