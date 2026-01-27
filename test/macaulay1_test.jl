using DynamicPolynomials, AlgebraicSolvers, LinearAlgebra

X = @polyvar x y

P = [2-x-x*y-x^2, y^2+ x*y+x^2-x+y-3]

Mc = Macaulay()

B= quot_basis(Mc,P)


N, L = tnf(Mc,P)
B = quot_basis(Mc,P)
M = mult_matrices(Mc,P)

Xi = AlgebraicSolvers.solve(:macaulay, P; verbose=true)
Er = rel_error(P,Xi,X)
println("-- Rel error: ", norm(Er,Inf));

Xi

