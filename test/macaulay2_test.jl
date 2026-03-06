using DynamicPolynomials, LinearAlgebra, AlgebraicSolvers

X = @polyvar x y
P = [x^2+x*y-2*x, y^2]

P = [2y + 2x + x*y + x^2, y^2]

Mc = Macaulay()

R, L = res_matrix(P, Mc)

N, IB = nullspace(R)

N, L, IB = tnf(P, Mc)

B = quot_basis(P, Mc)

M = mult_matrices(P, variables(P), Mc)

Xi, ms = AlgebraicSolvers.solve(P, Mc; verbose=true)

Er = rel_error(P,Xi,X)
println("-- Rel error: ", norm(Er,Inf));
println("-- Sol mult: ",ms)

@polyvar z
Ph = [x^2-2*z^2, y^3-z^3]

Xih, msh = AlgebraicSolvers.solve(Ph, Mc; verbose=true)

Erh = rel_error(Ph,Xih)
println("-- Rel error: ", norm(Erh,Inf));
println("-- Mult sols: ", msh)


