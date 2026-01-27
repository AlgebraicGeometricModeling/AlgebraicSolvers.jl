using DynamicPolynomials, AlgebraicSolvers, LinearAlgebra

X = @polyvar x y
P = [x^2-2, y^2-1]

Mc = Macaulay()

R, L = res_matrix(Mc,P)

N, IB = nullspace(R)

N, L, IB = tnf(Mc,P)

B = quot_basis(Mc, P)

M = mult_matrices(Mc,P)


Xi = AlgebraicSolvers.solve(Mc, P; verbose=true)
Er = rel_error(P,Xi,X)
println("-- Rel error: ", norm(Er,Inf));
println("-- sol: ",Xi)

@polyvar z
Ph = [x^2-2*z^2, y^3-z^3]

Xih = AlgebraicSolvers.solve(Mc, Ph; verbose=true)
Erh = rel_error(Ph,Xih,[x,y,z])
println("-- Rel error: ", norm(Erh,Inf));
println("-- sol: ",Xih)

