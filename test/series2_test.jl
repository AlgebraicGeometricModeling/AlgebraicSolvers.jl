using DynamicPolynomials, LinearAlgebra, AlgebraicSolvers

X = @polyvar x y z
r  = 4
w0 = rand(r)
A0 = rand(3,r)

T0 = series(w0, A0, X, d)
w, A = decompose(T0, eps_rkf(1.e-10))

T = series(w, A, X, d)
@assert norm(T-T0)< 1.e-6
