using DynamicPolynomials, AlgebraicSolvers

X = @polyvar x0 x1 x2
n = length(X)-1

d  = 5
M  = DynamicPolynomials.monomials(X,d)
s  = length(M)

P = (2*rand(n,s)-fill(1.0,n,s))*M

Mth  = Macaulay()
R, L = res_matrix(Mth,P)
N, _ = LinearAlgebra.nullspace(R)

B    = quot_basis(Mth,P)
N, L = tnf(Mth, P)
    

Xi = AlgebraicSolvers.solve(Mth,P;verbose=true)

Er = rel_error(P,Xi)
println("-- Rel. error: ", norm(Er,Inf))
