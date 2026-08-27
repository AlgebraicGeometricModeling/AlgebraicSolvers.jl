using DynamicPolynomials, AlgebraicSolvers

X = @polyvar x0 x1 x2
n = length(X)-1

d  = 5
M  = DynamicPolynomials.monomials(X,d)
s  = length(M)

P = (2*rand(n,s)-fill(1.0,n,s))*M

P = subs.(P,x0=>1)
Mth  = Macaulay()
R, L = res_matrix(P,Mth)
N, _ = LinearAlgebra.nullspace(R)

B    = quot_basis(P,Mth)
N, L = tnf(P, Mth)
    

Xi, ms = AlgebraicSolvers.solve(P,Mth;verbose=true)

println("-- Mult sols: ", ms);
