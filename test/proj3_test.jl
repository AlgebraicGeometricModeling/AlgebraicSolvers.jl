using DynamicPolynomials, AlgebraicSolvers

X = @polyvar x0 x1 x2
n = length(X)-1

d  = 5
M  = monomials(X,d)
s  = length(M)

P = (2*rand(n,s)-fill(1.0,n,s))*M

Xi = solve_macaulay(P;verbose=true)
Er = rel_error(P,Xi)
println("-- Rel. error: ", norm(Er,Inf))
