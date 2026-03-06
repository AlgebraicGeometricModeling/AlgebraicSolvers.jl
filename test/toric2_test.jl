using DynamicPolynomials, AlgebraicSolvers

X = @polyvar x y
n = length(X)

A0 = DynamicPolynomials.monomials(1+x+x^2+y+x*y)
A1 = DynamicPolynomials.monomials(1+x+y+x*y+x^2+x^2*y)
A2 = DynamicPolynomials.monomials(1+x+x^2)

#= 
p1= rand(length(A1))'*A1
p2= rand(length(A2))'*A2
=#
p1= rand(length(A0))'*A0
p2= rand(length(A0))'*A0

P = [p1,p2]

Mth = Toric()
R, L  = res_matrix(P,Mth)
N, _ = LinearAlgebra.nullspace(R)
B =  quot_basis(P,Mth)

Xi, ms = AlgebraicSolvers.solve(P,Toric(); verbose=true)

#println("-- sol ", Xi)
Er = rel_error(P,Xi)
println("-- Rel error: ", norm(Er,Inf));
println("-- Mult sols: ", ms);
