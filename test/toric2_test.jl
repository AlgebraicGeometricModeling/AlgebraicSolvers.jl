using DynamicPolynomials, AlgebraicSolvers

X = @polyvar x y z
n = length(X)

A =  monomials(1+x+y+z+x*y+x*z+y*z+x*y*z)

#= 
p1= rand(length(A1))'*A1
p2= rand(length(A2))'*A2
=#

p1= rand(length(A))'*A
p2= rand(length(A))'*A
p3= rand(length(A))'*A

P = [p1,p2,p3]

Xi = solve_toric(P; verbose=true)

#println("-- sol ", Xi)
Er = rel_error(P,Xi,X)
println("-- Rel error: ", norm(Er,Inf));

