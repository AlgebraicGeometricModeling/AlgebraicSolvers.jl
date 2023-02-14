using AlgebraicSolvers, LinearAlgebra

X = @Ring x y
n = length(X)

A0 = [one(X),x,x^2,y,x*y]
A1 = [one(X),x,y,x*y,x^2,x^2*y]
A2 = [one(X),x,x^2]


#= 
p1= rand(length(A1))'*A1
p2= rand(length(A2))'*A2
=#
p1= rand(length(A0))'*A0
p2= rand(length(A0))'*A0

P = [p1,p2]

Xi = solve_toric(P,X)

println("-- sol ", Xi)
Er = rel_error(P,Xi,X)
println("-- Rel error: ", norm(Er,Inf));

