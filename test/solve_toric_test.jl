using DynamicPolynomials, AlgebraicSolvers

X = @polyvar x y

P = [1-3y+2x*y,x*y-x+y-1]

R, L = toric_matrix(P) 

Xi = solve(Toric(), P; verbose=true)

Er = rel_error(P,Xi,X)
println("-- Rel error: ", norm(Er,Inf));

Xi
