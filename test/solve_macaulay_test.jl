using DynamicPolynomials, AlgebraicSolvers

X = @polyvar x y

P = [2-x-x*y-x^2, y^2+ x*y+x^2-x+y-1]

Xi = solve_macaulay(P,X; verbose=false)

Er = rel_error(P,Xi,X)
println("-- Rel error: ", norm(Er,Inf));
