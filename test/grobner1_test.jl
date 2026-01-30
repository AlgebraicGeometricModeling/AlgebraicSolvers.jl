using AlgebraicSolvers, DynamicPolynomials, Groebner, LinearAlgebra

Mth = Grobner((P,X) -> Groebner.groebner(P, ordering = Groebner.DegRevLex(X)),
            (p,G) -> Groebner.normalform(p,G),
             G -> Groebner.quotient_basis(G)
            )

X = @polyvar x y

P = [-2+y-y^2+x^2*y, 1-3x+y+x*y^2]
P = [x^2+1, y^2-2]

B0 = quot_basis(Mth,P)
N, L = tnf(Mth,P)
M = mult_matrices(Mth,P)

Xi, ms, G, B  = AlgebraicSolvers.solve(Mth, P; verbose=true)

Er = rel_error(P,Xi)
println("-- Rel error: ", norm(Er,Inf));
println("-- Mult sols: ", ms);

 
