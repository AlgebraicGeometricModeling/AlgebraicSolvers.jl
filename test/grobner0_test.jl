using AlgebraicSolvers, DynamicPolynomials, Groebner, LinearAlgebra

X = @polyvar x y

Mth = Grobner((P,X) -> Groebner.groebner(P, ordering = Groebner.DegRevLex(X)),
            (p,G) -> Groebner.normalform(p,G),
             G -> Groebner.quotient_basis(G)
            )

P = [-2+y-y^2+x^2*y, 1-3x+y+x*y^2]

B0 = quot_basis(Mth,P)
N, L = tnf(Mth,P)
M = mult_matrices(Mth,P)

Xi, G, B  = AlgebraicSolvers.solve(Mth, P; verbose=true)

Er = rel_error(P,Xi,X)
println("-- Rel error: ", norm(Er,Inf));

Xi
 
