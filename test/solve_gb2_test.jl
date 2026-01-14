using AlgebraicSolvers, DynamicPolynomials, Groebner

Mth = Grobner((P,X) -> Groebner.groebner(P, ordering = Groebner.DegRevLex(X)),
            (p,G) -> Groebner.normalform(p,G),
             G -> Groebner.quotient_basis(G)
            )

using AbstractAlgebra

R, (x,y,z) = polynomial_ring(AbstractAlgebra.QQ,[:x,:y, :z], internal_ordering = :degrevlex )

P = [-2+y-y^2+x^2*y+x*z, 1-3x+y+x*y^2, x*z^2+z^2+2*x*y-1]

Xi, G, B  = AlgebraicSolvers.solve(Mth, P;verbose=true)
 
