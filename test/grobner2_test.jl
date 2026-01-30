    using Groebner, AlgebraicSolvers, DynamicPolynomials
GB = Grobner((P,X) -> Groebner.groebner(P, ordering = Groebner.DegRevLex(X)),
             (p,G) -> Groebner.normalform(p,G),
              G -> Groebner.quotient_basis(G)
             )

N=4
println("\nKatsura ",N)
P = Groebner.Examples.katsuran(N)
Xi, G, B  = AlgebraicSolvers.solve(GB,P,verbose=true);


P1 = convert_DP(P)
X = DynamicPolynomials.variables(P1)
Er = rel_error(P1, Xi, X)
println("-- Rel error: ", norm(Er,Inf));
