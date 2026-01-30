using Groebner, AlgebraicSolvers
GB = Grobner((P,X) -> Groebner.groebner(P, ordering = Groebner.DegRevLex(X)),
             (p,G) -> Groebner.normalform(p,G),
              G -> Groebner.quotient_basis(G)
             )

N=4
println("\nKatsura ",N)
P = Groebner.Examples.katsuran(N)
Xi, ms, G, B  = AlgebraicSolvers.solve(GB,P,verbose=true);

Er = rel_error(P, Xi)
println("-- Rel error: ", norm(Er,Inf));
println("-- Mult sols: ", ms);
