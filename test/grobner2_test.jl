using Groebner, AlgebraicSolvers
GB = Grobner((P,X) -> Groebner.groebner(P, ordering = Groebner.DegRevLex(X)),
             (p,G) -> Groebner.normalform(p,G),
              G -> Groebner.quotient_basis(G)
             )

N=5
println("\nKatsura ",N)
P = Groebner.Examples.katsuran(N)
Xi, G, B  = AlgebraicSolvers.solve(GB,P,verbose=true)
Xi
