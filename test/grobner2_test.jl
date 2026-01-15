using Groebner, AlgebraicSolvers
GB = Grobner((P,X) -> Groebner.groebner(P, ordering = Groebner.DegRevLex(X)),
              (p,G) -> Groebner.normalform(p,G),
              G -> Groebner.quotient_basis(G)
             )

for N in 4:8
    println("\nKatsura ",N)
    global P = Groebner.Examples.katsuran(N)
    global Xi, G, B  = AlgebraicSolvers.solve(GB,P,verbose=true)
    println("Nb roots: ",size(Xi,2))
end
