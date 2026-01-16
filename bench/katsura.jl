using Groebner, DynamicPolynomials, AlgebraicSolvers, HomotopyContinuation

GB = Grobner((P,X) -> Groebner.groebner(P, ordering = Groebner.DegRevLex(X)),
             (p,G) -> Groebner.normalform(p,G),
              G -> Groebner.quotient_basis(G)
            )

for N in 4:6
    println("\nKatsura ",N)
    global P = Groebner.Examples.katsuran(N)
    t0 = @elapsed Xi, G, B  = AlgebraicSolvers.solve(GB,P,verbose=true)
    println("Nb roots: ",size(Xi,2))
    println("Groebner: ",t0, "(s)\n")

    X = (@polyvar x[1:N+1])[1]
    global P1 = [as_polynomial(p,X,Float64) for p in P]
    
    P2 = System(P1)
    t2 = @elapsed global Xi2 = HomotopyContinuation.solve(P2)
    println("Homotopy: ",t2, "(s)\n")


    t1 = @elapsed Xi1  = AlgebraicSolvers.solve(Macaulay(),P1,verbose=true)
    println("Nb roots: ",size(Xi1,2))
    println("Macaulay: ",t1, "(s)\n")

    #=
    t3 = @elapsed Xi3  = AlgebraicSolvers.solve(Toric(),P1,verbose=true)
    println("Nb roots: ",size(Xi2,2))
    println("Total: ",t2, "(s)\n")
    =#
end
