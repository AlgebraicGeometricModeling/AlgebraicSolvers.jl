using LinearAlgebra, AbstractAlgebra, AlgebraicSolvers, Groebner

GB = Grobner((P,X) -> Groebner.groebner(P, ordering = Groebner.DegRevLex(X)),
              (p,G) -> Groebner.normalform(p,G),
              G -> Groebner.quotient_basis(G)
             )

R, (x,y) = QQ["x","y"]

n = 2
d = 3

M =  AbstractAlgebra.monomials((1+x+y)^d)

#P = [x1^2+1.0, x2^2-2.0]
P = [
    sum(m*rand(Int64) for m in M),
    sum(m*rand(Int64) for m in M)
]

Xi,G,B = solve(GB,P)
#println("-- sol ", Xi)

#Er = rel_error(P,Xi,X)
#println("-- Rel error: ", norm(Er,Inf));

