using LinearAlgebra, AbstractAlgebra, AlgebraicSolvers, Groebner

GB = Grobner((P,X) -> Groebner.groebner(P, ordering = Groebner.DegRevLex(X)),
              (p,G) -> Groebner.normalform(p,G),
              G -> Groebner.quotient_basis(G)
             )

R, (x,y) = QQ["x","y"]

n = 2
d = 3

M =  AbstractAlgebra.monomials((1+x+y)^d)

P = [x^2+1.0, y^2-2.0]

#P = [sum(m*rand(Int64) for m in M), sum(m*rand(Int64) for m in M) ]

Xi, ms, G,B = AlgebraicSolvers.solve(GB,P; verbose=true)

#println("-- sol ", Xi)

#Er = rel_error(P,Xi,[x,y])
#println("-- Rel error: ", norm(Er,Inf));

