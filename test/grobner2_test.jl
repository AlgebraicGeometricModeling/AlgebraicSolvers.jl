using LinearAlgebra, AbstractAlgebra, AlgebraicSolvers, Groebner

GB = GbSolver(
          P->Groebner.groebner(P, ordering = Groebner.DegRevLex(variables(P))),
          Groebner.normalform,
          Groebner.quotient_basis
 )

R, (x,y) = QQ["x","y"]

n = 2
d = 3

M =  AbstractAlgebra.monomials((1+x+y)^d)

P = [x^2+1.0, y^2-2.0]

#P = [sum(m*rand(Int64) for m in M), sum(m*rand(Int64) for m in M) ]

Xi, ms = AlgebraicSolvers.solve(P,GB; verbose=true)

println("-- Mult sols: ", ms);

