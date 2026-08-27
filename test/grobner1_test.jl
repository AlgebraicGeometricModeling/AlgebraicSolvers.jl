using AlgebraicSolvers, DynamicPolynomials, Groebner, LinearAlgebra

GB = GbSolver(
          P->Groebner.groebner(P, ordering = Groebner.DegRevLex(variables(P))),
          Groebner.normalform,
          Groebner.quotient_basis
 )

X = @polyvar x y

P = [-2+y-y^2+x^2*y, 1-3x+y+x*y^2]
#P = [x^2+1, y^2-2]

X = variables(P)

G =GB.grobner_basis(P)

B0 = quot_basis(P,GB)

N, L = tnf(P,GB)
M = mult_matrices(P, variables(P), GB)

Xi, ms  = AlgebraicSolvers.solve(P, GB; verbose=true)
println("-- Mult sols: ", ms);
 
