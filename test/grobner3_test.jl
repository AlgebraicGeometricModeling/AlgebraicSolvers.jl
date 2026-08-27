using Groebner, AlgebraicSolvers, LinearAlgebra

GB = GbSolver()

N=5
println("\nKatsura ",N)
P = Groebner.Examples.katsuran(N)
Xi, ms  = AlgebraicSolvers.solve(P,GB;verbose=true);

println("-- Mult sols: ", ms);
