using Groebner, AlgebraicSolvers, LinearAlgebra

GB=Grobner(
    Groebner.DegRevLex,
    Groebner.groebner,
    Groebner.normalform,
    Groebner.quotient_basis
)

N=5
println("\nKatsura ",N)
P = Groebner.Examples.katsuran(N)
Xi, ms, G, B  = AlgebraicSolvers.solve(P,GB;verbose=true);

Er = rel_error(P, Xi)
println("-- Rel error: ", norm(Er,Inf));
println("-- Mult sols: ", ms);
