using LinearAlgebra, AbstractAlgebra, AlgebraicSolvers, Groebner

GB =Grobner(
    Groebner.DegRevLex,
    Groebner.groebner,
    Groebner.normalform,
    Groebner.quotient_basis
)
R, (x,y) = QQ["x","y"]

n = 2
d = 3

M =  AbstractAlgebra.monomials((1+x+y)^d)

P = [x^2+1.0, y^2-2.0]

#P = [sum(m*rand(Int64) for m in M), sum(m*rand(Int64) for m in M) ]

Xi, ms, G,B = AlgebraicSolvers.solve(P,GB; verbose=true)

#println("-- sol ", Xi)

Er = rel_error(P,Xi)
println("-- Rel error: ", norm(Er,Inf));
println("-- Mult sols: ", ms);

