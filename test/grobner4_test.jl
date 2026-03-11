using AlgebraicSolvers, DynamicPolynomials, Groebner, LinearAlgebra


GB=Grobner(
    Groebner.DegRevLex,
    Groebner.groebner,
    Groebner.normalform,
    Groebner.quotient_basis
)

X = @polyvar x y z

P = [
    8*x^2*y^2+5*x*y^3+3*x^3*z+x^2*y*z,
    x^5+2*y^2*z^2+13*y^2*z^3+5*y*z^4,
    8*x^3+12*y^3+x*z^2+3,
    7*x^2*y^4+18*x*y^3*z^2+y^3*z^3
]

Xi, ms, G, B =  AlgebraicSolvers.solve(P, GB);
