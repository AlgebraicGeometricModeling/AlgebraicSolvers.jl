using AlgebraicSolvers, DynamicPolynomials, Groebner, LinearAlgebra


GB=Grobner(
    Groebner.DegRevLex,      # function of the variables defining the monomial ordering
    Groebner.groebner,       # function for computing Grobner basis
    Groebner.normalform,     # function that computes the normal form of a polynomial
                             # w.r.t a Grobner basis 
    Groebner.quotient_basis  # function that computes a basis of the quotient algebra
)

X = @polyvar x y

P = [-2+y-y^2+x^2*y, 1-3x+y+x*y^2]
#P = [x^2+1, y^2-2]

X = variables(P)

G =GB.grobner_basis(P)

B0 = quot_basis(P,GB)

N, L = tnf(P,GB)
M = mult_matrices(P, variables(P), GB)

Xi, ms, G, B  = AlgebraicSolvers.solve(P, GB; verbose=true)

Er = rel_error(P,Xi)
println("-- Rel error: ", norm(Er,Inf));
println("-- Mult sols: ", ms);
 
