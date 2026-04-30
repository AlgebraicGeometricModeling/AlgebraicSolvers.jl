using DynamicPolynomials, AlgebraicSolvers

X = @polyvar x1 x2

#f = x2^2 - x1^2*x2

f = x2^2 - x1^3

F = [differentiate(f,x1), differentiate(f,x2)] 

D, B, h, Nu = invsys(F;verbose=true)
D
