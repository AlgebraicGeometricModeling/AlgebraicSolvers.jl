using DynamicPolynomials, AlgebraicSolvers

X = @polyvar x1 x2

F = [x1 - x2 + x1^2*1.0, x1 - x2 + x2^2*1.0]

D, B, h, Nu = invsys(F,[0,0];verbose=true)
D
