using AlgebraicSolvers, DynamicPolynomials

X = @polyvar x y


P = [2-x*y-x^2, x*y-x+y-1]

Xi = solve_macaulay(P,X)
