using AlgebraicSolvers, DynamicPolynomials

X = @polyvar x y

P = [2-x*y,x+y-2]

solve_toric(P,X)
