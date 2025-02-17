using DynamicPolynomials, AlgebraicSolvers

X = @polyvar x y

P = [1-3y+2x*y,x*y-x+y-1]

R, L = matrix_toric(P) 

solve_toric(P; verbose=false)
