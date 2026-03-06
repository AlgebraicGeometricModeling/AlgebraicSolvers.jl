using DynamicPolynomials, LinearAlgebra, AlgebraicSolvers

X = @polyvar X0 X1 X2 X3 X4 X5

F = 24*X0^3 + 70*X0^2*X1 + 75*X0^2*X2 + 70*X0^2*X3 + 180*X0^2*X4 + 10*X0^2*X5 + 10*X0*X1^2+ 70*X0*X2^2 + 360*X0*X2*X3 + 120*X0*X2*X4 + 60*X0*X3^2 + 60*X2^3 + 60*X2^2*X3

G = 6*X0^4 +70//3*X0^3*X1 + 25*X0^3*X2 +70//3*X0^3*X3 + 60*X0^3*X4 + 10//3*X0^3*X5 + 5*X0^2*X1^2 + 35*X0^2*X2^2+ 180*X0^2*X2*X3 + 60*X0^2*X2*X4 + 30*X0^2*X3^2 + 60*X0*X2^3 + 60*X0*X2^2*X3 + 5*X2^4

DP = DynamicPolynomials

L = DP.monomials(F)
m = L[1]
dF = dual(F,3)

L1 = reverse(DP.monomials(X,1))
L2 = reverse(DP.monomials(X,2))

H = hankel(dF,L1,L2)

s = series(H, L1, L2)

s1 = dual(subs(dual(dF),X0=>1))
D1 = invsys(s1)
M1, L1 = matrixof(D1)

@assert rank(Matrix{Float64}(M1)) == 7

dG = dual(G,4)
s2 = dual(subs(dual(dG),X0=>1))
D2 = invsys(s2)
M2, L2 = matrixof(D2)

@assert rank(Matrix{Float64}(M2)) == 6


4*X0*dG-dF


