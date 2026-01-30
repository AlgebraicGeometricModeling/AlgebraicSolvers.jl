using DynamicPolynomials, AlgebraicSolvers, LinearAlgebra, SparseArrays

X = @polyvar x0 x1 x2 x3 x4 

P=[ x0^2 + 2*x1^2 + 2*x2^2 + 2*x3^2 + 2*x4^2 - x0,
 2*x0*x1 + 2*x1*x2 + 2*x2*x3 + 2*x3*x4 - x1,
 x1^2 + 2*x0*x2 + 2*x1*x3 + 2*x2*x4 - x2,
 2*x1*x2 + 2*x0*x3 + 2*x1*x4 - x3,
 x0 + 2*x1 + 2*x2 + 2*x3 + 2*x4 - 1
]

DP = DynamicPolynomials

Mth = Macaulay()
verbose = true

rho = Mth.degree(P)
ish = Mth.is_homogeneous(P)

verbose && println("\033[96m-- Degrees = ", map(p->DP.maxdegree(p),P),"   rho = ", rho, "   Homogeneity = ", ish, "\033[0m")

X = DP.variables(P)

t = @elapsed R, L = res_matrix(Mth, P)
verbose && println("\033[96m-- Macaulay matrix ", size(R,1),"x",size(R,2),"  \033[0m", t, "(s)"); t0 = time()

N, _ = LinearAlgebra.nullspace(R)
verbose && println("\033[96m-- Null space ",size(N,1),"x",size(N,2), "   \033[0m",time()-t0, "(s)"); t0 = time()

Nt = N';
F = qr(Nt)
IB = column_basis(F.R)

verbose && println("\033[96m-- Basis ", length(IB), "   \033[0m",time()-t0, "(s)"); t0 = time(); 
    
M = _mult_matrices(F.R, L, IB, X)

verbose && println("\033[96m-- Mult matrices  \033[0m",time()-t0, "(s)"); t0 = time()

Xi = eigdiag(M)
#Xi, E, Info = MultivariateSeries.diagonalization(M)

verbose && println("\033[96m-- Eigen diag",  "  \033[0m ",time()-t0, "(s)"); t0 = time()

#=
R, L = res_matrix(Mc,P)

N, L = tnf(Mc,P)

B = quot_basis(:macaulay,P)

N, _ = LinearAlgebra.nullspace(R)

F = qr!(N')
U = F.R 
IB = column_basis(U)

M = mult_matrices(Mc,P)

Xi = eigdiag(M)

Xi = AlgebraicSolvers.solve(Mc, P; verbose=true)
=#

Er = rel_error(P, Xi, X)
println("-- Rel error: ", norm(Er,Inf));

#println("-- Sol: ", Xi)
