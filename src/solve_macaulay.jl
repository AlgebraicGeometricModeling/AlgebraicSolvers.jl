export res_matrix,  tnf, quot_basis, solve, qr_basis, is_not_homogeneous

import LinearAlgebra, DynamicPolynomials

DP = DynamicPolynomials

export Macaulay

"""
Structure for the construction of Macaulay resultant solvers. It stores
  - `degree :: Function P->` degree of regularity +1 (default: `P ->  sum(DynamicPolynomials.maxdegree(P[i])-1 for i in 1:length(P)) + 1`) 
  -  `is_homogeneous :: Function P->` boolean testing if the system is homogeneous or not (default: `P -> !any(AlgebraicSolvers.is_not_homogeneous, P)`)

"""
struct Macaulay
    degree :: Function 
    is_homogeneous :: Function 
end


function Macaulay()
    Macaulay(P ->  sum(DP.maxdegree(P[i])-1 for i in 1:length(P)) + 1,
             P -> !any(is_not_homogeneous, P))
end


function is_not_homogeneous(p)
    L = [DP.maxdegree(t) for t in DP.monomials(p)]
    maximum(L) != minimum(L)
end


"""
    R, L = res_matrix(Mth::Macaulay, P, X, rho, ish = false)

 - `P` polynomial system
 - `X` (optional) array of variables
 - `rho` (optional) maximal degree of the multiples of P
 - `ish` (optional) set to true if the polynomials are homogeneous

It outputs 
 - `R` the transpose of Sylvester matrix of all monomial multiples mi*pi in degree ≤ rho.
 - `L` array of monomials indexing the columns of `R`

"""
function res_matrix(Mth::Macaulay, P, Q = nothing)
    X=DP.variables(P)
    rho = Mth.degree(P)
    ish = Mth.is_homogeneous(P)
    
    if ish
        L = [m for m in DP.monomials(X, rho)]
        Q = [DP.monomials(X,rho-DP.maxdegree(P[i])) for i in 1:length(P)]
    else
        L = [m for m in DP.monomials(X, 0:rho)]
        Q = [DP.monomials(X,0:rho-DP.maxdegree(P[i])) for i in 1:length(P)]
    end

    M = []
    for i in 1:length(P)
        for m in Q[i]
            push!(M,P[i]*m)
        end
    end
    sparse_matrix(M,idx(L)), L
end

function qr_basis(N, L, ish = false)
    Idx= idx(L)
    if ish
        L0 = filter(m->(m.z[1]>0), L)
    else
        d  = maximum([DP.maxdegree(m) for m in L])
        L0 = filter(m->DP.maxdegree(m)<d,L)
    end
    N0 = fill(zero(N[1,1]), size(N,2),length(L0))
    for i in 1:length(L0)
        for j in 1:size(N,2)
            N0[j,i]= N[get(Idx,L0[i],0),j]
        end
    end

    # F= qrfact(N0, Val(true))
    F = LinearAlgebra.qr(N0,Val(true))
    B = []
    if ish
        for i in 1:size(N,2)
            m = copy(L0[F.p[i]])
            m.z[1]-=1
            push!(B, m)
            # should test if the diag. coeff. is not small
        end
    else
        for i in 1:size(N,2)
            push!(B, L0[F.p[i]])
        # should test if the diag. coeff. is not small
        end
    end
    B, N*F.Q
end

function res_matrix(::Val{:macaulay}, P) res_matrix(Macaulay(),P) end
function tnf(::Val{:macaulay}, P) tnf(Macaulay(),P) end
function quot_basis(::Val{:macaulay}, P) quot_basis(Macaulay(),P) end
function solve(::Val{:macaulay}, P; verbose::Bool = false ) solve(Macaulay(),P; verbose=verbose) end


"""
    Xi, ms = solve(Macaulay(), P)

 - `P` polynomial system


Solve the system P=[p1, ..., pn], building Sylvester matrix of all monomial multiples mi*pi in degree ≤ ρ.
It outputs the solutions `Xi` as a matrix of points, one per column, and the vector of their multiplicities `ms`.

The default value for ρ is ∑ deg(pi) - n + 1.

Example
-------
```
using AlgebraicSolvers, DynamicPolynomials

X = @polyvar x y

P = [2-x*y+x^2,y^2+x-2]

Xi = solve(Macaulay(), P)

```
"""
function solve(Mth::Macaulay, P;
               verbose::Bool = false )

    rho = Mth.degree(P)
    ish = Mth.is_homogeneous(P)

    verbose && println("\033[96m-- Degrees = ", map(p->DP.maxdegree(p),P),"   rho = ", rho, "   Homogeneity = ", ish, "\033[0m")

    X = DP.variables(P)

    t = @elapsed R, L = res_matrix(Mth, P)
    verbose && println("\033[96m-- Macaulay matrix ", size(R,1),"x",size(R,2),"  \033[0m", t, "(s)"); t0 = time()
    if ish 
        L = reverse(L) # graded reverse lex order if default order is Graded{LexOrder}
        R = R[:,end:-1:1] 
    end
    
    N, _ = LinearAlgebra.nullspace(R)
    verbose && println("\033[96m-- Null space ",size(N,1),"x",size(N,2), "   \033[0m",time()-t0, "(s)"); t0 = time()

    Nt = N';
    F = qr!(Nt)
    IB = column_basis(F.R)

    verbose && println("\033[96m-- Basis ", length(IB), "   \033[0m",time()-t0, "(s)"); t0 = time(); 
    
    if ish
        M = _mult_matrices_h(F.R, L, IB, X, X[1])
    else
        M = _mult_matrices(F.R, L, IB, X)
    end
    
    verbose && println("\033[96m-- Mult matrices  \033[0m",time()-t0, "(s)"); t0 = time()

    Xi, ms = schur_dcp(M)
    verbose && println("\033[96m-- Eigen diag",  "  \033[0m ",time()-t0, "(s)"); t0 = time()
    if (ish)
        for i in 1:size(Xi,2) Xi[:,i]/=LinearAlgebra.norm(Xi[:,i]) end
    end
    Xi, ms
    
end





