export matrix_macaulay,  tnf_macaulay, qr_basis, solve_macaulay, is_not_homogeneous

import LinearAlgebra, DynamicPolynomials

DP = DynamicPolynomials

function is_not_homogeneous(p)
    L = [DP.maxdegree(t) for t in DP.monomials(p)]
    maximum(L) != minimum(L)
end

"""
    R, L = matrix_macaulay(P, X, rho, ish = false)

 - `P` polynomial system
 - `X` (optional) array of variables
 - `rho` (optional) maximal degree of the multiples of P
 - `ish` (optional) set to true if the polynomials are homogeneous

It outputs 
 - `R` the transpose of Sylvester matrix of all monomial multiples mi*pi in degree ≤ rho.
 - `L` array of monomials indexing the columns of `R`

"""
function matrix_macaulay(P, X=DP.variables(P), rho =  sum(DP.maxdegree(P[i])-1 for i in 1:length(P)) + 1,  ish = false )
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
    matrix(M,idx(L)), L
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



"""
    N, L = tnf_macaulay(P, rho)

Compute the Truncated Normal Form of P=[p1, ..., pn], using Macaulay matrix of all monomial multiples mi*pi in degree ≤ ρ.

The default value for ρ is ∑ deg(pi) - n + 1.

"""
tnf_macaulay = function(P,  rho = sum(DP.maxdegree(P[i])-1 for i in 1:length(P)) + 1)
    
    ish = !any(is_not_homogeneous, P)
    X=DP.variables(P)
    R, L = matrix_macaulay(P, X, rho, ish)
    N = LinearAlgebra.nullspace(R)
    return N, L
end

"""
    Xi = solve_macaulay(P, rho ; verbose = false)

 - `P` polynomial system
 - `rho` (optional) degree of regularity for the Sylvester matrix construction (optional)

Solve the system P=[p1, ..., pn], building Sylvester matrix of all monomial multiples mi*pi in degree ≤ ρ.

The default value for ρ is ∑ deg(pi) - n + 1.

Example
-------
```
using AlgebraicSolvers, DynamicPolynomials

X = @polyvar x y

P = [2-x*y+x^2,y^2+x-2]

Xi = solve_macaulay(P)

```
"""
solve_macaulay = function(P, rho =  sum(DP.maxdegree(P[i])-1 for i in 1:length(P)) + 1;
                          verbose::Bool = false )
    
    verbose && println("-- Degrees ", map(p->DP.maxdegree(p),P))
    ish = !any(is_not_homogeneous, P)
    verbose && println("-- Homogeneity ", ish)
    t0 = time()
    #println("-- Monomials ", length(L), " degree ", rho,"   ",time()-t0, "(s)"); t0 = time()

    X = DP.variables(P)
    R, L = matrix_macaulay(P, X, rho, ish)
    
    verbose && println("-- Macaulay matrix ", size(R,1),"x",size(R,2),"   rho ",rho,"   ", time()-t0, "(s)"); t0 = time()

    N = LinearAlgebra.nullspace(R)
    verbose && println("-- Null space ",size(N,1),"x",size(N,2), "   ",time()-t0, "(s)"); t0 = time()

    B, Nr = qr_basis(N, L, ish)
    verbose && println("-- Qr basis ",  length(B), "   ",time()-t0, "(s)"); t0 = time()

    M = mult_matrix(B, X, Nr, idx(L), ish)
    verbose && println("-- Mult matrices ",time()-t0, "(s)"); t0 = time()

    Xi = eigdiag(M)
    verbose && println("-- Eigen diag",  "   ",time()-t0, "(s)"); t0 = time()
    if (!ish)
        for i in 1:size(Xi,2) Xi[:,i]/=Xi[1,i] end
        Xi = Xi[2:size(Xi,1),:]
    else
        for i in 1:size(Xi,2) Xi[:,i]/=norm(Xi[:,i]) end
    end
    Xi
    
end
