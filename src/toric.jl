export support, matrix_toric, solve_toric, tnf_toric

import DynamicPolynomials

DP = DynamicPolynomials

function support(p::DynamicPolynomials.Polynomial)
    sum(monomials(p))
end

"""
```
R, L = matrix_toric(P, A)
```
where
 - `P` polynomial system
 - `A` array of supports of `Pi`
 
It outputs 
 - `R` transpose of the Sylvester matrix of all monomial multiples mi*pi for mi in supp(∏_{j != i} pj).
 - `L` the list of monomials indexing the colums of `R`

"""
function matrix_toric(P, A = support.(P))
    M = typeof(P[1])[]
    mult = one(A[1])
    for i in 1:length(P)
        mult = one(A[1])
        for j in 1:length(A)
            if j!= i
                mult*=A[j]
            end
        end
        for m in reverse(DynamicPolynomials.monomials(mult))
            push!(M,P[i]*m)
        end
    end
    mult *= A[end]
    L = (DynamicPolynomials.monomials(mult))
    R = matrix(M,idx(L))
    R, L
end

"""
    N, L = tnf_toric(P, A = support.(P))

Compute the Truncated Normal Form of P=[p1, ..., pn], using toric resultant matrix of all monomial multiples mi*pi in degree ≤ ρ.

The default value for ρ is ∑ deg(pi) - n + 1.

"""
tnf_toric = function(P,  A = support.(P))

    R, L = matrix_toric(P, A)
    N = LinearAlgebra.nullspace(R)
    return N, L
end


"""

    solve_toric(P; verbose = false)

 - `P` polynomial system
 - `X` array of variables

Solve the system `P=[p1, ..., pn]`, building Sylvester matrix of all monomial multiples mi*pi for mi in supp(∏_{j != i} pj).


Example
=======
```
using AlgebraicSolvers, DynamicPolynomials

X = @polyvar x y

P = [y - x*y + x^2,  1 + y + x + x^2]

Xi = solve_toric(P)

```
"""
function solve_toric(P; verbose::Bool=false)
    t0 = time()
    #A = [support(p) for p in P]
    R, L = matrix_toric(P)
    verbose && println("-- Toric matrix ", size(R,1),"x",size(R,2),  "   ",time()-t0, "(s)"); t0 = time()
    #println("-- L ", L)
    N = LinearAlgebra.nullspace(R)
    verbose && println("-- Null space ",size(N,1),"x",size(N,2), "   ",time()-t0, "(s)"); t0 = time()

    X = DynamicPolynomials.variables(P);
    B = tnf_basis(N, L, X)
    if B == nothing
        return nothing
    end
    verbose && println("-- Basis ", B, "  ", time()-t0, "(s)"); t0 = time()

    
    M = mult_matrix(B, X, N, idx(L), false)
    verbose && println("-- Mult matrices ",time()-t0, "(s)"); t0 = time()

    Xi = eigdiag(M)
    verbose && println("-- Eigen diag",  "   ",time()-t0, "(s)"); t0 = time()

    for i in 1:size(Xi,2) Xi[:,i]/=Xi[1,i] end
    Xi = Xi[2:size(Xi,1),:]
    Xi
end
