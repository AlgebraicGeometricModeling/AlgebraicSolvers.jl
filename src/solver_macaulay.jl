export res_matrix,  tnf, quot_basis, solve, qr_basis, is_not_homogeneous

import LinearAlgebra, DynamicPolynomials

DP = DynamicPolynomials

export Macaulay

"""
Structure for the construction of Macaulay resultant solvers. It stores
  - `degree :: Function P->` degree ≤ ρ (default: `P ->  sum(DynamicPolynomials.maxdegree(P[i])-1 for i in 1:length(P)) + 1`) 
  -  `is_homogeneous :: Function P->` boolean testing if the system is homogeneous or not (default: `P -> !any(AlgebraicSolvers.is_not_homogeneous, P)`)

The default value for ρ is ∑ deg(pi) - n + 1.


This specifies
 - the resultant matrix construction associated to Macaulay construction `res_matrix(P,Macaulay())`
 - the TNF which is the transpose of the kernel of the resultant matrix `tnf(P,Macaulay())`
 - the solver which is use Macaulay TNF to solve the system  `solve(P,Macaulay())`

It can be used as an argument of [`solve`](@ref), [`tnf`](@ref), [`quot_basis`](@ref).

## Constructors
```
Macaulay()
```
The degree function is `P ->  sum(DynamicPolynomials.maxdegree(P[i])-1 for i in 1:length(P)) + 1`

```
Macaulay(rho::Int64)
```
The degree function is `P -> rho` 


## Example

```
using AlgebraicSolvers, DynamicPolynomials
X = @polyvar x y
P = [2-x*y-x^2, y^2+x-2, x^2-y^2+x-1]
Xi, ms = solve(P, Macaulay(3))
```
This gives the unique solution `Xi`:

```
2×1 Matrix{ComplexF64}:
 1.0 + 0.0im
 1.0 + 0.0im
```
"""
struct Macaulay
    degree :: Function 
    is_homogeneous :: Function 
end

function Macaulay()
    Macaulay(P ->  sum(DP.maxdegree(P[i])-1 for i in 1:length(P)) + 1,
             P -> !any(is_not_homogeneous, P))
end


function Macaulay(rho::Int64)
    Macaulay(P -> rho,
             P -> !any(is_not_homogeneous, P))
end

function is_not_homogeneous(p)
    L = [DP.maxdegree(t) for t in DP.monomials(p)]
    maximum(L) != minimum(L)
end


"""
```
R, L = res_matrix(P, Mth::Macaulay)
```
where `P` the polynomial system.

It outputs 
 - `R` the transpose of Sylvester matrix of all monomial multiples mi*pi in degree  ≤ ρ.
 - `L` array of monomials indexing the columns of `R`

The matrix `R` has a sparse representation.

"""
function res_matrix(P, Mth::Macaulay)
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

function res_matrix(P::AbstractVector, ::Val{:macaulay})
    res_matrix(P, Macaulay())
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






