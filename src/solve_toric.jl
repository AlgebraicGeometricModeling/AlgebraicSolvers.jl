export support, solve, tnf, quo_basis, mult_matrices

import DynamicPolynomials

DP = DynamicPolynomials

export Toric

"""
Structure for the construction of Toric resultant solvers. It store
  - `supports :: Function  P->` supports of the polynomials `P` (default: `P->AlgebraicSolvers.support.(P)` )  
"""
struct Toric
   supports::Function
end

function support(p::DynamicPolynomials.Polynomial)
    sum(monomials(p))
end

function Toric()
    Toric(P->A = support.(P))
end

"""
```
R, L = res_matrix(P, Mth::Toric)
```
where
 - `P` polynomial system
 - `A` array of supports of `Pi`
 
It outputs 
 - `R` transpose of the Sylvester matrix of all monomial multiples mi*pi for mi in supp(∏_{j != i} pj).
 - `L` the list of monomials indexing the colums of `R`

"""
function res_matrix(P::AbstractVector, Mth::Toric)
    M = typeof(P[1])[]
    A = Mth.supports(P)
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
    R = sparse_matrix(M,idx(L))
    R, L
end








function res_matrix(P::AbstractVector, ::Val{:toric})  res_matrix(P, Toric()) end
function tnf(P::AbstractVector, ::Val{:toric})  tnf(P, Toric()) end
function quot_basis(P::AbstractVector, ::Val{:toric})  quot_basis(P, Toric()) end
function solve(P::AbstractVector, ::Val{:toric}; verbose::Bool = false )  solve(P, Toric(); verbose=verbose) end
