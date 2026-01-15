export support, solve, tnf, quotient_basis, mult_matrices

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
R, L = res_matrix(Mth::Toric, P)
```
where
 - `P` polynomial system
 - `A` array of supports of `Pi`
 
It outputs 
 - `R` transpose of the Sylvester matrix of all monomial multiples mi*pi for mi in supp(∏_{j != i} pj).
 - `L` the list of monomials indexing the colums of `R`

"""
function res_matrix(Mth::Toric, P)
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
    R = matrix(M,idx(L))
    R, L
end
s
function res_matrix(::Val{:toric}, P)  res_matrix(Toric(),P) end

function tnf(::Val{:toric}, P)  tnf(Toric(),P) end

function solve(::Val{:toric}, P; verbose::Bool = false )  solve(Toric(),P; verbose=verbose) end
