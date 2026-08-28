export support, solve, tnf, quo_basis, mult_matrices

import DynamicPolynomials

DP = DynamicPolynomials

export Toric

"""
Structure for the construction of Toric resultant solvers. It store
  - `supports :: Function  P->` supports of the polynomials `P` (default: `P->AlgebraicSolvers.support.(P)` )

This specifies
 - the resultant matrix  associated to the Toric construction `res_matrix(P, Toric())`
 - the TNF which is the transpose of the kernel of the resultant matrix `tnf(P,Toric())`
 - the solver which is using the Toric TNF to solve the system  `solve(P,Toric())`

"""
struct Toric
   supports::Function
end

function support(p::DynamicPolynomials.Polynomial)
    monomials(p)
end

function Toric()
    Toric(P->A = support.(P))
end

function support_diff(L, x)
    monomials(differentiate(sum(L),x))
end
export support_diff

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
    A = sum.(Mth.supports(P))
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

function filter_basis(N::AbstractMatrix, L::AbstractVector, Mth::Toric)
    X = variables(L)
    M = intersect(L, [support_diff(L,x) for x in X]...)
    return [findfirst(m -> m==M[i],L) for i in 1:length(M)]
    
    #column_basis(N)
end

