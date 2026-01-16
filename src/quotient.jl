export solve, res_matrix, tnf, quo_basis, mult_matrices

function res_matrix(mth::Symbol, P::Vector)
      res_matrix(Val(mth),P)  
end


"""
```
    B = quo_basis(Mth, P)
```
Computes a basis of the quotient by the ideal P, using `res_matrix(Mth,P)`

"""
function quo_basis(Mth, P)
    R, L = res_matrix(Mth,P)
    N = LinearAlgebra.nullspace(R)
    X = DynamicPolynomials.variables(P);
    B = tnf_basis(N, L, X)
    return B
end

function quo_basis(mth::Symbol, P::Vector)
      quo_basis(Val(mth),P)  
end



"""
    N, L = tnf(Mth, P)

Compute the Truncated Normal Form `N` of `P=[p1, ..., pn]`, using `res_matrix(Mth,P)`.

The list `L` is the list of monomials indexing the colmuns of `N`.

"""
function tnf(Mth, P)
    R, L = res_matrix(Mth, P)
    N = LinearAlgebra.nullspace(R)
    return N', L
end

function tnf(mth::Symbol, P::Vector)
      tnf(Val(mth),P)  
end


export mult_matrices

"""
```
    M = mult_matrices(Mth, P, A = DynamicPolynomials.variables(P))
```
Computes the vector of multiplication matrices `M=[1, M2, ...]` by the variables in a basis `B` of the quotient by the ideal (`P`), using `res_matrix(Mth,P)`.
"""
function mult_matrices(Mth, P, A = DynamicPolynomials.variables(P))
    R, L = res_matrix(Mth,P)
    N = LinearAlgebra.nullspace(R)

    X = DynamicPolynomials.variables(P);
    B = tnf_basis(N, L, X)
    if B == nothing
        return nothing
    end
    M = mult_matrix(B, A, N, idx(L), false)
    M0i=inv(M[1])
    [M0i*m for m in M[2:end]]
end


"""

    solve(Mth, P; verbose = false)

 - `P` polynomial system

Solve the system or polynomials `P=[p1, ..., pn]`, building Sylvester matrix `res_matrix(Mth,P)`.


Example
=======
```
using AlgebraicSolvers, DynamicPolynomials

X = @polyvar x y

P = [y - x*y + x^2+ y^2,  1 + y + x + x^2]

Xi = solve(Macaulay(), P; verbose=true)
```
"""
function solve(Mth, P; verbose::Bool=false)
    t0 = time()
    #A = [support(p) for p in P]
    R, L = res_matrix(Mth,P)
    verbose && println("\033[96m-- Resultant matrix ", size(R,1),"x",size(R,2),  "   \033[0m",time()-t0, "(s)"); t0 = time()
    #println("-- L ", L)
    N = LinearAlgebra.nullspace(R)
    verbose && println("\033[96m-- Null space ",size(N,1),"x",size(N,2), "   \033[0m",time()-t0, "(s)"); t0 = time()

    X = DynamicPolynomials.variables(P);
    B = tnf_basis(N, L, X)
    if B == nothing
        return nothing
    end
    verbose && println("\033[96m-- Basis ", B, "  \033[0m", time()-t0, "(s)"); t0 = time()

    
    M = mult_matrix(B, X, N, idx(L), false)
    verbose && println("\033[96m-- Mult matrices \033[0m",time()-t0, "(s)"); t0 = time()

    Xi = eigdiag(M)
    verbose && println("\033[96m-- Eigen diag",  "   \033[0m",time()-t0, "(s)"); t0 = time()

    for i in 1:size(Xi,2) Xi[:,i]/=Xi[1,i] end
    Xi = Xi[2:size(Xi,1),:]
    Xi
end

function solve(mth::Symbol, P::Vector; verbose::Bool = false)
      solve(Val(mth),P;verbose=true)  
end


