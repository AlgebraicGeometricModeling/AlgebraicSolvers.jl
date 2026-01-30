export solve, res_matrix, tnf, quot_basis, mult_matrices

function res_matrix(mth::Symbol, P::Vector)
      res_matrix(Val(mth),P)  
end


"""
```
    B = quot_basis(Mth, P)
```
Computes a basis of the quotient by the ideal P, using `res_matrix(Mth,P)`

"""
function quot_basis(Mth, P)

    R, L = res_matrix(Mth,P)

    N, IB = nullspace(R)
    
    F = qr(N')
    
    return L[column_basis(F.R)]
end

function quot_basis(mth::Symbol, P::Vector)
      quot_basis(Val(mth),P)  
end


"""
    N, L, B = tnf(Mth, P)

Compute the Truncated Normal Form `N` of `P=[p1, ..., pn]`, using `res_matrix(Mth,P)`.

The list `L` is the list of monomials indexing the colmuns of `N`.

"""
function tnf(Mth, P)
    R, L = res_matrix(Mth, P)
    N, IB = LinearAlgebra.nullspace(R)
    return N', L, IB
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
function mult_matrices(Mth, P, X = DynamicPolynomials.variables(P))
    R, L = res_matrix(Mth,P)
    N, _ = LinearAlgebra.nullspace(R)
    F  = qr(N')
    IB = column_basis(F.R)
    _mult_matrices(F.R,L,IB,X)

end

export _mult_matrices
function _mult_matrices(N, L, IB, X)
    M0i = inv(N[:, IB])
    
    Idx = idx(L)

    M = Matrix{typeof(N[1,1])}[]
    
    for x in X
        J = [get(Idx,L[i]*x,0) for i in IB]
        if  findfirst(x-> x==0, J) == nothing
            push!(M, (M0i*N[:,J])')
        else
            @error "-- Basis*X not in L"
        end
    end
    M
end

function _mult_matrices_h(N, L, IB, X, v0)

    #println("---> ",v0)
    M0i = inv(N[:, IB])
    
    Idx = idx(L)

    M = []
    x0 = monomial(v0)
    for x in X
        J = [get(Idx,div(L[i]*x,x0),0) for i in IB]
            #println("---> ",[div(L[i]*x,x0) for i in IB])
        if  findfirst(x-> x==0, J) == nothing
            push!(M, M0i*N[:,J])
        else
            @error "-- Basis*X not in L"
        end
    end
    M
    
end


"""

    Xi, ms = solve(Mth, P; verbose = false)

 - `P` polynomial system

Solve the system or polynomials `P=[p1, ..., pn]`, building Sylvester matrix `res_matrix(Mth,P)`.
It outputs the solutions `Xi` as a matrix of points, one per column, and the vector of their multiplicities `ms`.

Example
=======
```
using AlgebraicSolvers, DynamicPolynomials

X = @polyvar x y

P = [y - x*y + x^2+ y^2,  1 + y + x + x^2]

Xi, ms = solve(Macaulay(), P; verbose=true)
```
"""
function solve(Mth, P; verbose::Bool=false)

    X = DynamicPolynomials.variables(P)
    
    t0 = time()
    R, L = res_matrix(Mth,P)
    verbose && println("\033[96m-- Resultant matrix ", size(R,1),"x",size(R,2),  "   \033[0m",time()-t0, "(s)"); t0 = time()

    N, IB = LinearAlgebra.nullspace(R)
    verbose && println("\033[96m-- Null space ", size(N,1),"x",size(N,2), "   \033[0m",time()-t0, "(s)"); t0 = time()

    Nt = N'
    F = qr!(Nt)
    IB = column_basis(F.R)

    verbose && println("\033[96m-- Basis ", length(IB), "  \033[0m",time()-t0, "(s)"); t0 = time()

    M = _mult_matrices(F.R, L, IB, X)
    verbose && println("\033[96m-- Mult matrices \033[0m",time()-t0, "(s)"); t0 = time()

    #Xi = eigdiag(M)
    Xi, ms = schur_dcp(M)
    verbose && println("\033[96m-- Eigen diag",  "   \033[0m",time()-t0, "(s)"); t0 = time()

    Xi, ms
end

function solve(mth::Symbol, P::Vector; verbose::Bool = false)
      solve(Val(mth),P;verbose=verbose)  
end


