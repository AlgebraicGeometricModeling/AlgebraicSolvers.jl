export solve, res_matrix, tnf, quot_basis, mult_matrices

function res_matrix( P::AbstractVector, mth::Symbol)
      res_matrix(P,Val(mth))  
end


"""
```
    B = quot_basis(P, Mth)
```
Computes a basis of the quotient by the ideal P, using `res_matrix(P,Mth)`

"""
function quot_basis(P, Mth)

    R, L = res_matrix(P, Mth)

    N, IB = nullspace(R)
    
    F = qr(N')
    
    return L[column_basis(F.R)]
end

function quot_basis(P::AbstractVector, mth::Symbol)
      quot_basis(P, Val(mth))  
end


"""
    N, L, B = tnf(P, Mth)

Compute the Truncated Normal Form `N` of `P=[p1, ..., pn]`, using `res_matrix(P,Mth)`.

The list `L` is the list of monomials indexing the colmuns of `N`.

"""
function tnf(P::AbstractVector, Mth)
    R, L = res_matrix(P, Mth)
    N, IB = LinearAlgebra.nullspace(R)
    return N', L, IB
end

function tnf(P::AbstractVector, mth::Symbol)
      tnf(P, Val(mth))  
end


export mult_matrices
"""
```
    M = mult_matrices(P, DynamicPolynomials.variables(P), Mth)
```
Computes the vector of multiplication matrices `M=[M1, M2, ...]` by the variables in a basis `B` of the quotient by the ideal (`P`), using `res_matrix(P,Mth)`.
"""
function mult_matrices(P::AbstractVector, X, Mth)
    R, L = res_matrix(P, Mth)
    N, _ = LinearAlgebra.nullspace(R)
    F  = qr(N')
    IB = column_basis(F.R)
    mult_matrices(F.R,L,IB,X)
end

"""
Computes the vector of multiplication matrices `M=[M1, M2, ...]` by the variables `X` in a basis `B = L[IB]` from the Truncated Normal Form `N`, assuming `N[IB,IB]=Id`.
"""
function mult_matrices(N::Matrix, L::AbstractVector, IB::Vector, X)

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
"""
Computes the vector of multiplication matrices `M=[M1, M2, ...]` by the variables `X` in a basis `B = L[IB]` from the Truncated Normal Form `N`, assuming `N[IB,IB]=Id`.
"""
function mult_matrices(N::Matrix, L::AbstractVector, IB::Vector, X, v0)

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

    Xi, ms = solve(P, Mth; verbose = false)

 - `P` polynomial system
 - `Mth` class specifying the solver.

Solve the system of polynomials `P=[p1, ..., pn]`, using Sylvester matrix `res_matrix(P,Mth)`.
It outputs the solutions `Xi` as a matrix of points, one per column, and the vector of their multiplicities `ms`.
"""
function solve(P::AbstractVector, Mth; verbose::Bool=false)

    X = DynamicPolynomials.variables(P)
    
    t0 = time()
    R, L = res_matrix(P, Mth)
    verbose && println("\033[96m-- Resultant matrix ", size(R,1),"x",size(R,2),  "   \033[0m",time()-t0, "(s)"); t0 = time()

    N, IB = LinearAlgebra.nullspace(R)
    verbose && println("\033[96m-- Null space ", size(N,1),"x",size(N,2), "   \033[0m",time()-t0, "(s)"); t0 = time()

    Nt = N'
    F = qr!(Nt)
    IB = column_basis(F.R)

    verbose && println("\033[96m-- Basis ", length(IB), "  \033[0m",time()-t0, "(s)"); t0 = time()

    M = mult_matrices(F.R, L, IB, X)
    verbose && println("\033[96m-- Mult matrices \033[0m",time()-t0, "(s)"); t0 = time()

    #Xi = eigdiag(M)
    Xi, ms = schur_dcp(M)
    verbose && println("\033[96m-- Eigen diag",  "   \033[0m",time()-t0, "(s)"); t0 = time()

    Xi, ms
end

function solve(P::AbstractVector, mth::Symbol; verbose::Bool = false)
      solve(P,Val(mth);verbose=verbose)  
end


