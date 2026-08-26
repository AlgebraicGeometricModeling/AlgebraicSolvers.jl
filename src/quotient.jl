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

    N = nullspace(R)
    
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
function tnf(P::AbstractVector, Mth; verbose = false)

    t = @elapsed R, L = res_matrix(P, Mth)
    verbose && println("\033[96m-- Resultant matrix ", size(R,2),"x",size(R,1),  "   \033[0m",t, "(s)"); 

    rho = maxdegree(L)
    
    t = @elapsed N = LinearAlgebra.nullspace(R,tol=1.e-4)'
    verbose && println("\033[96m-- Cokernel ", size(N,1),"x",size(N,2), "   \033[0m $t(s)"); t0 = time()

    I0 = filter(i->maxdegree(L[i])<rho,1:length(L))
    N0 = N[:,I0]
    
    t = @elapsed F = qr(N0,Val(true))
    verbose && println("\033[96m-- QR fact.  \033[0m $t(s)");

    Jb= F.p[1:size(N0,1)]
    Ib = I0[Jb]
    
    #verbose && println("\033[96m-- Column basis\033[0m $(Ib)")
    
    return N, L, Ib
end

function tnf(P::AbstractVector, mth::Symbol)
      tnf(P, Val(mth))  
end


export mult_matrices
"""
```
    M = mult_matrices(P, DynamicPolynomials.variables(P), Mth)
```
Computes the vector of multiplication matrices `M=[M1, M2, ...]` by the variables in a basis `B` of the quotient by the ideal (`P`), using `tnf(P,Mth)`.

If the option `verbose=true` is set, it prints the different substeps and their timing.
"""
function mult_matrices(P::AbstractVector, X, Mth; verbose = false)
    N, L, IB = tnf(P, Mth; verbose = verbose)
    M = mult_matrices(N,L,IB,X)
    if length(M)==0
        @warn "Cannot find a good basis; Method not adapted "
        return Float64[], Int64[]
    else
        return M
    end
    
end

"""
```
    M = mult_matrices(N::Matrix, L::AbstractVector, IB::Vector, X)
```
Computes the vector of multiplication matrices `M=[M1, M2, ...]` by the variables `X` in the basis `B = L[IB]` from the Truncated Normal Form `N`, assuming `N[IB,IB]=Id`.
"""
function mult_matrices(N::AbstractMatrix, L::AbstractVector, IB::Vector, X)

    M0i = inv(N[:, IB])
    Idx = idx(L)

    M = Matrix{typeof(N[1,1])}[]
    
    for x in X
        J = [get(Idx,L[i]*x,0) for i in IB]
        if  findfirst(x-> x==0, J) == nothing
            push!(M, (M0i*N[:,J]))
        else
            @error "-- Basis*$x not in L"
        end
    end
    M
end
"""
```
    M = mult_matrices(N::Matrix, L::AbstractVector, IB::Vector, X, v0)
```
Computes the vector of multiplication matrices `M=[M1, M2, ...]` by the variables `X/v0` in a basis `B = L[IB]` from the Truncated Normal Form `N`, assuming `N[IB,IB]=Id`.
"""
function mult_matrices(N::Matrix, L::AbstractVector, IB::Vector, X, v0)

    M0i = inv(N[:, IB])
    
    Idx = idx(L)

    M = []
    x0 = monomial(v0)
    for x in X
        J = [get(Idx,div(L[i]*x,x0),0) for i in IB]
        if  findfirst(x-> x==0, J) == nothing
            push!(M, M0i*N[:,J])
        else
            @error "-- Basis*X not in L"
        end
    end
    M
    
end






