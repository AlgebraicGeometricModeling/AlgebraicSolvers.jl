export matrix, sparse_matrix, tnf_basis, mult_matrix, eigdiag, kernel, rel_error, idx

import DynamicPolynomials: coefficients, monomials

using LinearAlgebra, SparseArrays

function idx(L)
    return Dict{typeof(L[1]), Int64}(L[i]=> i for i in 1:length(L))
end

function kernel(A::Matrix)
    U,S,V = svd(A)
    r=1;
    while r<=min(size(A,1),size(A,2)) && abs(S[r])> 1.e-4
        r+=1
    end
    r-=1
    #P = F[:p]
    #Pi = fill(0, length(P))
    #for i in 1:length(P) Pi[P[i]]= i end
    V[:,r+1:end]
end

export nullspace
function LinearAlgebra.nullspace(R::AbstractSparseMatrix; tol::Float64=1.e-10)

    s = size(R,2)

    Rr = R[:,end:-1:1]

    F = qr(Rr)
    l = rank(F.R; tol = tol)
    r = s-l
    

    IB = [i for i in s:-1:1][F.pcol]
    IB = IB[l+1:end]

    U = F.R
    G = Array(U[1:l,l+1:end])
    N = U[1:l,1:l]\G

    N = vcat(-N, Matrix(I,r,r))
    N = N[reverse(invperm(F.pcol)),:]

#=
    F = qr(R)
    
    l = rank(F.R)

    r = s-l

    R0 = (R[end:-1:1,F.pcol])[:,1:l]
    
    F0 = lu(R0)

    #IB = [i for i in s:-1:1][F0.q]
    IB = [i for i in s:-1:1][F0.p]
    IB = IB[l+1:end]

    #U = F0.U
    U = spzeros(size(R0,2), size(R0,1))
    transpose!(U,F0.L)

    G  = Array(U[1:l,l+1:end])
    N = U[1:l,1:l]\G
    
    N = vcat(-N, Matrix(I,r,r))

    N = (F0.Rs[F0.p]).*N
   
    #    N = N[:,reverse(invperm(F0.p))]
    N = N[reverse(invperm(F0.p)),:]
=#
    return N, IB
end

#=
function LinearAlgebra.nullspace(R::AbstractSparseMatrix, L::Vector)
    Rt= R'
    s = size(R,2)
    F = qr(Rt)
    l = rank(F.R)
    r = s-l

    R0 = (Rt[end:-1:1,F.pcol])[:,1:l]
    
    F0 = lu(R0)

    U = spzeros(size(R0,2), size(R0,1))
    transpose!(U,F0.L)
    
    G  = Array(U[1:l,l+1:end])
    
    N = U[1:l,1:l]\G
    
    N = vcat(-N, Matrix(I,r,r))

    N = (F0.Rs[F0.p]).*N
   
    #    N = N[:,reverse(invperm(F0.p))]
    N = N[reverse(invperm(F0.p)),:]
    
#    B = L[F.pcol[l+1:end]]

    return N, B
end
=#
function matrix(P::Vector, M::Dict)
    
    A = fill(zero(coefficients(P[1])[1]), length(P), length(M))
    j = 1
    for j in 1:length(P)
        p = P[j]
        m = monomials(p)
        c = coefficients(p)
        for k in 1:length(m)
            i = get(M, m[k], 0)
            if i != 0
                A[j,i] = c[k]
            end
        end
    end
    A
end

export sparse_matrix
function sparse_matrix(P::Vector, M::Dict)
    I = Int64[]
    J = Int64[]
    V = Float64[]
    
    for i in 1:length(P)
        for (m,c) in zip(monomials(P[i]),coefficients(P[i]))
            j = get(M, m, 0)
            if i != 0
                push!(I,i);
                push!(J,j), push!(V,c)
            end
        end
    end
    sparse(I,J,V)
end

function issmall(x)
    return abs(x) < 10^(-6)
end

function _rowechelon_basis(N, eps = 1.e-10)
    r = size(N,1)
    Idx = fill(0,r)
    for i in 1:r
        Idx[i] = findfirst(x->abs(x)>eps, N[i,:])
    end
    return Idx
end

export column_basis
function column_basis(N, eps = 1.e-3)
    r = size(N,1)
    Idx = fill(0,r)
    N0 = N[:,1:1]
    Idx[1]=1
    for i in 2:r
        j = Idx[i-1]
        while j <= size(N,2)
            s = N0\N[:,j]
            err =  norm(N0*s-N[:,j])
            #println("-- ", j," ", err)
            if err > eps
                Idx[i] = j
                N0 =hcat(N0,N[:,j])
                break
            end
            j += 1
        end
    end
    return Idx
end

function tnf_basis(N, L::AbstractVector, X)
    r = size(N,2)
    Idx = idx(L)
    I0 = Int[]
    for i in 1:length(L)
        I = map(t->get(Idx,t,0), [v*L[i] for v in X])
        if minimum(I)!=0
            push!(I0,i)
        end
    end
    N0 = transpose(N[I0, :])
    if size(N0,1) > size(N0,2)
        @error "Not enough monomials in the sum of the supports to compute a basis"
        return nothing
    end
    F = LinearAlgebra.qr(N0, Val(true))
    # println(I0, "   ",F.R)
    if issmall(F.R[r,r])
        @error "N0 not of maximal rank r=$r" F.R[r,r]
        return nothing
    end
    B = []
    for i in 1:size(N,2)
        push!(B,L[I0[F.p[i]]])
    end
    sort(B)
end

function mult_matrix(B, X, N, Nidx, ish = false)
    R = []
    Idx = idx(B)
    if !ish
        M = fill(0.0, length(B), size(N,2) )
        for t in Idx
            m = t.first
            i = t.second
            k = get(Nidx, m, 0)
            if k != 0
                for j in 1:size(N,2)
                    M[i,j] = N[k,j]
                end
            end
        end
        push!(R,M)
    end
    for v in X
        M = fill(0.0, length(B), size(N,2) )
        for t in Idx
            m = t.first
            i = t.second
            k = get(Nidx, m*v, 0)
            if k != 0
                for j in 1:size(N,2)
                    M[i,j] = N[k,j]
                end
            end
        end
        push!(R,M)
    end
    R
end


export multiplicities

function multiplicities(v)

    r =size(v,1)
    ms = [Int64[] for i in 1:r]

    Used = Dict{Int64, Bool}()
    for i in 1:r
        if !get(Used,i,false)
            push!(ms[i],i)
        end
        for j = i+1:r
            if !get(Used,j,false) && norm(v[j]-v[i])<5.e-3
                push!(ms[i],j)
                Used[j] = true
                #println("-> ", ms, " ", i, " ", j)
            end
        end
    end
    I = findall(m -> length(m)>0, ms)
    return ms[I]
end

export schur_dcp
function schur_dcp(M)
    lbd = randn(length(M))
    lbd /= LinearAlgebra.norm(lbd)

    M0 = sum(M[i]*lbd[i] for i in 1:length(M))
    
    T, Z, v = schur(ComplexF64.(M0))
    #println("... eig   ", t, "(s)"); t0=time()
    ms = multiplicities(v)
    #println("... ",ms)
    
    n = length(M)
    r = length(ms)

    Xi = fill(zero(eltype(v)),n,r)
    Tr = [Z'*M[i]*Z for i in 1:n]

    for k in 1:r
          for i in 1:n
            Xi[i,k] = 1/length(ms[k])*tr(Tr[i][ms[k],ms[k]])
        end
    end
    Xi, length.(ms), Z
end


function eigdiag(M)

    M0 = sum(M[i]*rand() for i in 1:length(M))

    #I0 = inv(M0)
    #Mg = I0*M[1]

    t = @elapsed E  = LinearAlgebra.eigvecs(M0)
    #println("... eig   ", t, "(s)"); t0=time()
    Z  = inv(E) #E\I0

    #t0 = time()
    #F = schurfact(Mg)
    #println("... schur ", time()-t0, "(s)"); t0=time()
    # E = F[:vectors]
    # Z = E'

    X = fill(eltype(E)(0),length(M),size(M0,1))
    for j in 1:length(M)
        Yj = Z*M[j]*E
        for i in 1:size(M0,1)
            X[j,i] = Yj[i,i] #(Y[:,i]\Yj[:,i])[1] #D[i,i]
        end
    end
    X
end



function (p::DynamicPolynomials.Polynomial)(x::Vector)
    return p(x...)
end

export res_matrix
"""
     R, L =  res_matrix(P::AbstractVector, M::AbstractVector)

Compute the matrix `R` of coefficients of the polynomials `P[i]*M[j]` and the liste `L` of monomials
which appear in these polynomials. A *row* of `R` represents a polynomial as a vector of coefficients with respect to the monomials in `L`.

"""
function res_matrix(P::AbstractVector, M::AbstractVector)
    X = DP.variables(P)
    RM = typeof(P[1])[]
    L = Set{typeof(DP.Monomial(X[1]))}()
    for i in 1:length(P)
        for m in M[i]
            RM = vcat(RM,P[i]*m)
            mns = DP.monomials(P[i]*m)
            L = union(L,Set(mns))
        end
    end
    L = sort(collect(L))
    sparse_matrix(RM, idx(L)), L
end

"""
Vector of relative errors of P at the points X
"""
function rel_error(P, Xi::Matrix, X = DP.variables(P))
    r = fill(0.0, length(P), size(Xi,2))
    n = size(Xi,2)

    for i in 1: size(Xi,2)
        for j in 1:length(P)
            V = Xi[:,i]
            r[j,i]= norm(DP.coefficients(DP.subs(P[j],X=>V)))
            r[j,i]/=norm(DP.coefficients(P[j]))
        end
    end
    r
end


function rel_error(P::Vector{AbstractAlgebra.Generic.MPoly{C}}, Xi::Matrix) where {C}
    P1 = as_polynomial_DP(P)
    return rel_error(P1, Xi)
end
