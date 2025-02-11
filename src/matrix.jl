export matrix, smatrix, tnf_basis, mult_matrix, eigdiag, kernel, rel_error

import DynamicPolynomials: coefficients, monomials

using LinearAlgebra
using SparseArrays

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

function nullspace(A::AbstractSparseMatrix)

    F = lufact(A')
    U = F[:U]
    r = 1
    while r<=min(size(A,1),size(A,2)) && abs(U[r,r])> 1.e-4
         r+=1
    end
    r-= 1
    L = F[:L]'

    L0= L[1:r,1:r]
    K = L[1:r,r+1:end]
    N = cat(1, - L0\K, eye(size(A,2)-r))

    P = fill(0, size(A,2))
    for i in 1:size(A,2) P[F[:p][i]]= i end
    return (F[:Rs] .* N[P,:])

end

function matrix(P::Vector, M::MonomialIdx)
    
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

function smatrix(P::Vector, M::MonomialIdx)
    I = Int64[]
    J = Int64[]
    T = coeftype(P[1])
    V = T[]
    for (p,j) in zip(P,1:length(P))
        for t in p
            i = get(M, t.x, 0)
            if i != 0
                push!(I,i); push!(J,j), push!(V,t.α)
            end
        end
    end
    sparse(J,I,V)
end

function issmall(x)
    return abs(x) < 10^(-6)
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
        @error "N0 not of maximal rank" F.R[r,r]
        return nothing
    end
    B = []
    for i in 1:size(N,2)
        push!(B,L[I0[F.p[i]]])
    end
    B
end

function mult_matrix(B, X, N, Nidx, ish = false)
    R = []
    Idx = idx(B)
    if !ish
        M = fill(0.0, length(B), size(N,2) )
        for (m,i) in Idx.terms
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
        for (m,i) in Idx.terms
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

function eigdiag(M)
    M0 = sum(M[i]*rand() for i in 1:length(M))
    #t0=time()
    I0 = inv(M0)
    #println("... inv   ", time()-t0, "(s)"); t0=time()
    Mg = I0*M[1]

    E  = LinearAlgebra.eigvecs(Mg)
    #println("... eig   ", time()-t0, "(s)"); t0=time()
    Z  = E\I0

    #t0 = time()
    #F = schurfact(Mg)
    #println("... schur ", time()-t0, "(s)"); t0=time()
    # E = F[:vectors]
    # Z = E'

    X = fill(Complex{Float64}(0.0),length(M),size(M0,1))
    for j in 1:length(M)
        Yj = Z*M[j]*E
        # D = Y\Yj
        for i in 1:size(M0,1)
            X[j,i]= Yj[i,i] #(Y[:,i]\Yj[:,i])[1] #D[i,i]
        end
    end
    X
end

function (p::DynamicPolynomials.Polynomial)(x::Vector)
    return p(x...)
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
            s = 0.0
            m = DP.monomials(P[j])
            c = DP.coefficients(P[j])
            s = dot(norm.(c),norm.(DP.subs(m, X => V)))
            r[j,i]/=s
        end
    end
    r
end
