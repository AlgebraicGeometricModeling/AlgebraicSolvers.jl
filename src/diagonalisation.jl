
# norm of off diagonal terms of a square matrix
export norm_off
function norm_off(M)
    if size(M[1],1)>1
        return sqrt(sum(abs2(M[i,j]) + abs2(M[j,i]) for i in 1:size(M,1) for j in i+1:size(M,1)))
    else
        return 0.0
    end
end

function diagonalization_iter(D)
    n = size(D[1],1)
    s = length(D)
    
    X = fill(zero(D[1][1,1]),n,n)
    Y = fill(zero(D[1][1,1]),n,n)

    A = fill(zero(D[1][1,1]),s,2)
    b = fill(zero(D[1][1,1]),s)
    for i in 1:n
        for j in 1:n
            if i != j
                for k in 1:s
                    A[k,1] = D[k][i,i]
                    A[k,2] = D[k][j,j]
                    b[k]   = -D[k][i,j]
                end
                v = A\b
                X[i,j] =  v[1]
                Y[i,j] =  v[2]
            end
        end
    end
    for i in 1:n
        X[i,i]=1
        Y[i,i]=1
    end
    return X, Y
end


export diagonalization
"""
```
Xi, E, Info =  diagonalization(M::Vector{Matrix{C}}; maxiter::Int64 = 10, epsiter::Float64 = 1.e-3)
```

Compute the joint diagonalizaion of the matrices in `M` (assuming it exists).
`maxiter` is the maximum number of iterations and `epsiter` is the threshold error to stop the iterations.
It outputs

    - `Xi` the diagonals row by row, the ith row for the matrix  `M[i]`
    - `E`  the matrix of common eigenvectors so that `M[i] = E*diagm(Xi[i,:])*inv(E)`

"""
function diagonalization(M::Vector{Matrix{C}};
                         maxiter::Int64 = 10,
                         epsiter::Float64 = 1.e-3) where {C}
    n  = length(M)
    r  = size(M[1],1)

    N   = maxiter 
    eps = epsiter 

    M1 = sum(M[i]*randn(Float64) for i in 1:n)
    E  = eigvecs(M1)

    F  = inv(E)
    
    D  = vcat([Matrix{C}(I,r,r)],[F*M[i]*E for i in 1:length(M)])
    err = sum(norm_off.(D))
    delta = sum(norm.(D))
    #println("diag off: ", err)

    Info = Dict{String,Any}()
    Info["maxIter"] = N
    Info["epsIter"] = eps
    Info["diagErr"] = err
    
    nit = 0

    if err/delta > 5.e-2
        delta = err
        while nit < N && delta > eps
            err0 = err
            X,Y = diagonalization_iter(D)
            D = [Y*D[i]*X for i in 1:length(D)]
            E = E*X
            F = Y*F
            nit+=1
            err = sum(norm_off.(D))
            delta = err0-err
            #println("Off", nit,": ", err, "   delta: ", delta)
        end
        Info["d*"]= err
    end
    Info["nIter"] = nit
    
    Xi = fill(zero(E[1,1]),n,r)
    for i in 1:r
    	for j in 1:n
	    Xi[j,i] = D[j+1][i,i]/D[1][i,i]
            #Xi[j,i] =(E[:,i]\(M[j]*E[:,i]))[1]
	end
    end
    return Xi, E, Info
end


function eigdiag(M)

    M0 = sum(M[i]*rand() for i in 1:length(M))

    E  = LinearAlgebra.eigvecs(M0)

    Z  = inv(E)

    X = fill(eltype(E)(0),length(M),size(M0,1))
    for j in 1:length(M)
        Yj = Z*M[j]*E
        for i in 1:size(M0,1)
            X[j,i] = Yj[i,i] #(Y[:,i]\Yj[:,i])[1] #D[i,i]
        end
    end
    X, E
end


export multiplicities
"""
   ms = multiplicities(v::AbstractVector, eps = 1.e-5)

Group the values in v by clusters of radius < eps and output the list
of lists of indices of the clusters

## Example
julia> multiplicities([0.,3.,0.0000000001])
2-element Vector{Vector{Int64}}:
 [1, 3]
 [2]
"""
function multiplicities(v::AbstractVector, eps = 1.e-5)

    r =size(v,1)
    ms = [Int64[] for i in 1:r]

    Used = Dict{Int64, Bool}()
    for i in 1:r
        if !get(Used,i,false)
            push!(ms[i],i)
        end
        for j = i+1:r
            if !get(Used,j,false) && norm(v[j]-v[i])<eps
                push!(ms[i],j)
                Used[j] = true
            end
        end
    end
    I = findall(m -> length(m)>0, ms)
    return ms[I]
end

export schur_dcp

function schur_dcp(M::AbstractVector, eps::Float64=1.e-5)
    lbd = randn(length(M))
    lbd /= LinearAlgebra.norm(lbd)

    M0 = sum(M[i]*lbd[i] for i in 1:length(M))
    
    T, Z, v = schur(ComplexF64.(M0))
    #println("... eig   ", t, "(s)"); t0=time()
    ms = multiplicities(v, eps)
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
    Xi, ms, Z, Tr
end


