using LinearAlgebra, DynamicPolynomials

export invsys
export matrixof, ideal
import Base: (&), (|)



function invsys(s::Series{C, M}, X = variables(s), d::Int = maxdegree(s)) where {C, M <: AbstractMonomial }

    D = monomials(X,0:d)
    Ip = Series{C,M}[]
    z0 = zero(s)
    for m in D
        sn = m*s
        if sn != z0
            push!(Ip,sn)
        end
    end
    Ip
end




#----------------------------------------------------------------------
function row_echelon_with_pivots!(A::Matrix{T}, ɛ = T <: Union{Rational,Integer} ? 0 : 100*eps(norm(A,Inf))) where T
    nr, nc = size(A)
    pivots = Vector{Int64}()
    i = j = 1
    while i <= nr && j <= nc
        (m, mi) = findmax(abs.(A[i:nr,j]))
        mi = mi+i - 1
        if m <= ɛ
            if ɛ > 0
                A[i:nr,j] .= zero(T)
            end
            j += 1
        else
            for k=j:nc
                A[i, k], A[mi, k] = A[mi, k], A[i, k]
            end
            d = A[i,j]
            for k = j:nc
                A[i,k] /= d
            end
            for k = 1:nr
                if k != i
                    d = A[k,j]
                    for l = j:nc
                        A[k,l] -= d*A[i,l]
                    end
                end
            end
            append!(pivots,j)
            i += 1
            j += 1
        end
    end
    return A, pivots
end


function matrixof(L::Vector{Series{T,Mon}}, M::AbstractVector) where {T, Mon}
    Idx = Dict{Monomial{true}, Int64}()
    for (m,i) in zip(M,1:length(M)) Idx[m]=i end

    A = fill(zero(T), length(L), length(Idx))
    for i in 1:length(L)
        for (m,c) in L[i]
            A[i,Idx[m]] = c
        end
    end
    A
end

#----------------------------------------------------------------------
function invsys!(F::Vector,
                  B0::Vector,
                  IB::Vector,
                  D0::Vector{Series{C,M}},
                  ID::Vector{Series{C,M}},
                  Nu::Vector{Vector{Int64}},
                  X,
                  Xi,
                  theta=1.e-6;
                  verbose = false) where {C,M}
    
    n = length(X)

    MIdx=Dict{Vector{Int64},Int64}()

    J = Int64[]
    N = 0
    for j in 1:length(D0)
        for k in 1:n
            m = B0[j]*X[n+1-k]
            if !in(m,B0)
                N+=1;
                MIdx[[j,n+1-k]] = N
                push!(J,n*(j-1)+k)
            end
        end
    end

    A = fill(zero(C), 0, N)

    # commutation relations
    for k in 1:n
        for l in 1:n
            if l != k
                for m in B0 
                    a = fill(zero(C),N)
                    for i in 1:length(D0)
                        if (p = get(MIdx,[i,l],0)) !=0
                            a[p] = dot(X[k]*D0[i],m)
                        end
                        if (p = get(MIdx,[i,k],0)) !=0
                            a[p] = - dot(X[l]*D0[i],m)
                        end
                    end
                    if a != fill(zero(C),N)
                        A = vcat(A, a')
                    end
                end
            end
        end
    end

    
    # orthogonality to B0
    for m in B0
        a = [dot(ID[i], m) for i in J]
        if a != fill(zero(C),N)
            A = vcat(A, a')
        end
    end

    # A, piv = row_echelon_with_pivots!(A)
    # r = length(piv)
    
    # orthogonality to F
    for f in F
        a = [(ID[j]|f) for j in J]
        if a != fill(zero(C),N)
            A = vcat(A, a')
        end
    end

    if size(A,1) != 0
        A, piv = row_echelon_with_pivots!(A,theta)
        r = length(piv)
    end
    
    # if A is null
    if size(A,1) == 0 || r==0
        n = length(IB)
        h = 0 
        for i in 1:N
            m = IB[i]
            push!(B0,m)
            lambda = dual(m*one(C))
            push!(D0,lambda)
            h += 1
            j = length(D0)
            for k in 1:n
                m = B0[j]*X[n+1-k]
                Ilambda = integrate(lambda,X,n+1-k)
                if !in(m,B0)
                    push!(ID, Ilambda)
                    push!(IB, m)
                    #Idx[[j,n+1-k]] = length(IB)
                end
            end
        end
        return h
    end
    
    #print(eps(norm(A,Inf)), " ")
    #println("\nA: ",A, "\n",piv)
    A = inv(A[1:r,piv])*A[1:r,:]
    #println("\nA: ",A)

    npiv = filter(i->!in(i,piv), collect(1:N))
    #println("Pivots: ", piv, " ", r, " " , npiv)
    #println("mu: ", length(D0), " ", J[piv])
    
    # Basis of the null space of A
    for i in npiv
        lambda = ID[J[i]]- sum(A[j,i]*ID[J[piv[j]]] for j in 1:r)
        push!(D0, lambda)
        push!(B0, IB[J[i]])
        j = length(D0)
        for l in J[piv]
            push!(Nu,[j,n-(l-1)%n,div(l,n)+1])
        end
        for k in 1:n
            m = B0[j]*X[n+1-k]
            Ilambda = integrate(lambda,X,n+1-k)
            #println(">> ", lambda, " ", X[n+1-k], "-> ", Ilambda, " ", m)
            if !in(m,B0)
                push!(ID, Ilambda)
                push!(IB, m)
            end
        end
    end
    length(npiv)
 end

"""

    D,B,h,Nu = invsys(F::Vector, Xi::Vector)

Compute the inverse system of the polynomial system `F` at the point `Xi`. It outputs

   - `D` a basis of the inverse system as series in the dual basis of the monomial basis.
   - `B` the dual monomial basis
   - `h` the dimensions of the spaces in each degree
   - `Nu` the indices of the unknowns at each integration


"""
function invsys(F::Vector, Xi = fill(0,length(variables(F)));
                verbose = false, theta::Float64=1.e-6)
    X = variables(F)
    n = length(X)

    C = coefficienttype(F[1])
    C = promote_type(C, eltype(Xi), Rational{BigInt})
    B = [one(F[1]*one(C))]
    D = [dual(B[1])]
    ID = typeof(D[1])[]
    IB = []
    Nu = Vector{Int64}[] #Dict{Vector{Int64},Int64}()
    H =  Int64[]

    for k in 1:n
        m = B[1]*X[n+1-k]
        Ilambda = integrate(D[1],X,n+1-k)
        push!(ID, Ilambda)
        push!(IB, m)
    end
    mu = 0;
    h  = 1
    while h != 0
        verbose && print(h," ")
        push!(H,h)
        mu += h
        h = invsys!(F,B,IB,D,ID,Nu,X,Xi,theta; verbose = verbose)
    end
    verbose && println(": mu = ",mu)
    D, B, H, Nu
end

function diff(F,D)
    R = typeof(F[1])[]
    for d in D
        for f in F
            push!(R, d|f)
        end
    end
    R
end

function Jacobian(F,D,B,X,Xi,Nu)
    n = length(X)
    r = length(B)
    s = length(Nu)

    JF = zeros(n+s,n+s)

    for i in 1:s
        dD = typeof(D[1])[]
        for k in 1:r



        end
    end
end

function iter_newton!(F,I,X,Xi)
    B, D, h = invsyst(F,X,Xi,1.e-3)
    G = diff(F,D)[I]
    Gxi = map(f-> f(Xi,X), G)
    Jxi = map(f-> f(Xi,X), differentiate(G,X))
    println("Gxi: ",Gxi)

    Xi  -= Jxi\Gxi
    println("Gxi: ", map(f-> f(Xi,X), G))
    B, D, Xi
end
