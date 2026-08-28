export alpha_beta, jacobian, newton_iter, newton_improve!

#=
function der(p::Polynomial{C,T}, v::PolyVar{C}) where {C,T}
    r = zero(p)
    X = variables(p)
    i = findfirst(v,X)
    if i>0
        for t in p
            if t.x.z[i] >0
                e = copy(t.x.z);
                d = e[i]
                e[i] = d-1;
                r+= Monomial(X,e)*t.α*d
            end
        end
    end
    r
end
=#

function jacobian(P,X)
    [differentiate(P[i],X[j]) for i in 1:length(P), j in 1:length(X)]
end

function newton_iter(F, J, X, Xi; verbose = false)
    J0 = fill(zero(Xi[1]), size(J,1), size(J,2))
    for i in 1:size(J,1)
        for j in 1:size(J,2)
            if J[i,j] != zero(J[i,j])
                J0[i,j] = subs(J[i,j], X=>Xi)
            end
        end
    end
    F0 = fill(zero(Xi[1]), length(F))
    for i in 1:length(F)
        if F[i] != zero(F[i])
                F0[i] = cst(subs(F[i], X=>Xi))
        end
    end
    err = norm(F0)
    Xi -= J0\F0
    Xi, err
end

"""
```
newton_improve!(Xi::Matrix, P, X=variables(P), eps::Float64=1.e-12, Nit::Int64 = 20)
``` 

Improve the roots `Xi` of the system `P` by Newton iteration.

 - `Xi` matrix of n x r roots where n is the number of coordinates of the roots and r the number of roots
 - `P` is the (square) system  of polynomials 
 - `X` the array of variables
 - `eps` threshold for stoping the iteration when the relative error is smaller.
 - `Nit` is the maximal number of iterations per root.

"""
function newton_improve!(Xi::Matrix, P, X=variables(P), eps::Float64=1.e-12, Nit::Int64 = 20;verbose = false)
    J = jacobian(P,X)
    if verbose
        L0err = [norm([cst(subs(p,X=>Xi[:,i])) for p in P]) for i in 1:size(Xi,2)]
        err0 = max(L0err...)
        println("\033[96m-- Init. max residual\033[0m: $err0")
    end
    Lerr = Float64[]
    for j in 1:size(Xi,2)
        i = 1
        err = 1.0
        V = Xi[:,j]
        while err>eps && i< Nit
            V, err = newton_iter(P, J, X, V)
            i+=1
        end
        if verbose && i==Nit && err>eps
            println("\033[93m nwt: stop iteration for root$j res.: $err  \033[0m",)
        end
        Xi[:,j] = V
        push!(Lerr,err)
    end
    mxerr= max(Lerr...)
    verbose && println("\033[96m-- Final max residual:\033[0m $(max(Lerr...))")
    Xi
end


"""
alpha, beta quantities for Newton convergence to an approximate zero.

 - If alpha < 0.125, then the approximate zero is within 2*beta from Xi and Newton methods converges to it from Xi quadratically.
 - If alpha < 0.02, then Newton method converges from all points in the ball of center Xi and radius 2*beta.

"""
function alpha_beta(P::Vector, Xi::Vector)
    X  = variables(P[1])
    for i in 2:length(P) X = union(X,variables(P[i])) end
    F0 = fill(zero(Xi[1]),length(P))
    for (p,i) in zip(P,1:length(P))
        f=subs(p,X=>Xi)
        if f != zero(f)
            F0[i]=(f[1]).α
        end
    end
    J0 = fill(zero(Xi[1]), length(P), length(X))
    for i in 1:length(P)
        for j in 1:length(X)
            dP = subs(differentiate(P[i],X[j]),X=>Xi)
            if dP != zero(dP)
                J0[i,j] = (dP[1]).α
            end
        end
    end
    beta = norm(J0\F0)

    s  =  sqrt(1 + norm(Xi)^2)
    mu =  sqrt(sum(norm(p,degree(p))^2 for p in P))
    mu *= norm(J0\diagm(0 => [sqrt(degree(p))*s^(degree(p)-1) for p in P]))

    d  = maximum([degree(p) for p in P])
    gamma = 0.5*d*sqrt(d)*mu

    beta*gamma, 2*beta
end


"""
    rel_error(P, Xi::Matrix, X = variables(P))

Vector of relative errors of P at the columns of the matrix `Xi`.
"""
function rel_error(P, Xi::Matrix)
    X = DP.variables(P)
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

"""
    rel_error(P, Xi::Vector, X = variables(P))

Vector of relative errors of P at the points Xi (which is a vector of vectors)
"""
function rel_error(P, Xi::Vector, X = variables(P))
    X = DP.variables(P)
    r = fill(0.0, length(P), length(Xi))
    n = length(Xi[1])

    for i in 1: length(Xi)
        for j in 1:length(P)
            V = Xi[i]
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
