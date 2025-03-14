#AA = AbstractAlgebra

export matrix_mult, solve_groebner,
    convert_poly, convert_coeff, reduce_by,
    prolong, border, normform

function as_monomial(p)
    DP.monomials(p)[end]
end

function convert_poly(p, X, C::Type=Rational{BigInt})

    Exp   = collect(AbstractAlgebra.exponent_vectors(p))
    Coeff = [C(c) for c in AbstractAlgebra.coefficients(p)]
    
    Mon = [prod(X[k]^e[k] for k in 1:length(e)) for e in Exp]
    pol = dot(Coeff, Mon)

    return pol
end

function convert_coeff(p, C::Type=Float64)

    Mon = (monomials(p))
    Coeff = [C(c) for c in coefficients(p)]
    
    pol = dot(Coeff, Mon)

    return pol
end


function prolong(B,X)
    Bp = union(Set(B), [B*x for x in X]...)
    Bp = sort([m for m in Bp])
end

function border(B,X)
    Bp = union(Set(B), [B*x for x in X]...)
    Bb = setdiff(Bp,B)
    Bb = sort([m for m in Bb])
end

function first_divisible_by(m,L)
    findfirst(t->DynamicPolynomials.divides(t,m),L)
end

function last_divisible_by(m,L)
    findlast(t->DynamicPolynomials.divides(t,m),L)
end

function _reduced_by(p,G)
    L  = [DynamicPolynomials.leading_monomial(g) for g in G]
    m  = DynamicPolynomials.leading_monomial(p)
    ir = _is_divisible_by(m,L)
    r  = zero(p)
    while length(DynamicPolynomials.coefficients(p))>0
        m  = DynamicPolynomials.leading_monomial(p)
        ir = _is_divisible_by(m,L)
        if ir != nothing
            p -= DynamicPolynomials.leading_coefficient(p)*div(m,L[ir])*G[ir]
        else
            t = DynamicPolynomials.leading_term(p)
            p -= t
            r += t
        end
    end
    r
end

function reduce_by(p,G)
    Groebner.normalform(G,p)
end

function normform(G::AbstractVector, B)
    r = length(B)
    Mnx = Dict{typeof(B[1]),Int64}([B[i] => i for i in 1:r]...)

    dB = border(B, variables(G))

    for i in 1:length(dB)  Mnx[dB[i]]=r+i  end

    N = fill(0.0, r, r+length(dB))
    for i in 1:r  N[i,i] = 1.0 end
    
    L  = [DynamicPolynomials.leading_monomial(g) for g in G]

    for g in G
        m = DP.leading_monomial(g)
        mns = DP.monomials(g)
        cfs = DP.coefficients(g)
        i = Mnx[m] 
        for k in 1:length(mns)-1
            N[Mnx[mns[k]],i] = -cfs[k]
        end
    end
        
    for i in 1:length(dB)
        alpha = dB[i]
        ir = last_divisible_by(alpha,L)
        #println(i," ",alpha," ",ir, " ",L[ir], " ", div(alpha,L[ir]))
        mns = B*div(alpha,L[ir])
        cfs = N[:,Mnx[L[ir]]]
        N[:,r+i] = sum(cfs[j]*N[:,Mnx[mns[j]]] for j in 1:length(cfs))
        push!(L,alpha)
    end
    N, Mnx
end




"""
```
M = matrix_mult(p, G::AbstractVector, Idx::Dict)
```
Compute the matrix of multication by `p` modulo `g` in the basis associated to the basis dictionary `Idx`. It is assumed that `g` is a Groebner basis and that the quotient is finite dimensional.

"""
function matrix_mult(p, G::AbstractVector, Idx::Dict)
    delta = length(Idx)
    M = fill(zero(first(coefficients(G[1]))),delta,delta)        
    for key in Idx
        j = key.second
        m = key.first*p
        if (k = get(Idx,m,0)) == 0
            r = reduce_by(key.first*p,G)
            for (cr,mr) in zip(coefficients(r),monomials(r))
                #println("NF ",m,"  ", r)
                M[Idx[mr],j] = cr
            end
        else
            M[k,j] = one(M[1,1])
        end
    end
    M
end


"""
```
M = matrix_mult(p, G::AbstractVector, B::AbstractVector)
```
Compute the matrix of multication by `p` modulo `G` in the basis `B`. It is assumed that `G` is a Groebner basis and that the quotient is finite dimensional.

"""
function matrix_mult(p, G::AbstractVector, B)
    Idx = Dict{typeof(B[1]),Int64}([B[i] => i for i in 1:length(B)]...)
    return matrix_mult(p, G, Idx)
end



function matrix_mult(X, N, B::AbstractVector, Idx::Dict)
    M = typeof(N)[]
    for x in X
        Bx = B*x
        II = [Idx[m] for m in Bx]
        push!(M, N[:,II])
    end
    return M
end

#=
function roots(M::AbstractVector)
    n = length(M)
    r = size(M[1],1)
    Mrnd = sum(M[i]*randn(Float64) for i in 1:n)
    Sch  = LinearAlgebra.Schur{Complex}(LinearAlgebra.schur(Mrnd))
    E = Sch.vectors
    D = [E'*m*E for m in M]
    Xi = fill(zero(D[1][1,1]),n,r)
    for i in 1:n
        for j in 1:r
            Xi[i,j]=D[i][j,j]
        end
    end
    Xi
end
=#

"""
```
Xi, G, B = solve_groebner(P::Vector; verbose = false)
```
Solve the system of polynomials `P`. It outputs:

 -  `Xi` the complex solution points, one per column of `Xi`
 -  `G` the computed Groebner basis
 -  `B` the basis of the quotient by the ideal of the equations

Example:
========
```
using DynamicPolynomials, AlgebraicSolvers

X = @polyvar x y 

P = [x^2-y, x^2*y-4]

Xi, G, B = solve_groebner(P)

```
"""
function solve_groebner(P::AbstractVector; verbose=false)

    X = variables(P)
    
    verbose && print("\033[96mComputing Groebner basis \033[0m")
    t = @elapsed G = Groebner.groebner(P)
    verbose && println(t, "s")

    verbose && print("\033[96mComputing quotient basis \033[0m")
    t = @elapsed B = sort(as_monomial.(Groebner.quotient_basis(G)))
    verbose && println(t, "s")

    Gf = [convert_coeff(g,Float64) for g in G]

    verbose && print("\033[96mComputing normal form    \033[0m")
    t = @elapsed N, II = normform(Gf,B)
    verbose && println(t, "s")

    verbose && print("\033[96mComputing mult matrices  \033[0m")
    t = @elapsed  M = matrix_mult(X,N,B,II)
    verbose && println(t, "s")
    
    verbose && print("\033[96mJoint diagonalisation    \033[0m")
    t = @elapsed Xi, E, Info = MultivariateSeries.diagonalization(M)
    verbose && println(t, "s\n")

    Xi, G, B
end


