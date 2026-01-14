#AA = AbstractAlgebra
import AbstractAlgebra: gens, exponent_vectors

export matrix_mult, solve_groebner,
    as_polynomial, convert_coeff, reduce_by, as_monomial,
    prolong, border, normform

function as_monomial(p)
    DynamicPolynomials.monomials(p)[end]
end


function as_monomial(p, X)
    e = first(collect(exponent_vectors(p)))
    return prod(X[k]^e[k] for k in 1:length(e))
end

function as_polynomial(p, X, C::Type=Rational{BigInt})

    Exp   = collect(AbstractAlgebra.exponent_vectors(p))
    Coeff = [C(c) for c in AbstractAlgebra.coefficients(p)]
    
    Mon = [prod(X[k]^e[k] for k in 1:length(e)) for e in Exp]

    pol = sum(Coeff[i]*Mon[i] for i in 1:length(Coeff))

    return pol
end

function as_polynomial(p::DynamicPolynomials.Polynomial, R)

    Exp   = DynamicPolynomials.exponents.(DynamicPolynomials.monomials(p))
    Coeff = [c for c in DynamicPolynomials.coefficients(p)]

    X = gens(R)
    
    Mon = [prod(X[k]^e[k] for k in 1:length(e)) for e in Exp]

    pol = sum(Coeff[i]*Mon[i] for i in 1:length(Coeff))

end

function convert_coeff(p, C::Type=Float64)

    Mon = (monomials(p))
    Coeff = [C(c) for c in coefficients(p)]
    
    pol = sum(Coeff[i]*Mon[i] for i in 1:length(Coeff))

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

export Grobner
"""
Structure which defines

  - grobner_basis::Function
  - reduce:: Function
  -  quotient_basis:: Function
  -  ordering
"""
struct Grobner
    grobner_basis::Function
    reduce:: Function
    quotient_basis:: Function
end

function reduce_by(M,p,G)
    Groebner.normalform(G,p)
end

function normform(M, G::AbstractVector, B)
    r = length(B)
    Mnx = Dict{typeof(B[1]),Int64}([B[i] => i for i in 1:r]...)

    dB = border(B, variables(G))

    for i in 1:length(dB)  Mnx[dB[i]]=r+i  end

    N = fill(0.0, r, r+length(dB))

    for i in 1:r  N[i,i] = 1.0 end
    
    L  = [DynamicPolynomials.leading_monomial(g) for g in G]

    for g in G
        lm = DP.leading_monomial(g)
        mns = DP.monomials(g)
        cfs = DP.coefficients(g)
        i = Mnx[lm] 
        for k in 1:length(mns)-1
            if Mnx[mns[k]] >r
                println(">> ", lm, " ", mns[k], "   ",g)
            end
            
            if mns[k] != lm
                N[Mnx[mns[k]],i] = -cfs[k]
            end
        end
    end
        
    for i in 1:length(dB)
        alpha = dB[i]
        ir = last_divisible_by(alpha,L)
        mns = B*div(alpha,L[ir])
        cfs = N[:,Mnx[L[ir]]]
        N[:,r+i] = sum(cfs[j]*N[:,Mnx[mns[j]]] for j in 1:length(cfs))
        push!(L,alpha)
    end
    N, Mnx
end

"""
```
M = matrix_mult(M::Grobner, p, G::AbstractVector, Idx::Dict)
```
Compute the matrix of multication by `p` modulo `g` in the basis associated to the basis dictionary `Idx`. It is assumed that `g` is a Groebner basis and that the quotient is finite dimensional.

"""
function matrix_mult(Mth::Grobner, p, G::AbstractVector, Idx::Dict)
    delta = length(Idx)
    M = fill(zero(first(coefficients(G[1]))),delta,delta)        
    for key in Idx
        j = key.second
        m = key.first*p
        if (k = get(Idx,m,0)) == 0
            r = Mth.reduce(key.first*p,G)
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
M = mult_matrix(p, G::AbstractVector, B::AbstractVector)
```
Compute the matrix of multication by `p` modulo `G` in the basis `B`. It is assumed that `G` is a Groebner basis and that the quotient is finite dimensional.

"""
function mult_matrix(Mth::Grobner, p, G::AbstractVector, B)
    Idx = Dict{typeof(B[1]),Int64}([B[i] => i for i in 1:length(B)]...)
    return mult_matrix(Mth, p, G, Idx)
end


function mult_matrix(Mth, X, N, B::AbstractVector, Idx::Dict)
    M = typeof(N)[]
    for x in X
        Bx = B*x
        II = [Idx[m] for m in Bx]
        push!(M, N[:,II])
    end
    return M
end



"""
```
Xi, G, B = solve(Mth::Grobner, P::Vector; verbose = false)
```
Solve the system of polynomials `P`. It outputs:

 -  `Xi` the complex solution points, one per column of `Xi`
 -  `G` the computed Grobner basis
 -  `B` the basis of the quotient by the ideal of the equations

If `verbose = true`, the timing of the different steps is printed.

Example:
========
```
using AlgebraicSolvers, DynamicPolynomials, Groebner

Mth = Grobner((P,X) -> Groebner.groebner(P, ordering = Groebner.DegRevLex(X)),
              (p,G) -> Groebner.normalform(p,G),
              G -> Groebner.quotient_basis(G)
             )

X = @polyvar x y 

P = [-2+y-y^2+x^2*y, 1-3x+y+x*y^2]

Xi, G, B = AlgebraicSolvers.solve(Mth, P; verbose = true)


using AbstractAlgebra

R ,(x,y) =  polynomial_ring(AbstractAlgebra.QQ, [:x,:y])

P = [x^2*y-y^2, x*y^2-3]

Xi, G, B = AlgebraicSolvers.solve(Mth, P)

```

"""
function solve(Mth::Grobner, P::Vector{DynamicPolynomials.Polynomial{T,O,C}}; verbose=false) where {T,O,C}
    _solve_grobner_DP(Mth, P; verbose=verbose)
end

function solve(Mth::Grobner, P::AbstractVector; verbose=false)
    _solve_grobner_AA(Mth, P; verbose=verbose)
end
  
function _solve_grobner_DP(Mth::Grobner, P; verbose=false)
    
    X = DynamicPolynomials.variables(P)

    verbose && print("\033[96mComputing Grobner basis \033[0m")
    t = @elapsed G = Mth.grobner_basis(P,X)
    verbose && println(t, "s")

    verbose && print("\033[96mComputing quotient basis \033[0m")
    t = @elapsed B = sort(as_monomial.(Mth.quotient_basis(G))); 
    verbose && println(t, "s")

    Gf = [convert_coeff(g, Float64) for g in G]

    verbose && print("\033[96mComputing normal form    \033[0m")
    t = @elapsed N, II = normform(Mth, Gf, B)
    verbose && println(t, "s")

    verbose && print("\033[96mComputing mult matrices  \033[0m")
    t = @elapsed  M = mult_matrix(Mth, X, N, B, II)
    verbose && println(t, "s")
    
    verbose && print("\033[96mJoint diagonalisation    \033[0m")
    t = @elapsed Xi, E, Info = MultivariateSeries.diagonalization(M)
    verbose && println(t, "s\n")
    
    return Xi, G, B
end


function _solve_grobner_AA(Mth, P::AbstractVector; verbose=false)

    R = parent(P[1])
    n = length(AbstractAlgebra.gens(R))

    X = (DynamicPolynomials.@polyvar x[1:n] monomial_order = Graded{Reverse{InverseLexOrder}})[1]
    C = typeof(first(AbstractAlgebra.coefficients(P[1])))

    P1 = [as_polynomial(p, X, C) for p in P]
    
    Xi, G1, B = _solve_grobner_DP(Mth, P1, verbose=verbose) 
end

