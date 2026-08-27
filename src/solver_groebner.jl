#AA = AbstractAlgebra
import AbstractAlgebra: gens, exponent_vectors

export mult_matrices, solve_groebner, reduce_by, prolong, border, tnf


"""
Structure which defines

  - `grobner_basis :: Function (P) -> ` grobner basis of `P` with variables `X`
  - `reduce :: Function (p, G) -> ` normal form of the polynomial `p` modulo G 
  - `quot_basis :: Function P -> ` basis of the quotient by the ideal(`P`)

It can be used as an argument of [`solve`](@ref), [`tnf`](@ref), [`quot_basis`](@ref).

## Constructors:
```
    GB = GbSolver(gbasis, normform, qbasis)
```
Construct a structure of type `GbSolver` where

 - `gbasis` is the function which computes the GbSolver basis
 - `normform` is the function to reduce modulo a GbSolver basis
 - `qbasis` is the function which provides a basis of the quotient algebra.

```
    GB = GbSolver()
```
Construct a GbSolver. If a Groebner engine is available, from a package such as `Groebner`, it uses it otherwise it returns `nothing`.


## Example:
```
using AlgebraicSolvers, DynamicPolynomials, Groebner

GB = GbSolver(
        P->Groebner.groebner(P, ordering = Groebner.DegRevLex(variables(P))),
        Groebner.normalform,
        Groebner.quotient_basis
    )

X = @polyvar x y
P = [x-y+x^2,x-y+y^2]
Xi, ms = solve(P, GB)

```
"""
mutable struct GbSolver
    grobner_basis ::Function
    reduce ::Function
    quot_basis ::Function
    defined::Bool
end

export GbSolver


function GbSolver(gbasis, normform, qbasis)
    GbSolver(gbasis, normform, qbasis,true)
end

function GbSolver()
     if isdefined(Main,:Groebner)
         return GbSolver(P -> Main.Groebner.groebner(P, ordering = Main.Groebner.DegRevLex(variables(P))),
                         Main.Groebner.normalform,
                         Main.Groebner.quotient_basis
                         )
    else
        @error "No Groebner engine defined: add e.g. using Groebner"
        return  nothing #GbSolver(X->(), P->(), P->(), false)
    end

end



function check_gb(Mth::GbSolver)
 if isdefined(Main,:Groebner)
     if !Mth.defined
         Mth.grobner_basis = P ->Main.Groebner.groebner(P, ordering = Main.Groebner.DegRevLex(variables(P)))
         Mth.reduce = Main.Groebner.normalform
         Mth.quot_basis = Main.Groebner.quotient_basis
     end
     return Mth
    else
        @error "No Groebner engine available: add e.g. using Groebner"
        return nothing
    end
end

function reduce_by(M,p,G)
    Groebner.normalform(G,p)
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

#=
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
=#

export mult_matrix
"""
```
M = mult_matrix(p, G::AbstractVector, Idx::Dict, M::GbSolver)
```
Compute the matrix of multication by `p` modulo `g` in the basis associated to the basis dictionary `Idx`. It is assumed that `g` is a Groebner basis and that the quotient is finite dimensional.

"""
function mult_matrix(p, G::AbstractVector, Idx::Dict, Mth::GbSolver)
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


#==
```
M = _mult_matrix(p, G::AbstractVector, B::AbstractVector)
```
Compute the matrix of multication by `p` modulo `G` in the basis `B`. It is assumed that `G` is a Groebner basis and that the quotient is finite dimensional.

==#
function _mult_matrix(p, G::AbstractVector, B, Mth::GbSolver)
    Idx = Dict{typeof(B[1]),Int64}([B[i] => i for i in 1:length(B)]...)
    return mult_matrix(p, G, Idx, Mth)
end


function _mult_matrices(X, N, B::AbstractVector, Idx::Dict, Mth::GbSolver)
    M = typeof(N)[]
    for x in X
        Bx = B*x
        II = [Idx[m] for m in Bx]
        push!(M, N[:,II])
    end
    return M
end

export quot_basis
"""
```
    B = quot_basis(P, Mth::GbSolver; verbose = false)
```
Computes the basis `B` of the quotient by the ideal (P), formed by the monomials which are not in the inital of (P).
"""
function quot_basis(P, Mth::GbSolver; verbose = false)
    
    t = @elapsed G = Mth.grobner_basis(P)
    verbose && println("\033[96m-- Groebner basis \033[0m",t, "(s)")

    B = sort(as_monomial.(Mth.quot_basis(G))); 
end

"""
```
    B = tnf(P, Mth::GbSolver; verbose = false)
```
Computes the Truncated Normal Form on `B^+` where `B` is the basis of the quotient by the ideal (`P`), formed by the monomials which are not in the inital of (`P`).
"""
function tnf(P::AbstractVector, Mth::GbSolver; verbose = false)

    X = DynamicPolynomials.variables(P)
    t = @elapsed G = Mth.grobner_basis(P)
    verbose && println("\033[96m-- Groebner basis \033[0m",t, "(s)")

    B = sort(as_monomial.(Mth.quot_basis(G)))
    Gf = [convert_coeff(g, Float64) for g in G]
    t = @elapsed N, II = tnf_gb(Gf, B, Mth)
    verbose && println("\033[96m-- TNF ", size(N,1),"x",size(N,2)," \033[0m",t, "(s)")
    
    l = 1; m0 =1
    for (m,i) in II
        l = max(l,i)
    end
    L = fill(DynamicPolynomials.monomial(X[1]^0),l)
    for (m,i) in II
        L[i] = m
    end

    return N, L, [II[b] for b in B]
end

function tnf_gb(G::AbstractVector, B, Mth::GbSolver; verbose = false)
  
    r = length(B)
    Mnx = Dict{typeof(B[1]),Int64}([B[i] => i for i in 1:r]...)

    dB = border(B, union(DynamicPolynomials.variables.(G)...))

    for i in 1:length(dB)  Mnx[dB[i]]=r+i  end

    N = fill(0.0, r, r+length(dB))

    for i in 1:r  N[i,i] = 1.0 end
    
    L  = [DynamicPolynomials.leading_monomial(g) for g in G]

    for l in 1:length(G)
        g= G[l]
        lm = L[l]
        mns = DP.monomials(g)
        cfs = DP.coefficients(g)
        i = Mnx[lm] 
        for k in 1:length(mns)-1
            if Mnx[mns[k]] > r
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


#=
function _mult_matrices(P::AbstractVector, X, Mth::GbSolver)
    X = DynamicPolynomials.variables(P)
    G = Mth.grobner_basis(P)
    B = sort(as_monomial.(Mth.quot_basis(G)))
    Gf = [convert_coeff(g, Float64) for g in G]
    N, II = tnf(Gf, B, Mth)
    M = mult_matrices(X, N, B, II, Mth)
end
=#


function _solve(P::Vector{DynamicPolynomials.Polynomial{T,O,C}}, Mth::GbSolver; verbose=false) where {T,O,C}
    _solve_grobner_DP(P, Mth; verbose=verbose)
end

function _solve(P::AbstractVector, Mth::GbSolver; verbose=false)
    _solve_grobner_AA(P, Mth; verbose=verbose)
end
  
function _solve_grobner_DP(P, Mth::GbSolver; verbose=false)
    
    X = DynamicPolynomials.variables(P)

    t = @elapsed G = Mth.grobner_basis(P)
    verbose && println("\033[96m-- Groebner basis \033[0m",t, "(s)")

    t = @elapsed B = sort(as_monomial.(Mth.quot_basis(G))); 
    verbose && println("\033[96m-- Computing quotient basis ", length(B)," \033[0m",t, "(s)")

    Gf = [convert_coeff(g, Float64) for g in G]

    t = @elapsed N, II = tnf(Gf, B, Mth)
    verbose && println("\033[96m-- Computing normal form ", size(N,1), "x", size(N,2),"   \033[0m", t, "(s)")

    
    t = @elapsed  M = mult_matrices(X, N, B, II, Mth)
    verbose && println("\033[96m-- Computing mult matrices  \033[0m", t, "(s)")
    
    #t = @elapsed Xi = eigdiag(M)
    t = @elapsed Xi, ms = schur_dcp(M)
    verbose && println("\033[96m-- Schur  decomposition    \033[0m", t, "(s)")
    
    return Xi, ms, G, B, N
end


function _solve_grobner_AA(P::AbstractVector, Mth::GbSolver; verbose=false)

    R = parent(P[1])
    n = length(AbstractAlgebra.gens(R))

    X = (DynamicPolynomials.@polyvar x[1:n] monomial_order = Graded{Reverse{InverseLexOrder}})[1]
    C = typeof(first(AbstractAlgebra.coefficients(P[1])))

    P1 = [as_polynomial(p, X, C) for p in P]
    
    Xi, G1, B = _solve_grobner_DP(P1, Mth; verbose=verbose) 
end


