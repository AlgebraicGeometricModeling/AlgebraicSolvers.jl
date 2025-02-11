import Groebner, MultivariateSeries, AbstractAlgebra, DynamicPolynomials, LinearAlgebra

export multiplication, quotient_basis, solve_groebner
    
function _is_divisible_by(m,L)
    max(AbstractAlgebra.is_divisible_by.(m,L)...)
end

"""
```
B, Idx = quotient_bais(g::AbstractVector)
```
Compute the monomial basis of the quotient by `g`, as the complementary of the initial ideal of `g`, assuming that `g` is a Groebner basis and it is finite.
It outputs: 

   - `B` the basis of monomlials  in increasing monomial order.
   - `Idx` the dictionary associting to a monomial its index in the basis.   

"""
function quotient_basis(g::AbstractVector)
    R = AbstractAlgebra.parent(g[1])
    n = AbstractAlgebra.number_of_variables(R)
    L = AbstractAlgebra.leading_monomial.(g)
    D = [max(AbstractAlgebra.degree.(L,i)...) for i in 1:n]
    p = prod((1+AbstractAlgebra.gen(R,i))^(D[i]-1) for i in 1:n)
    B0 = reverse([m for m in AbstractAlgebra.monomials(p)])
    B  = typeof(B0[1])[]
    Idx = Dict{typeof(B0[1]),Int64}()
    c = 1
    for m in B0
        if !_is_divisible_by(m,L)
            push!(B,m)
            Idx[m] = c
            c += 1
        end
    end
    B, Idx
end

"""
```
M = multiplication(p, g::AbstractVector, Idx::Dict)
```
Compute the matrix of multiplication by `p` modulo `g` in the basis associated to the basis dictionary `Idx`. It is assumed that `g` is a Groebner basis and that the quotient is finite dimensional.

"""
function multiplication(p, G::AbstractVector, Idx::Dict)
    delta = length(Idx)
    M = fill(zero(first(AbstractAlgebra.coefficients(G[1]))),delta,delta)        
    for key in Idx
        j = key.second
        r = normal_form(key.first*p, G)
        for (c,m) in zip(AbstractAlgebra.coefficients(r),AbstractAlgebra.monomials(r))
            # println(c," ",m, " ", Idx[m])
            M[Idx[m],j] = c
        end
    end
    M
end


function roots(M::AbstractVector)
    n = length(M)
    r = size(M[1],1)
    Mrnd = sum(M[i]*randn(Float64) for i in 1:n)
    Sch  = Schur{Complex}(LinearAlgebra.schur(Mrnd))
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

"""
```
Xi, G, B = solve_groebner(P::Vector)
```
Solve the system of polynomials `P`. It outputs:

    -  `Xi` the complex solutions, one per colmun of `Xi`
    -  `G` the computed Groebner basis
    -  `B` the basis dictionary for the quotient by the ideal of the equations

Example:
========
```
using AbstractAlgebra

R, (x,y) = QQ["x","y"]

I = [x^2-y, x^2*y-4]

Xi, G, B = solve_groebner(I)
```
"""
function solve_groebner(I::AbstractVector)
    G = Groebner.groebner(I)
    R = AbstractAlgebra.parent(G[1])
    B, Bidx = quotient_basis(G)
    M = [Float64.(multiplication(v, G, Bidx)) for v in AbstractAlgebra.gens(R)]
    Xi, E, Info = MultivariateSeries.diagonalization(M)
    Xi, G, B
end


