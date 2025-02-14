import Groebner, MultivariateSeries, AbstractAlgebra, LinearAlgebra

export matrix_mult, quotient_basis, solve_groebner
    
function _is_divisible_by(m,L)
    max(AbstractAlgebra.is_divisible_by.(m,L)...)
end

"""
```
B, Idx = quotient_basis(g::AbstractVector)
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
    t = @elapsed p = prod((1+AbstractAlgebra.gen(R,i))^(D[i]-1) for i in 1:n)
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
M = matrix_mult(p, g::AbstractVector, Idx::Dict)
```
Compute the matrix of multication by `p` modulo `g` in the basis associated to the basis dictionary `Idx`. It is assumed that `g` is a Groebner basis and that the quotient is finite dimensional.

"""
function matrix_mult(p, G::AbstractVector, Idx::Dict)
    delta = length(Idx)
    M = fill(zero(first(AbstractAlgebra.coefficients(G[1]))),delta,delta)        
    for key in Idx
        j = key.second
        r = Groebner.normalform(G,key.first*p)
        for (c,m) in zip(AbstractAlgebra.coefficients(r),AbstractAlgebra.monomials(r))
            # println(c," ",m, " ", Idx[m])
            M[Idx[m],j] = c
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
using AbstractAlgebra

R, (x,y) = QQ["x","y"]

I = [x^2-y, x^2*y-4]

Xi, G, B = solve_groebner(I)
```
"""
function solve_groebner(P::AbstractVector; verbose=false)
    verbose && print("Computing Groebner basis ")
    t = @elapsed G = Groebner.groebner(P)
    verbose && println(t, "s")

    R = AbstractAlgebra.parent(G[1])
    verbose && print("Computing quotient basis ")
    t = @elapsed B, Bidx = quotient_basis(G)
    verbose && println(t, "s")

    verbose && print("Computing matrices of multiplication ")
    t = @elapsed M = [Float64.(matrix_mult(v, G, Bidx)) for v in AbstractAlgebra.gens(R)]
    verbose && println(t, "s")

    verbose && print("Joint triangularization ")
    t = @elapsed Xi, E, Info = MultivariateSeries.diagonalization(M)
    verbose && println(t, "s\n")

    Xi, G, B
end


