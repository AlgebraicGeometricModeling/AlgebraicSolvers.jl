export deg, monoms, exponent, matrixof, prodvec, prodset

import DynamicPolynomials: maxdegree, monomials

using LinearAlgebra
#import LinearAlgebra: norm, dot

function buildpolvar(::Type{PV}, arg, var) where PV
    :($(esc(arg)) = $var)
end

#----------------------------------------------------------------------
function Base.one(::Type{DynamicPolynomials.Monomial{true}})
    Monomial{true}()
end

#----------------------------------------------------------------------
"""
```
exponent(m::Monomial) -> Array{Int64,1}
```
Get the exponent of a monomial as an array of Int64
"""
function Base.exponent(m::DynamicPolynomials.Monomial)
    return m.z
end

function Base.exponent(t::Term{B,T}) where {B,T}
    return t.x.z
end

# function DynamicPolynomials.monomial(m::Monomial)
#     return m
# end

# function DynamicPolynomials.monomial(t::Term{B,T}) where {B,T}
#     return t.x
# end
#----------------------------------------------------------------------
"""
```
 inv(m :: Monomial{true})
```
 return the inverse monomial with opposite exponents.
"""
function Base.inv(m::DynamicPolynomials.Monomial{true})
    Monomial(m.vars,-m.z)
end

function Base.inv(v::DynamicPolynomials.Variable{T}) where T
    inv(Monomial(v))
end

function inv!(m::DynamicPolynomials.Monomial{true})
    m.z=-m.z
end
#----------------------------------------------------------------------
function isprimal(m::DynamicPolynomials.Monomial{true})
    return !any(t->t<0, m.z)
end
#-----------------------------------------------------------------------
#=
"""
Evaluate a polynomial p at a point x;

## Example
```
julia> X = @polyvar x1 x2;

julia> p = x1^2+x1*x2;

julia> p([1.0,0.5])
1.5
```
""" 
function (p::DynamicPolynomials.Polynomial{B,T})(x::AbstractVector, X=variables(p)) where {B,T}
    return subs(p,[X[i]=>x[i] for i in 1:length(x)]...)

end
=#

#----------------------------------------------------------------------
function LinearAlgebra.norm(p::AbstractPolynomial, x::Float64)
    if (x == Inf)
        r = - Inf
        for t in p
            r = max(r, abs(t.α))
        end
    else
        r = Inf
        for t in p
            r = min(r, abs(t.α))
        end
    end
    r
end

function LinearAlgebra.norm(pol::AbstractPolynomial, p::Int64=2)
    r=sum(abs(t.α)^p for t in pol)
    exp(log(r)/p)
end

#----------------------------------------------------------------------

"""
    prodvec(X,Y)

Product-wise vector of X by Y ordered by row: [x1*y1, x1*y2, ..., x2*y1, ....] 
"""
function prodvec(X, Y)
    res = typeof(X[1])[]
    [x*y for x in X for y in Y]
end

"""
    prodsect(X,Y)

Product-wise vector of X by Y ordered by row with no repetition: 
    if xi*yj appears as xi'*yj' with i'<i or i'=i and j'<j it is not inserted.
"""
function prodset(X, Y)
    res = typeof(X[1])[]
    for x in X
        for y in Y
            m = x*y
            if !in(m,res) push!(res, m) end
        end
    end
    res
end

export vdm
"""
  Vandermonde matrix of polynomials in L evaluated at the points Xi
"""
function vdm(L::AbstractVector, Xi::AbstractMatrix)
    X = variables(L)
    r  = size(Xi,2)
    un = one(L[1])
    V = fill(zero(Xi[1,1]),length(L),r)
    for i in 1:r
        V[:,i] = coefficient.(subs.(L, X=> Xi[:,i]),un)
    end
    return V
end

function cst(c::Number)
    return c
end

function cst(p::AbstractPolynomial)
    if p==0 return 0
    else
       return coefficient(p, one(monomials(p)[1]))
    end
end
export cst

#----------------------------------------------------------------------
function random_dense(n, d, RND = randn )

    X = (@polyvar x[1:n])[1]

    L = monomials(X,0:d)
    s = length(L)

    P = RND(n,s)*L
end
export random_dense
