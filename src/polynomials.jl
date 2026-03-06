export deg, monoms, exponent, matrixof, prodvec, prodset

import DynamicPolynomials: maxdegree, monomials

using LinearAlgebra
#import LinearAlgebra: norm, dot

function buildpolvar(::Type{PV}, arg, var) where PV
    :($(esc(arg)) = $var)
end

#=
"""
```
@ring args...
```
Defines the arguments as variables and output their array.

Example
-------
```
X = @ring x1 x2
```
"""
macro ring(args...)
    X = DynamicPolynomials.PolyVar{true}[DynamicPolynomials.PolyVar{true}(string(arg)) for arg in args]
    V = [buildpolvar(PolyVar{true}, args[i], X[i]) for i in 1:length(X)]
    push!(V, :(TMP = $X) )
    Base.reduce((x,y) -> :($x; $y), V; init = :() )
end
=#

#----------------------------------------------------------------------
#=
"""
```
deg(p:Polynomial) -> Int64
```
Degree of a polynomial
"""
function deg(p::DynamicPolynomials.Polynomial{B,T}) where {B,T}
    maxdegree(p.x)
end

#----------------------------------------------------------------------
function deg(t::Term{B,T})  where {B,T}
    deg(t.x)
end
#----------------------------------------------------------------------
function deg(m::Monomial{C}) where C
    sum(m.z)
end
#----------------------------------------------------------------------
function deg(v::DynamicPolynomials.PolyVar{T}) where T
    1
end
#----------------------------------------------------------------------
function MultivariatePolynomials.variables(m::Monomial{C}) where C
    m.vars
end
#----------------------------------------------------------------------
function coeff(t::Term{B,T}) where {B,T}
    exponent(t.α)
end
=#
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
