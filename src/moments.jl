export moment, exp, log, series, dual, sparse_pol, sparse_decompose, binomial;


#----------------------------------------------------------------------
function Base.exp(A::Matrix{U}, pt::Vector{T}) where {U,T}
    N = size(A)
    R = fill(zero(eltype(pt)), N[1], N[2])
    for i  in 1:N[1]
        for j in 1:N[2]
            R[i,j] = pt[j]^A[i,j]
        end
    end
    R
end

#----------------------------------------------------------------------
function Base.log(A::Array{U}, pt::Vector{T}) where {U,T}
    N = size(A)
    R = fill(zero(eltype(pt)), N[1], N[2])
    for i  in 1:N[1]
        for j in 1:N[2]
            R[i,j] = real(log(A[i,j])/log(pt[i]))
        end
    end
    R
end

#------------------------------------------------------------------------
"""
```
dual(p::Polynomial) -> Series{C}
```
Compute the series associated to the polynomial p, replacing
the variables xi by their dual variables dxi. C is the type of coefficients 
of the polynomial p and M its type of monomials.
"""
function dual(p::AbstractPolynomial)

    C = coefficients(p)
    M = monomials(p)

    return series([M[i]=> C[i] for i in 1:length(C) ])
    
end

"""
```
binomial(p::Polynomial) -> Series{C}
```
Multi-index bonomial coefficients.
"""
function Base.binomial(d::Int64, alpha::Vector{Int64})
  r = binomial(d, alpha[1])
  for i in 2:length(alpha)
      d -= alpha[i-1]
      r *= binomial(d, alpha[i])
  end
  r
end

""""
```
dual(p::AbstractPolynomial, d:: Int64) -> Series{C}
```
Compute the series associated to the tensor p in degree d, using the apolar product.
The coefficients are of type `C`, which is the promoted type between the type of the coefficients of `p` and `Rational{Int}`.
The coefficients of the monomials `m` of  `p` are divided by the `binomial(d, exponents(m))`.
"""
function dual(p::AbstractPolynomial, d::Int)

    C = promote_type(coefficient_type(p),Rational{Int})
    c = coefficients(p)
    m = monomials(p)

    return series([m[i]=> C(c[i])/Base.binomial(d,exponents(m[i])) for i in 1:length(c) ])
end

"""
Recover the polynomial from the series by duality.
"""
function dual(s::Series)
    dot(coefficients(s), monomials(s))
end

"""
Recover the polynomial from the series, using the apolar duality i.e. multiplying the coefficients by a binomial.
"""
function dual(s::Series, d::Int)
    Mm = monomials(s)
    Cf = coefficients(s)
    return sum([m*(c*Base.binomial(d,exponents(m))) for (c,m) in zip(Cf,Mm)])
end

#----------------------------------------------------------------------
"""
```
moment(w::Vector{C}, P::Matrix{C}) -> Vector{Int64} -> C
```
Compute the moment function ``α -> ∑_{i} ω_{i} P_{i}^α``
associated to the sequence P of r points of dimension n, which is a matrix
of size r*n and the weights w.
"""
moment(w, P) = function(α)
  res = 0
  for j in 1:length(w)
      m = 1
      for i in 1:length(α)
        m *= P[i,j]^α[i]
      end
      res+=m*w[j];
  end
  res
end

#------------------------------------------------------------------------
"""
```
series(f::Function, L::Vector{M}) -> Series{C,M}
```
Compute the generating series ``\\sum_{x^{α} \\in L} f(α) z^α``
for a function  ``f: \\mathbb{N}^n \\rightarrow C`` and a sequence L of monomials.
"""
function series(f::Function, L::AbstractVector)

   return series([m => f(exponents(m)) for m in L])

end



#----------------------------------------------------------------------
"""
```
moment(p::Polynomial, zeta::Vector{C}) -> Vector{Int64} -> C
```
Compute the moment function ``α \\rightarrow p(ζ^α)``.
"""
moment(p::AbstractPolynomial, zeta::AbstractVector) =  function(V::Vector{Int})
    U = [zeta[i]^V[i] for i in 1:length(V)]
    p(U)
end

#----------------------------------------------------------------------
"""
```
series(w:: Vector{C}, P::Matrix{C}, L::Vector{M}) -> Series{C,M}
```
Compute the series of the moment sequence ``∑_{i} ω_{i} P_{i}^α`` for ``α \\in L``.
"""
function series(w:: AbstractVector, P::AbstractMatrix, L::AbstractVector) 
   series(moment(w,P), L)
end

#----------------------------------------------------------------------
"""
```
series(w::AbstractVector, P::AbstractMatrix, X, d::Int64) -> Series{C,M}
```
Compute the series of the moment sequence ``∑_i ω_{i} P_{i}^α`` for ``|α| \\leq d``.
"""
function series(w::AbstractVector, P::AbstractMatrix, X, d::Int64) 
    h = moment(w,P)
    L = monomials(X,0:d)
   series(h,L)
end


#----------------------------------------------------------------------
"""
```
series(p::Polynomial, zeta, X, d::Int64) -> Series
```
Compute the series of moments ``p(ζ^α)`` for ``|α| \\leq d``.
"""
function series(p::DynamicPolynomials.Polynomial, zeta, X, d::Int64)
    h = moment(p,zeta)
    L = monomials(X,seq(0:d))
    series(h,L)
end

#----------------------------------------------------------------------
"""
```
sparse_pol(w, E, X) -> Polynomial{true,C}
```
Compute the polynomial ``∑ ωᵢ X^E[i,:]`` with coefficients ``ωᵢ`` and monomial exponents ``E[i,:]``.
"""
function sparse_pol(w::Vector, E::Matrix{Int64}, X)
    sum(w[j]*prod(X[i]^E[j,i] for i in 1:size(E,2)) for j in 1:length(w))
end



#------------------------------------------------------------------------
function sparse_decompose(f, zeta, X, d:: Int64)
    sigma = series(moment(f,zeta), monomials(X,seq(0:d)))
    w, Xi = decompose(sigma)
    w, log(Xi,zeta)
end

