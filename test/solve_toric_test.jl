using AlgebraicSolvers
using DynamicPolynomials, LinearAlgebra

#=
export coeftype
  coeftype(::Type{Polynomial{C, T}}) where {C, T} = T
  coeftype(p::Polynomial{C, T}) where {C, T} = T

  Base.one(X::Vector{PolyVar{true}}) = monomials(X,0)[1]

include("../src/mindex.jl")
include("../src/matrix.jl")
include("../src/toric.jl")
include("../src/macaulay.jl")
=#

X = @polyvar x y

P = [1-3y+2x*y,x*y-x+y-1]


solve_toric(P)
