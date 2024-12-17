module AlgebraicSolvers

  using MultivariatePolynomials
  using DynamicPolynomials

  degree = DynamicPolynomials.maxdegree


#=
export coeftype
  coeftype(::Type{DynamicPolynomials.Polynomial{C, T}}) where {C, T} = T
  coeftype(p::DynamicPolynomials.Polynomial{C, T}) where {C, T} = T

  Base.one(X::Vector{DynamicPolynomials.Variable}) = monomials(X,0)[1]
=#

  include("mindex.jl")
  include("matrix.jl")
  include("macaulay.jl")
  include("newton.jl")
  include("toric.jl")
  include("CannyEmiris.jl")
#  include("groebner.jl")


#=
  function buildpolyvar(::Type{PV}, arg, var) where PV
    :($(esc(arg)) = $var)
  end

  export @Ring
  macro Ring(args...)
    X = DynamicPolynomials.PolyVar{true}[DynamicPolynomials.PolyVar{true}(string(arg)) for arg in args]
    V = [buildpolyvar(DynamicPolynomials.PolyVar{true}, args[i], X[i]) for i in 1:length(X)]
    push!(V, :(TMP = $X) )
    reduce((x,y) -> :($x; $y), V; init = :() )
  end
=#

end

