module AlgebraicSolvers

#  using MultivariatePolynomials
  import DynamicPolynomials

#  degree = DynamicPolynomials.maxdegree


  include("mindex.jl")
  include("matrix.jl")
  include("macaulay.jl")
  include("newton.jl")
  include("toric.jl")
  include("CannyEmiris.jl")
  include("groebner.jl")


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

