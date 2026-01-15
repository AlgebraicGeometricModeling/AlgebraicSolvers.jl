module AlgebraicSolvers

  using DynamicPolynomials  #import DynamicPolynomials as DP
  import LinearAlgebra, MultivariateSeries,  AbstractAlgebra

  include("matrix.jl")
  include("quotient.jl")
  include("newton.jl")

  include("solve_macaulay.jl")
  include("solve_toric.jl")
  include("solve_groebner.jl")
  include("res_ce.jl")

end

