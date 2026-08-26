module AlgebraicSolvers

  using DynamicPolynomials  #import DynamicPolynomials as DP
  import LinearAlgebra, AbstractAlgebra


  include("matrix.jl")
  include("convert.jl")


  include("series.jl")
  include("polynomials.jl")
  include("moments.jl")
  include("hankel.jl")
  include("invsys.jl")

  include("quotient.jl")
  include("diagonalisation.jl")

  include("newton.jl")

  include("decompose.jl")

  include("solve.jl")
  include("solver_macaulay.jl")
  include("solver_toric.jl")
  include("solver_groebner.jl")
  include("res_ce.jl")

end

