module AlgebraicSolvers

  using DynamicPolynomials
#import DynamicPolynomials as DP
  import LinearAlgebra, Groebner, MultivariateSeries,  AbstractAlgebra



  include("matrix.jl")
  include("newton.jl")

  include("solve_macaulay.jl")
  include("solve_toric.jl")
  include("CannyEmiris.jl")
  include("solve_groebner.jl")



end

