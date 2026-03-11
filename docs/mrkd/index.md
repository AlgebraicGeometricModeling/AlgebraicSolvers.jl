# AlgebraicSolvers

The package provides several functions for computing in Artinian Algebras. 
These are the quotient of a polynomial ring by a zero-dimensional ideal. The main features include:

 - solution of polynomial systems, via multiplication operators and eigen computation.
 - resultant matrix constructions
 - Truncated Normal Form (TNF) computation 
 - duality and multivariate series
 - decomposition of polynomial exponential series
 - inverse systems

## Table of contents
```@contents
Pages = map(file -> joinpath("code", file), filter(x ->endswith(x, "md"), readdir("code"))) 
```


## Dependencies

  - `DataStructures .jl`
  - `DynamicPolynomials.jl`
  - `AbstractAlgebra.jl`
