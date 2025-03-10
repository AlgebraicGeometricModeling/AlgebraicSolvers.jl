# AlgebraicSolvers

Package for the solution of polynomial systems, using eigen computation.
It outputs all the complex solutions if the system is zero-dimensional.

It proceeds as follows:

 - Compute a Truncated Normal Form (using either resultant constructions or Groebner basis computation)
 - Compute operators of multiplication by the variables in a basis of the quotient algebra.
 - Compute a joint triangularization of the multplication matrices to obtain the (complex) solutions, with their multiplicity.
 

```@contents
Pages = map(file -> joinpath("code", file), filter(x ->endswith(x, "md"), readdir("code"))) 
```


## Dependencies

For the solver based on Groebner basis computation, the package `Grobner.jl` is used.

