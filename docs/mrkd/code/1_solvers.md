# Solvers

Several types of solvers are available in the packages. 
They differ in the way the truncated normal form is computed.

```@docs 
AlgebraicSolvers.solve
```

## Solvers using resultant constructions

They corresponds to the following type
```@docs 
AlgebraicSolvers.Macaulay
AlgebraicSolvers.Toric
```

## Solvers using Groebner basis computation

```@docs 
AlgebraicSolvers.Grobner
```

## Quotient Algebra


Truncated Normal Forms, which kernel generates the ideal, are available:

```@docs 
AlgebraicSolvers.tnf
```

The operators of multiplications (by the variables or any polynomial) can be computed with the following function:
```@docs 
AlgebraicSolvers.mult_matrices
AlgebraicSolvers.mult_matrix
```

A basis of the quotient by the ideal (`P`) can be also computed:
```@docs 
AlgebraicSolvers.quo_basis

```









