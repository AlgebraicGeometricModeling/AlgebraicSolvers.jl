# Solvers

Several types of solvers are available in the packages. 
They differ in the way the truncated normal form is computed.


## Solvers using resultant constructions

```@docs 
AlgebraicSolvers.solve_macaulay
AlgebraicSolvers.solve_toric
```

## Solvers using Groebner basis computation

```@docs 
AlgebraicSolvers.solve_groebner
```


# Resultant matrices

Different constructions of resultant matrices are available, including projective or Macaulay resultant matrices, toric resultant matrices.

```@docs 
AlgebraicSolvers.matrix_macaulay
AlgebraicSolvers.matrix_toric
```

# Truncated  Normal Forms

```@docs 
AlgebraicSolvers.tnf_macaulay
AlgebraicSolvers.tnf_toric
```

# Linear algebra in quotient 

```@docs 
AlgebraicSolvers.matrix_mult
AlgebraicSolvers.quotient_basis
```



