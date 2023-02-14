# Numerical solutions

In order to improve the quality roots, computed by eigensolvers, one can
use Newton iterations for each root, assuming the roots have mutliplicity 1. Here are some functions to analyse and improve the numerical qauality of these roots. 
    
```@docs 
AlgebraicSolvers.rel_error
```

```@docs 
AlgebraicSolvers.alpha_beta
```

```@docs 
AlgebraicSolvers.newton_improve!
```
