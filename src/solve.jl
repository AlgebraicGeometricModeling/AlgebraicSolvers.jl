"""

    Xi, ms = solve(P, Mth; verbose = false)

 - `P` polynomial system
 - `Mth` class specifying the solver.

Solve the system of polynomials `P=[p1, ..., pn]`, using Sylvester matrix `res_matrix(P,Mth)`.
It outputs the solutions `Xi` as a matrix of points, one per column, and the vector of their multiplicities `ms`.

If the option `verbose=true` is set, it prints the different substeps and their timing.


Example
-------
```
using AlgebraicSolvers, DynamicPolynomials

X = @polyvar x y
P = [2-x*y+x^2,y^2+x-2]
Xi, ms = solve(P, Macaulay())

using AlgebraicSolvers, DynamicPolynomials, Groebner

X = @polyvar x y
P = [2-x*y+x^2,y^2+x-2]
Xi, ms = solve(P, GbSolver())

```
"""
function solve(P::AbstractArray, Mth; diag::Symbol = :schur,  verbose::Bool = false)
    X = variables(P)

    N, L, IB = tnf(P, Mth; verbose=verbose)

    if N == nothing
        return nothing, nothing 
    end
    t = @elapsed M = mult_matrices(N, L, IB, X)
    if length(M)==0
        @warn "Cannot find a good basis; Method $(typeof(Mth)) seems not adapted "
        return nothing, nothing
    end
    verbose && println("\033[96m-- Mult matrices \033[0m",t, "(s)"); t0 = time()


    if diag == :schur
        t = @elapsed Xi, ms = schur_dcp(M)
        verbose && println("\033[96m-- Schur dec",  "   \033[0m",t, "(s)"); t0 = time()

        return Xi, length.(ms)

    elseif diag == :diag
         t = @elapsed  Xi = eigdiag(M)
        verbose && println("\033[96m-- Eigen diag",  "   \033[0m",t, "(s)"); t0 = time()
        return Xi
    end
end

function solve(P::AbstractVector, mth::Symbol; verbose::Bool = false)
      solve(P,Val(mth);verbose=verbose)  
end

function solve(P::AbstractArray, Mth::Nothing; diag::Symbol = :schur,  verbose::Bool = false)
    @error "Unable to solve the system"
    return nothing, nothing
end
