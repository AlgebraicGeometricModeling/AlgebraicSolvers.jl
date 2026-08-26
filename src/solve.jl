"""

    Xi, ms = solve(P, Mth; verbose = false, diag=:schur, refine=true)

 - `P` polynomial system
 - `Mth` class specifying the solver.

Solve the system of polynomials `P=[p1, ..., pn]`, using Sylvester matrix `res_matrix(P,Mth)`.

Output:

 - `Xi` the solutions as a matrix of points, one per column
 - `ms` the vector of their multiplicities.

Options: 

- `verbose=::Bool` (default `false`)  prints information related to the different substeps and their timing.

- `diag=::Symbol` (default `:schur`, other possibilities `:diag` )  specify the method for joint diagonalization of the operators of multiplication.
 
- `refine=::Bool` (default `true`) applies Newton iterations to refine the roots.


Example
-------
```
using AlgebraicSolvers, DynamicPolynomials

X = @polyvar x y
P = [2-x*y+x^2,y^2+x-2]
Xi, ms = solve(P, Macaulay())

using AlgebraicSolvers, DynamicPolynomials, Groebner

X = @polyvar x y
P = [x-y+x^2,x-y+y^2]
Xi, ms = solve(P, GbSolver())

```
"""
function solve(P::AbstractArray, Mth;
               diag::Symbol = :schur,  verbose::Bool = false, refine::Bool = true)
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
    verbose && println("\033[96m-- Mult matrices \033[0m $t(s)"); 


    if diag == :schur
        t = @elapsed Xi, ms = schur_dcp(M)
        verbose && println("\033[96m-- Schur dec",  "   \033[0m",t, "(s)");
        ms = length.(ms)
    elseif diag == :diag
        t = @elapsed  Xi, E = eigdiag(M)
        verbose && println("\033[96m-- Eigen diag",  "   \033[0m",t, "(s)");
        ms = fill(1,size(Xi,2))
    elseif diag == :jointdiag
        t = @elapsed Xi, _ = diagonalization(M)
        verbose && println("\033[96m-- Joint diag",  "   \033[0m",t, "(s)");
        ms = fill(1,size(Xi,2))
    end

    if refine
        t = @elapsed newton_improve!(Xi,P; verbose = verbose)
        verbose && println("\033[96m-- Newton improve   \033[0m",t, "(s)"); 
    end
    
    return Xi, ms
end

function solve(P::AbstractVector, mth::Symbol; verbose::Bool = false)
      solve(P,Val(mth);verbose=verbose)  
end

function solve(P::AbstractArray, Mth::Nothing; diag::Symbol = :schur,  verbose::Bool = false)
    @error "Unable to solve the system"
    return nothing, nothing
end
