using AlgebraicSolvers
include("solver_dense.jl")

n = 3

Ld = [2,4,6,8,10]#,12,14,16,18,20] 
t  = [solver_dense(n,d) for d in Ld]

save("dense3d.res", length(X), Ld, t)




