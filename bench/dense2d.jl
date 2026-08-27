include("solver_dense.jl")


n = 2
Ld = [2,3,4,5,6,7,8,9,10,20]#,30,40,50,60,70,80,90,100] 

t = [solver_dense(n,d) for d in Ld]

save("dense2d.res", n, Ld,t)
