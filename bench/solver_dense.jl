using DynamicPolynomials, AlgebraicSolvers, LinearAlgebra

function solver_dense(n, d; verbose = true )

    println("n=$n d=$d")

    P = random_dense(n, d)

    t = @elapsed Xi, ms = solve(P,Macaulay(); verbose=verbose)

    Er = norm(rel_error(P,Xi),Inf)
    println(">> n=$n d=$d   $t(s) Rel error: $Er")
    println()
    return t
end

function save(file::String, n, dg, t)
    io=open(file,"w")
    println(io,"n=", n)
    println(io,"d=", dg)
    println(io,"t=", t)
    close(io)
    println("n=", n)
    println("d=", dg)
    println("t=", t)
end
