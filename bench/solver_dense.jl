using DynamicPolynomials, AlgebraicSolvers, LinearAlgebra

function generic_dense(n, d)

    println("n=$n d=$d")
    X = (@polyvar x[1:n])[1]

    L = monomials(X,0:d)
    s = length(L)

    P = randn(n,s)*L
end


function solver_dense(n, d; verbose = true )

    println("n=$n d=$d")
    X = (@polyvar x[1:n])[1]

    L = monomials(X,0:d)
    s = length(L)

    P = randn(n,s)*L

    @time Xi, ms = solve(P,Macaulay(); verbose=true)
    verbose && println("-- Number of solutions: ", size(Xi,2),"    ",t1,"(s)" )
    Er = rel_error(P,Xi)
    println("-- Rel error: ", norm(Er,Inf))
    println()
    return t1
end

function save(file::String, n, dg, t)
    io=open(file,"w")
    println(io,"n=", length(X))
    println(io,"d=", dg)
    println(io,"t=", t)
    close(io)
    println("n=", length(X))
    println("d=", dg)
    println("t=", t)
end
