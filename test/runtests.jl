using AlgebraicSolvers, DynamicPolynomials
F = filter(x ->endswith(x, "test.jl"), readdir("."))

for f in F
    include(f)
end
