using Documenter, AlgebraicSolvers

dir="mrkd"
Expl = map(file -> joinpath("expl", file), filter(x ->endswith(x, "md"), readdir(dir*"/expl")))
Code = map(file -> joinpath("code", file), filter(x ->endswith(x, "md"), readdir(dir*"/code")))

makedocs(
         sitename = "AlgebraicSolvers",
         authors = "B. Mourrain",
         modules = [AlgebraicSolvers],
         build = "AlgebraicSolvers.jl/docs",
         source = dir,
         pages = Any[
                     "Home" => "index.md",
                     "Functions" => Code,
                     "Examples" => Expl,             
         ],
         repo = Remotes.GitHub("AlgebraicGeometricModeling", "AlgebraicSolvers.jl"),
         doctest = false
         )

deploydocs(
           repo = Remotes.GitHub("AlgebraicGeometricModeling", "AlgebraicSolvers.jl"),
           target = "site"
           )

