include("CannyEmiris.jl")
Main.CannyEmiris

A = [[1,1] [1,1] [1,1]]

H = [[1,0] [0,1]]

CE, PM = CannyEmiris.Zonotopes(A,H)
