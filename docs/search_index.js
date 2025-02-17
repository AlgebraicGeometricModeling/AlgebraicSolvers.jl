var documenterSearchIndex = {"docs":
[{"location":"code/3_functions/#Numerical-solutions","page":"Numerical solutions","title":"Numerical solutions","text":"","category":"section"},{"location":"code/3_functions/","page":"Numerical solutions","title":"Numerical solutions","text":"In order to improve the quality roots, computed by eigensolvers, one can use Newton iterations for each root, assuming the roots have mutliplicity 1. Here are some functions to analyse and improve the numerical qauality of these roots. ","category":"page"},{"location":"code/3_functions/","page":"Numerical solutions","title":"Numerical solutions","text":"AlgebraicSolvers.rel_error","category":"page"},{"location":"code/3_functions/#AlgebraicSolvers.rel_error","page":"Numerical solutions","title":"AlgebraicSolvers.rel_error","text":"Vector of relative errors of P at the points X\n\n\n\n\n\n","category":"function"},{"location":"code/3_functions/","page":"Numerical solutions","title":"Numerical solutions","text":"AlgebraicSolvers.alpha_beta","category":"page"},{"location":"code/3_functions/#AlgebraicSolvers.alpha_beta","page":"Numerical solutions","title":"AlgebraicSolvers.alpha_beta","text":"alpha, beta quantities for Newton convergence to an approximate zero.\n\nIf alpha < 0.125, then the approximate zero is within 2*beta from Xi and Newton methods converges to it from Xi quadratically.\nIf alpha < 0.02, then Newton method converges from all points in the ball of center Xi and radius 2*beta.\n\n\n\n\n\n","category":"function"},{"location":"code/3_functions/","page":"Numerical solutions","title":"Numerical solutions","text":"AlgebraicSolvers.newton_improve!","category":"page"},{"location":"code/3_functions/#AlgebraicSolvers.newton_improve!","page":"Numerical solutions","title":"AlgebraicSolvers.newton_improve!","text":"newton_improve!(Xi::Matrix, P, X=variables(P), eps::Float64=1.e-12, Nit::Int64 = 20)\n\nImprove the roots Xi of the system P by Newton iteration.\n\nXi matrix of n x r roots where n is the number of coordinates of the roots and r the number of roots\nP is the (square) system  of polynomials \nX the array of variables\neps threshold for stoping the iteration when the relative error is smaller.\nNit is the maximal number of iterations per root.\n\n\n\n\n\n","category":"function"},{"location":"code/2_resultants/#Resultants","page":"Resultants","title":"Resultants","text":"","category":"section"},{"location":"code/2_resultants/","page":"Resultants","title":"Resultants","text":"Greedy implementation of the Canny-Emiris formula as it was presented in \"A Greedy Approach to the Canny-Emiris Formula\" with the introduction of type functions. The underlying system of equations is written as:","category":"page"},{"location":"code/2_resultants/","page":"Resultants","title":"Resultants","text":"F_i = sum_a in mathcalA_iu_iachi^a quad i = 0dotsn","category":"page"},{"location":"code/2_resultants/","page":"Resultants","title":"Resultants","text":"for some finite sets of supports mathcalA_i subset M in a lattice of rank n. This implementation treats the cases in which the Newton polytopes  Delta_i = conv(mathcalA_i) are zonotopes (sums of line segments) or products of simplices (which correspond to multihomogeneous systems of equations).","category":"page"},{"location":"code/2_resultants/#Zonotopes","page":"Resultants","title":"Zonotopes","text":"","category":"section"},{"location":"code/2_resultants/","page":"Resultants","title":"Resultants","text":"Let v_1dotsv_n in M be independent vectors generating an n-zonotope. The supports of our zonotopes are: ","category":"page"},{"location":"code/2_resultants/","page":"Resultants","title":"Resultants","text":"mathcalA_i^ = big sum_j = 1^n lambda_j v_j in mathbbZ^n mid quad lambda_j in mathbbZ quad 0  leq lambda_j leq a_ijbig","category":"page"},{"location":"code/2_resultants/","page":"Resultants","title":"Resultants","text":"for some a_ij quad i = 0dotsn quad j = 1dotss. Let V be the n times n matrix whose columns are v_1dotsv_n and let A be the matrix whose entries are the a_ij.","category":"page"},{"location":"code/2_resultants/","page":"Resultants","title":"Resultants","text":"The function CannyEmiris.Zonotopes receives these two matrices and returns two symbolic matrices mathcalH_mathcalGmathcalE_mathcalG which correspond to the rational formula for the sparse resultant:","category":"page"},{"location":"code/2_resultants/","page":"Resultants","title":"Resultants","text":"Res_mathcalA = big(fracdet(mathcalH_mathcalG)det(mathcalE_mathcalG)big)^det(V)","category":"page"},{"location":"code/2_resultants/","page":"Resultants","title":"Resultants","text":"which follows from Theorem 1.1 and Corollary 3.1 on the text. ","category":"page"},{"location":"code/2_resultants/","page":"Resultants","title":"Resultants","text":"Moreover, the program prints the polynomials chi^b-a(b)F_i(b) corresponding to all the lattice points in the greedy subset b in mathcalG, the size of the matrix and the degree of the resultant (which corresponds to the lattice points in mixed cells).","category":"page"},{"location":"code/2_resultants/","page":"Resultants","title":"Resultants","text":"Let's see an example of the use of this function which corresponds to the system of Example 1.1 in the text.","category":"page"},{"location":"code/2_resultants/","page":"Resultants","title":"Resultants","text":"\njulia> Using AlgebraicSolvers\n\njulia> A = [[1,1] [1,1] [1,1]]\n2×3 Matrix{Int64}:\n 1  1  1\n 1  1  1\n\njulia> H = [[1,0] [0,1]]\n2×2 Matrix{Int64}:\n 1  0\n 0  1\n\njulia> CE, PM = CannyEmiris.Zonotopes(A,H)\nThe rows of the Canny-Emiris matrix x^{b-a(b)}F_{i(b)} are: \n[0, 1]-> x^[0, 1]*F_2\n[0, 2]-> x^[0, 1]*F_1\n[1, 0]-> x^[1, 0]*F_2\n[1, 1]-> x^[1, 1]*F_2\n[1, 2]-> x^[0, 1]*F_0\n[2, 0]-> x^[1, 0]*F_1\n[2, 1]-> x^[1, 0]*F_0\n[2, 2]-> x^[1, 1]*F_1\n\nThe size of the greedy Canny-Emiris matrix is: 8\nThe degree of the resultant is: 6\n\nThe sparse resultant is the ratio of the determinants of the returned matrices to the power 1.0\n\n\njulia> CE\n8×8 Matrix{SymPy.Sym}:\n (u_{2, [0, 0]})  (u_{2, [0, 1]})                0  (u_{2, [1, 0]})  (u_{2, [1, 1]})                0                0                0\n (u_{1, [0, 0]})  (u_{1, [0, 1]})                0  (u_{1, [1, 0]})  (u_{1, [1, 1]})                0                0                0\n               0                0  (u_{2, [0, 0]})  (u_{2, [0, 1]})                0  (u_{2, [1, 0]})  (u_{2, [1, 1]})                0\n               0                0                0  (u_{2, [0, 0]})  (u_{2, [0, 1]})                0  (u_{2, [1, 0]})  (u_{2, [1, 1]})\n (u_{0, [0, 0]})  (u_{0, [0, 1]})                0  (u_{0, [1, 0]})  (u_{0, [1, 1]})                0                0                0\n               0                0  (u_{1, [0, 0]})  (u_{1, [0, 1]})                0  (u_{1, [1, 0]})  (u_{1, [1, 1]})                0\n               0                0  (u_{0, [0, 0]})  (u_{0, [0, 1]})                0  (u_{0, [1, 0]})  (u_{0, [1, 1]})                0\n               0                0                0  (u_{1, [0, 0]})  (u_{1, [0, 1]})                0  (u_{1, [1, 0]})  (u_{1, [1, 1]})\n\njulia> PM\n2×2 Matrix{SymPy.Sym}:\n (u_{2, [0, 0]})  (u_{2, [1, 1]})\n (u_{1, [0, 0]})  (u_{1, [1, 1]})","category":"page"},{"location":"code/2_resultants/","page":"Resultants","title":"Resultants","text":"The program will return an error if i) the matrix A or H do not have the desired dimensions, ","category":"page"},{"location":"code/2_resultants/","page":"Resultants","title":"Resultants","text":"julia> A = [[1,1] [1,1]]\n2×2 Matrix{Int64}:\n 1  1\n 1  1\n\njulia> H\n2×2 Matrix{Int64}:\n 1  0\n 0  1\n\njulia> CE, PM = CannyEmiris.Zonotopes(A,H)\nThe matrix of the a_{i,j} does not have the correct dimensions\n(Dict{Any, Any}(), Dict{Any, Any}())","category":"page"},{"location":"code/2_resultants/","page":"Resultants","title":"Resultants","text":"ii) if H does not correspond to an n-zonotope","category":"page"},{"location":"code/2_resultants/","page":"Resultants","title":"Resultants","text":"julia> H = [[1,0] [0,0]]\n2×2 Matrix{Int64}:\n 1  0\n 0  0\n\njulia> CE, PM = CannyEmiris.Zonotopes(A,H)\nThe vectors do not correspond to an n-zonotope.\n(Dict{Any, Any}(), Dict{Any, Any}())","category":"page"},{"location":"code/2_resultants/","page":"Resultants","title":"Resultants","text":"or iii) entries of A are not ordered in the sense  0  a_0j leq a_1j leq dots leq a_n-1j quad j = 1dotsn:","category":"page"},{"location":"code/2_resultants/","page":"Resultants","title":"Resultants","text":"julia> A = [[2,1] [1,1] [1,1]]\n2×3 Matrix{Int64}:\n 2  1  1\n 1  1  1\n\njulia> H = [[1,0] [0,1]]\n2×2 Matrix{Int64}:\n 1  0\n 0  1\n\njulia> CE, PM = CannyEmiris.Zonotopes(A,H)\nThe matrix of the a_{i,j} does not satisfy a_{i-1,j} <= a_{i,j}\n(Dict{Any, Any}(), Dict{Any, Any}())","category":"page"},{"location":"code/2_resultants/#Multihomogeneous-systems","page":"Resultants","title":"Multihomogeneous systems","text":"","category":"section"},{"location":"code/2_resultants/","page":"Resultants","title":"Resultants","text":"Let now N represent the vector of n1,...ns in a multihomogeneous system in mathbbP^n_1 times dots times mathbbP^n_s and let D be a matrix whose columns are the multidegrees of the polynomials of the system so that:","category":"page"},{"location":"code/2_resultants/","page":"Resultants","title":"Resultants","text":"mathcalA_i = big(b_jl)_l = 1dotss^j = 1dotsn_l in oplus_l=1^smathbbZ^n_l  b_jl geq 0 quad sum_j=0^n_lb_jl leq d_il  quad i = 0dotsn","category":"page"},{"location":"code/2_resultants/","page":"Resultants","title":"Resultants","text":"The function CannyEmiris.Multihomogeneous(D::Matrix,N::Vector) receives these two matrices D and N and returns and returns two symbolic matrices which are the Canny-Emiris matrix and its principal minor, as before. Let's show its use in the system corresponding to Example 4.1.","category":"page"},{"location":"code/2_resultants/","page":"Resultants","title":"Resultants","text":"julia> N = [2]\n1-element Vector{Int64}:\n 2\n\njulia> D = [2 2 1]\n1×3 Matrix{Int64}:\n 2  2  1\n\njulia> CE,PM = CannyEmiris.Multihomogeneous(D,N)\nThe rows of the Canny-Emiris matrix x^{b-a(b)}F_{i(b)} are: \n[2, 1]-> x^[2, 1]*F_2\n[3, 1]-> x^[3, 1]*F_2\n[1, 2]-> x^[1, 2]*F_2\n[2, 2]-> x^[2, 2]*F_2\n[1, 3]-> x^[1, 3]*F_2\n[4, 1]-> x^[2, 1]*F_1\n[3, 2]-> x^[1, 2]*F_1\n[2, 3]-> x^[2, 1]*F_0\n[1, 4]-> x^[1, 2]*F_0\n\nThe size of the greedy Canny-Emiris matrix is: 9\nThe degree of the resultant is: 8\n\n(SymPy.Sym[(u_{2, [0, 0]}) (u_{2, [1, 0]}) … 0 0; 0 (u_{2, [0, 0]}) … 0 0; … ; (u_{0, [0, 0]}) (u_{0, [1, 0]}) … (u_{0, [0, 2]}) 0; 0 0 … (u_{0, [1, 1]}) (u_{0, [0, 2]})], SymPy.Sym[(u_{2, \n[0, 0]});;])\n\njulia> CE\n9×9 Matrix{SymPy.Sym}:\n (u_{2, [0, 0]})  (u_{2, [1, 0]})                0  (u_{2, [0, 1]})                0                0                0                0                0\n               0  (u_{2, [0, 0]})                0                0                0  (u_{2, [1, 0]})  (u_{2, [0, 1]})                0                0\n               0                0  (u_{2, [0, 0]})  (u_{2, [1, 0]})  (u_{2, [0, 1]})                0                0                0                0\n               0                0                0  (u_{2, [0, 0]})                0                0  (u_{2, [1, 0]})  (u_{2, [0, 1]})                0\n               0                0                0                0  (u_{2, [0, 0]})                0                0  (u_{2, [1, 0]})  (u_{2, [0, 1]})\n (u_{1, [0, 0]})  (u_{1, [1, 0]})                0  (u_{1, [0, 1]})                0  (u_{1, [2, 0]})  (u_{1, [1, 1]})  (u_{1, [0, 2]})                0\n               0                0  (u_{1, [0, 0]})  (u_{1, [1, 0]})  (u_{1, [0, 1]})                0  (u_{1, [2, 0]})  (u_{1, [1, 1]})  (u_{1, [0, 2]})\n (u_{0, [0, 0]})  (u_{0, [1, 0]})                0  (u_{0, [0, 1]})                0  (u_{0, [2, 0]})  (u_{0, [1, 1]})  (u_{0, [0, 2]})                0\n               0                0  (u_{0, [0, 0]})  (u_{0, [1, 0]})  (u_{0, [0, 1]})                0  (u_{0, [2, 0]})  (u_{0, [1, 1]})  (u_{0, [0, 2]})\n\njulia> PM\n1×1 Matrix{SymPy.Sym}:\n (u_{2, [0, 0]})","category":"page"},{"location":"code/2_resultants/","page":"Resultants","title":"Resultants","text":"Again, this function will return nothing if the dimension conditions and the order. ","category":"page"},{"location":"code/2_resultants/#Other-functions-to-call","page":"Resultants","title":"Other functions to call","text":"","category":"section"},{"location":"code/2_resultants/#MultihomogeneousEmbedding","page":"Resultants","title":"MultihomogeneousEmbedding","text":"","category":"section"},{"location":"code/2_resultants/","page":"Resultants","title":"Resultants","text":"The function CannyEmiris.MultihomogeneousEmbedding(MULTI_A::Matrix{Int64}, MULTI_N::Vector{Int64}))  builds the embedding of the multihomogeneous system into a zonotope system. This system is given by the supports:","category":"page"},{"location":"code/2_resultants/","page":"Resultants","title":"Resultants","text":" overlinemathcalA_i = big(b_jl)_l = 1dotss^j=1dotsn_l in oplus_j=1^smathbbZ^n_j    0 leq sum_j = J^n_lb_jl leq d_ijquad l = 1dotss quad J = 1dotsn_l","category":"page"},{"location":"code/2_resultants/","page":"Resultants","title":"Resultants","text":"\njulia> N = [2]\n1-element Vector{Int64}:\n 2\n\njulia> D = [2 2 1]\n1×3 Matrix{Int64}:\n 2  2  1\n\njulia> CannyEmiris.MultihomogeneousEmbedding(D,N)\n([2 2 1; 2 2 1], [1 -1; 0 1])\n ","category":"page"},{"location":"code/2_resultants/#GenerateTypeFunctions","page":"Resultants","title":"GenerateTypeFunctions","text":"","category":"section"},{"location":"code/2_resultants/","page":"Resultants","title":"Resultants","text":"The function CannyEmiris.GenerateTypeFunctions(n::Int) constructs the iterator that produces all the type functions varphi1dotsn xrightarrow 0dotsn satisfying the condition on the greedy subset mathcalG:","category":"page"},{"location":"code/2_resultants/","page":"Resultants","title":"Resultants","text":"sum_i = 0^I-1 t_bi leq I quad forall I = 1dotsn where t_bi = varphi^-1(i)","category":"page"},{"location":"code/2_resultants/","page":"Resultants","title":"Resultants","text":"julia> p = CannyEmiris.GenerateTypeFunctions(2)\nMain.CannyEmiris.TypeFunctions{Vector{Int64}}([0, 1, 2], 2)\n\njulia> for x in p println(x) end\n[0, 1]\n[0, 2]\n[1, 0]\n[1, 1]\n[1, 2]\n[2, 0]\n[2, 1]\n[2, 2]\n ","category":"page"},{"location":"code/2_resultants/","page":"Resultants","title":"Resultants","text":"AlgebraicSolvers.CannyEmiris","category":"page"},{"location":"code/2_resultants/#AlgebraicSolvers.CannyEmiris","page":"Resultants","title":"AlgebraicSolvers.CannyEmiris","text":"Module which provides optimised resultant constructions for multi-homogeneous systems.\n\n\n\n\n\n","category":"module"},{"location":"code/2_resultants/","page":"Resultants","title":"Resultants","text":"AlgebraicSolvers.CannyEmiris.Multihomogeneous","category":"page"},{"location":"code/2_resultants/#AlgebraicSolvers.CannyEmiris.Multihomogeneous","page":"Resultants","title":"AlgebraicSolvers.CannyEmiris.Multihomogeneous","text":"Add some documentation: what are A, H, Coeff, ...\n\n\n\n\n\n","category":"function"},{"location":"code/2_resultants/","page":"Resultants","title":"Resultants","text":"AlgebraicSolvers.CannyEmiris.Zonotopes \nAlgebraicSolvers.CannyEmiris.NumCoeff","category":"page"},{"location":"code/2_resultants/#AlgebraicSolvers.CannyEmiris.Zonotopes","page":"Resultants","title":"AlgebraicSolvers.CannyEmiris.Zonotopes","text":"I will add some documentation: what are A, H, Coeff, ...\n\n\n\n\n\n","category":"function"},{"location":"code/2_resultants/#AlgebraicSolvers.CannyEmiris.NumCoeff","page":"Resultants","title":"AlgebraicSolvers.CannyEmiris.NumCoeff","text":"Coeff = NumCoeff(F,X)\n\nF is an array of polynomials\nX is the variable vector\n\nThe output Coeff: (i,E)-> Coeff(i,E) is a function, which returns the coefficient of the monomial X^E in the polynomial F[i].\n\n\n\n\n\n","category":"function"},{"location":"code/1_solvers/#Solvers","page":"Solvers","title":"Solvers","text":"","category":"section"},{"location":"code/1_solvers/","page":"Solvers","title":"Solvers","text":"Several types of solvers are available in the packages.  They differ in the way the truncated normal form is computed.","category":"page"},{"location":"code/1_solvers/#Solvers-using-resultant-constructions","page":"Solvers","title":"Solvers using resultant constructions","text":"","category":"section"},{"location":"code/1_solvers/","page":"Solvers","title":"Solvers","text":"AlgebraicSolvers.solve_macaulay\nAlgebraicSolvers.solve_toric","category":"page"},{"location":"code/1_solvers/#AlgebraicSolvers.solve_macaulay","page":"Solvers","title":"AlgebraicSolvers.solve_macaulay","text":"Xi = solve_macaulay(P, rho ; verbose = false)\n\nP polynomial system\nrho (optional) degree of regularity for the Sylvester matrix construction (optional)\n\nSolve the system P=[p1, ..., pn], building Sylvester matrix of all monomial multiples mi*pi in degree ≤ ρ.\n\nThe default value for ρ is ∑ deg(pi) - n + 1.\n\nExample\n\nusing AlgebraicSolvers, DynamicPolynomials\n\nX = @polyvar x y\n\nP = [2-x*y+x^2,y^2+x-2]\n\nXi = solve_macaulay(P)\n\n\n\n\n\n\n","category":"function"},{"location":"code/1_solvers/#AlgebraicSolvers.solve_toric","page":"Solvers","title":"AlgebraicSolvers.solve_toric","text":"solve_toric(P; verbose = false)\n\nP polynomial system\nX array of variables\n\nSolve the system P=[p1, ..., pn], building Sylvester matrix of all monomial multiples mi*pi for mi in supp(∏_{j != i} pj).\n\nExample\n\nusing AlgebraicSolvers, DynamicPolynomials\n\nX = @polyvar x y\n\nP = [y - x*y + x^2,  1 + y + x + x^2]\n\nXi = solve_toric(P)\n\n\n\n\n\n\n","category":"function"},{"location":"code/1_solvers/#Solvers-using-Groebner-basis-computation","page":"Solvers","title":"Solvers using Groebner basis computation","text":"","category":"section"},{"location":"code/1_solvers/","page":"Solvers","title":"Solvers","text":"AlgebraicSolvers.solve_groebner","category":"page"},{"location":"code/1_solvers/#AlgebraicSolvers.solve_groebner","page":"Solvers","title":"AlgebraicSolvers.solve_groebner","text":"Xi, G, B = solve_groebner(P::Vector; verbose = false)\n\nSolve the system of polynomials P. It outputs:\n\nXi the complex solution points, one per column of Xi\nG the computed Groebner basis\nB the basis of the quotient by the ideal of the equations\n\nExample:\n\nusing AbstractAlgebra\n\nR, (x,y) = QQ[\"x\",\"y\"]\n\nP = [x^2-y, x^2*y-4]\n\nXi, G, B = solve_groebner(P)\n\n\n\n\n\n","category":"function"},{"location":"code/1_solvers/#Resultant-matrices","page":"Solvers","title":"Resultant matrices","text":"","category":"section"},{"location":"code/1_solvers/","page":"Solvers","title":"Solvers","text":"Different constructions of resultant matrices are available, including projective or Macaulay resultant matrices, toric resultant matrices.","category":"page"},{"location":"code/1_solvers/","page":"Solvers","title":"Solvers","text":"AlgebraicSolvers.matrix_macaulay\nAlgebraicSolvers.matrix_toric","category":"page"},{"location":"code/1_solvers/#AlgebraicSolvers.matrix_macaulay","page":"Solvers","title":"AlgebraicSolvers.matrix_macaulay","text":"R, L = matrix_macaulay(P, X, rho, ish = false)\n\nP polynomial system\nX (optional) array of variables\nrho (optional) maximal degree of the multiples of P\nish (optional) set to true if the polynomials are homogeneous\n\nIt outputs \n\nR the transpose of Sylvester matrix of all monomial multiples mi*pi in degree ≤ rho.\nL array of monomials indexing the columns of R\n\n\n\n\n\n","category":"function"},{"location":"code/1_solvers/#AlgebraicSolvers.matrix_toric","page":"Solvers","title":"AlgebraicSolvers.matrix_toric","text":"R, L = matrix_toric(P, A)\n\nwhere\n\nP polynomial system\nA array of supports of Pi\n\nIt outputs \n\nR transpose of the Sylvester matrix of all monomial multiples mi*pi for mi in supp(∏_{j != i} pj).\nL the list of monomials indexing the colums of R\n\n\n\n\n\n","category":"function"},{"location":"code/1_solvers/#Truncated-Normal-Forms","page":"Solvers","title":"Truncated  Normal Forms","text":"","category":"section"},{"location":"code/1_solvers/","page":"Solvers","title":"Solvers","text":"AlgebraicSolvers.tnf_macaulay\nAlgebraicSolvers.tnf_toric","category":"page"},{"location":"code/1_solvers/#AlgebraicSolvers.tnf_macaulay","page":"Solvers","title":"AlgebraicSolvers.tnf_macaulay","text":"N, L = tnf_macaulay(P, rho)\n\nCompute the Truncated Normal Form of P=[p1, ..., pn], using Macaulay matrix of all monomial multiples mi*pi in degree ≤ ρ.\n\nThe default value for ρ is ∑ deg(pi) - n + 1.\n\n\n\n\n\n","category":"function"},{"location":"code/1_solvers/#AlgebraicSolvers.tnf_toric","page":"Solvers","title":"AlgebraicSolvers.tnf_toric","text":"N, L = tnf_toric(P, A = support.(P))\n\nCompute the Truncated Normal Form of P=[p1, ..., pn], using toric resultant matrix of all monomial multiples mi*pi in degree ≤ ρ.\n\nThe default value for ρ is ∑ deg(pi) - n + 1.\n\n\n\n\n\n","category":"function"},{"location":"code/1_solvers/#Linear-algebra-in-quotient","page":"Solvers","title":"Linear algebra in quotient","text":"","category":"section"},{"location":"code/1_solvers/","page":"Solvers","title":"Solvers","text":"AlgebraicSolvers.matrix_mult\nAlgebraicSolvers.quotient_basis","category":"page"},{"location":"code/1_solvers/#AlgebraicSolvers.matrix_mult","page":"Solvers","title":"AlgebraicSolvers.matrix_mult","text":"M = matrix_mult(p, g::AbstractVector, Idx::Dict)\n\nCompute the matrix of multication by p modulo g in the basis associated to the basis dictionary Idx. It is assumed that g is a Groebner basis and that the quotient is finite dimensional.\n\n\n\n\n\nM = matrix_mult(p, G::AbstractVector, B::AbstractVector)\n\nCompute the matrix of multication by p modulo G in the basis B. It is assumed that G is a Groebner basis and that the quotient is finite dimensional.\n\n\n\n\n\n","category":"function"},{"location":"code/1_solvers/#AlgebraicSolvers.quotient_basis","page":"Solvers","title":"AlgebraicSolvers.quotient_basis","text":"B, Idx = quotient_basis(g::AbstractVector)\n\nCompute the monomial basis of the quotient by g, as the complementary of the initial ideal of g, assuming that g is a Groebner basis and it is finite. It outputs: \n\nB the basis of monomlials  in increasing monomial order.\nIdx the dictionary associting to a monomial its index in the basis.   \n\n\n\n\n\n","category":"function"},{"location":"expl/3.Groebner/#Solving-using-Groebner-basis-computation","page":"Solving using Groebner basis computation","title":"Solving using Groebner basis computation","text":"","category":"section"},{"location":"expl/3.Groebner/","page":"Solving using Groebner basis computation","title":"Solving using Groebner basis computation","text":"We load the needed packages, define the ring of polynomials we will need, and the polynomial system we will solve:","category":"page"},{"location":"expl/3.Groebner/","page":"Solving using Groebner basis computation","title":"Solving using Groebner basis computation","text":"using AbstractAlgebra, Groebner, AlgebraicSolvers\n\nR, (x,y,z) = QQ[\"x\",\"y\",\"z\"]\n\n\nP = [\n    8*x^2*y^2+5*x*y^3+3*x^3*z+x^2*y*z,\n    x^5+2*y^2*z^2+13*y^2*z^3+5*y*z^4,\n    8*x^3+12*y^3+x*z^2+3,\n    7*x^2*y^4+18*x*y^3*z^2+y^3*z^3\n    ]","category":"page"},{"location":"expl/3.Groebner/","page":"Solving using Groebner basis computation","title":"Solving using Groebner basis computation","text":"4-element Vector{AbstractAlgebra.Generic.MPoly{Rational{BigInt}}}:\n 3*x^3*z + 8*x^2*y^2 + x^2*y*z + 5*x*y^3\n x^5 + 13*y^2*z^3 + 2*y^2*z^2 + 5*y*z^4\n 8*x^3 + x*z^2 + 12*y^3 + 3\n 7*x^2*y^4 + 18*x*y^3*z^2 + y^3*z^3","category":"page"},{"location":"expl/3.Groebner/","page":"Solving using Groebner basis computation","title":"Solving using Groebner basis computation","text":"We solve the system P using the Groebner solver:","category":"page"},{"location":"expl/3.Groebner/","page":"Solving using Groebner basis computation","title":"Solving using Groebner basis computation","text":"Xi, G, B = solve_groebner(P);","category":"page"},{"location":"expl/3.Groebner/","page":"Solving using Groebner basis computation","title":"Solving using Groebner basis computation","text":"and we get 6 solutions (the columns of Xi):","category":"page"},{"location":"expl/3.Groebner/","page":"Solving using Groebner basis computation","title":"Solving using Groebner basis computation","text":"Xi","category":"page"},{"location":"expl/3.Groebner/","page":"Solving using Groebner basis computation","title":"Solving using Groebner basis computation","text":"3×6 Matrix{ComplexF64}:\n       0.0+0.0im          …           0.0+0.0im\n -0.629961-1.664e-31im            0.31498+0.545562im\n       0.0+4.22682e-10im     -1.18437e-10-1.9944e-10im","category":"page"},{"location":"expl/3.Groebner/#How-does-it-work-?","page":"Solving using Groebner basis computation","title":"How does it work ?","text":"","category":"section"},{"location":"expl/3.Groebner/","page":"Solving using Groebner basis computation","title":"Solving using Groebner basis computation","text":"First we compute the Groebner basis of P (for the degree reverse lexicographic ordering):","category":"page"},{"location":"expl/3.Groebner/","page":"Solving using Groebner basis computation","title":"Solving using Groebner basis computation","text":"G = groebner(P)","category":"page"},{"location":"expl/3.Groebner/","page":"Solving using Groebner basis computation","title":"Solving using Groebner basis computation","text":"3-element Vector{AbstractAlgebra.Generic.MPoly{Rational{BigInt}}}:\n z^2\n y^3 + 1//4\n x","category":"page"},{"location":"expl/3.Groebner/","page":"Solving using Groebner basis computation","title":"Solving using Groebner basis computation","text":"Then we deduce the basis of quotient by the ideal (P):","category":"page"},{"location":"expl/3.Groebner/","page":"Solving using Groebner basis computation","title":"Solving using Groebner basis computation","text":"B, BIdx = quotient_basis(G); B","category":"page"},{"location":"expl/3.Groebner/","page":"Solving using Groebner basis computation","title":"Solving using Groebner basis computation","text":"6-element Vector{AbstractAlgebra.Generic.MPoly{Rational{BigInt}}}:\n 1\n z\n y\n y*z\n y^2\n y^2*z","category":"page"},{"location":"expl/3.Groebner/","page":"Solving using Groebner basis computation","title":"Solving using Groebner basis computation","text":"(Here Bidx is a dictionary of monomials giving their index in the basis B).","category":"page"},{"location":"expl/3.Groebner/","page":"Solving using Groebner basis computation","title":"Solving using Groebner basis computation","text":"Then we compute the matrices of multiplication by the variables in the basis Bof the quotient:","category":"page"},{"location":"expl/3.Groebner/","page":"Solving using Groebner basis computation","title":"Solving using Groebner basis computation","text":"M = [mult_matrix(v, G, B) for v in [x,y,z]]","category":"page"},{"location":"expl/3.Groebner/","page":"Solving using Groebner basis computation","title":"Solving using Groebner basis computation","text":"3-element Vector{Matrix{Rational{BigInt}}}:\n [0//1 0//1 … 0//1 0//1; 0//1 0//1 … 0//1 0//1; … ; 0//1 0//1 … 0//1 0//1; 0//1 0//1 … 0//1 0//1]\n [0//1 0//1 … -1//4 0//1; 0//1 0//1 … 0//1 -1//4; … ; 0//1 0//1 … 0//1 0//1; 0//1 0//1 … 0//1 0//1]\n [0//1 0//1 … 0//1 0//1; 1//1 0//1 … 0//1 0//1; … ; 0//1 0//1 … 0//1 0//1; 0//1 0//1 … 1//1 0//1]","category":"page"},{"location":"expl/3.Groebner/","page":"Solving using Groebner basis computation","title":"Solving using Groebner basis computation","text":"Then we triangularise them jointly in the same basis and deduce the points Xifrom the values on the diagonal of triangularised M_i.","category":"page"},{"location":"expl/3.Groebner/","page":"Solving using Groebner basis computation","title":"Solving using Groebner basis computation","text":"For that purpose, we compute the Schur factorization of a random combination of the matrices.   ","category":"page"},{"location":"expl/3.Groebner/","page":"Solving using Groebner basis computation","title":"Solving using Groebner basis computation","text":"This gives us the points (possibly repeated with their multiplicity):","category":"page"},{"location":"expl/3.Groebner/","page":"Solving using Groebner basis computation","title":"Solving using Groebner basis computation","text":"Xi","category":"page"},{"location":"expl/3.Groebner/","page":"Solving using Groebner basis computation","title":"Solving using Groebner basis computation","text":"3×6 Matrix{ComplexF64}:\n       0.0+0.0im          …           0.0+0.0im\n -0.629961-1.664e-31im            0.31498+0.545562im\n       0.0+4.22682e-10im     -1.18437e-10-1.9944e-10im","category":"page"},{"location":"#AlgebraicSolvers","page":"Home","title":"AlgebraicSolvers","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Package for the solution of polynomial systems, using eigen computation. It outputs all the complex solutions if the system is zero-dimensional.","category":"page"},{"location":"","page":"Home","title":"Home","text":"It proceeds as follows:","category":"page"},{"location":"","page":"Home","title":"Home","text":"Compute a Truncated Normal Form (using either resultant constructions or Groebner basis computation)\nCompute operators of multiplication by the variables in a basis\nCompute a Joint triangularization of the multplication matrices to obtain the (complex) solutions, with their multiplicity.","category":"page"},{"location":"","page":"Home","title":"Home","text":"Pages = map(file -> joinpath(\"code\", file), filter(x ->endswith(x, \"md\"), readdir(\"code\"))) ","category":"page"},{"location":"#Dependencies","page":"Home","title":"Dependencies","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"For the solver based on Groebner basis computation, we are use the package Grobner.jl.","category":"page"}]
}
