all: AlgebraicSolvers.jl
	julia --color=yes make.jl

AlgebraicSolvers.jl:
	git clone -b gh-pages https://github.com/AlgebraicGeometricModeling/AlgebraicSolvers.jl.git

publish:
	cd AlgebraicSolvers.jl; git add docs; git commit -am "update gh-pages"; git push
