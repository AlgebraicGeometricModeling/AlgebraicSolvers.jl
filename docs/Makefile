.SUFFIXES: .ipynb .md

MD := $(patsubst %.ipynb,%.md,$(wildcard mrkd/expl/*.ipynb))

all: $(MD) AlgebraicSolvers.jl
	julia --color=yes make.jl

.ipynb.md:
	jupyter nbconvert --to markdown $*.ipynb
	perl -pi -e "s/\\$$\\$$/\\$$/g" $*.md


AlgebraicSolvers.jl:
	git clone -b gh-pages https://github.com/AlgebraicGeometricModeling/AlgebraicSolvers.jl.git

publish:
	cd AlgebraicSolvers.jl; git add docs; git commit -am "update gh-pages"; git push

clean:
	rm ./mrkd/expl/*.md

