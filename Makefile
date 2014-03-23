all: doc.pdf

.PHONY: all clean doc.pdf

doc.pdf: doc.tex rest.pgf
	latexmk -pdf $<

rest.pgf: %.pgf: %.py
	python $<

20points.pgf:
	python poly.py $@

clean:
	latexmk -c
