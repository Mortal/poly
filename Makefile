all: doc.pdf

.PHONY: all clean doc.pdf

doc.pdf: doc.tex rest.pgf
	latexmk -pdf $<

rest.pgf: %.pgf: %.py
	$(PYTHON) $<

20points.pgf:
	$(PYTHON) poly.py $@

clean:
	latexmk -c
