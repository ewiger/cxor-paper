LATEX  = lualatex
FLAGS  = -interaction=nonstopmode -halt-on-error
JOB    = cxor-paper

all: $(JOB).pdf

$(JOB).pdf: paper.tex preamble.tex paper.bib $(wildcard paper/*.tex) $(wildcard fig/*.tex)
	$(LATEX) $(FLAGS) -jobname=$(JOB) paper.tex
	biber $(JOB)
	$(LATEX) $(FLAGS) -jobname=$(JOB) paper.tex
	$(LATEX) $(FLAGS) -jobname=$(JOB) paper.tex

check:
	python3 cta_spec.py

clean:
	rm -f $(JOB).aux $(JOB).bbl $(JOB).bcf $(JOB).blg $(JOB).log $(JOB).out \
	      $(JOB).run.xml $(JOB).toc $(JOB).fls $(JOB).fdb_latexmk $(JOB)-blx.bib

distclean: clean
	rm -f $(JOB).pdf

.PHONY: all check clean distclean
