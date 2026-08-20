LATEX  = lualatex
FLAGS  = -interaction=nonstopmode -halt-on-error
MAIN   = paper
JOB    = cxor-paper

# ------------------------------------------------------------------ archive
DIST     = dist
HAL_NAME = $(JOB)-source
HAL_DIR  = $(DIST)/$(HAL_NAME)
HAL_ZIP  = $(DIST)/$(HAL_NAME).zip

# Directories that never belong in a source archive.
HAL_PRUNE = -name .git -o -name .github -o -name $(DIST) -o -name build \
            -o -name .vscode -o -name .idea -o -name __pycache__
# Extensions LaTeX may need in order to compile the manuscript.  The .bbl is
# matched here too, so archival systems get the bibliography without biber.
HAL_KEEP  = -name '*.tex' -o -name '*.bib' -o -name '*.bbl' -o -name '*.sty' \
            -o -name '*.cls' -o -name '*.pdf' -o -name '*.png' -o -name '*.jpg' \
            -o -name '*.jpeg' -o -name '*.eps' -o -name '*.svg'
# Generated files that happen to share one of those extensions.
HAL_DROP  = ^$(JOB)\.pdf$$|^$(JOB)-blx\.bib$$
# Extensionless files worth shipping, if present.
HAL_ALSO  = $(wildcard PAPER-LICENSE)

# Real files, derived from the project layout rather than a hand-written list.
HAL_FILES = { find . \( $(HAL_PRUNE) \) -prune -o -type f \( $(HAL_KEEP) \) -print \
                | sed 's|^\./||' | grep -Ev '$(HAL_DROP)'; \
              for f in $(HAL_ALSO); do echo "$$f"; done; } | LC_ALL=C sort -u
# biblatex resolves the bibliography as <jobname>.bbl, and a plain
# "lualatex $(MAIN).tex" uses jobname $(MAIN).  Ship a second copy under that
# name so the archive also compiles where biber is not run.
HAL_ALIAS = if [ -f $(JOB).bbl ] && [ "$(JOB)" != "$(MAIN)" ]; then echo $(MAIN).bbl; fi
HAL_LIST  = { $(HAL_FILES); $(HAL_ALIAS); } | LC_ALL=C sort -u

all: $(JOB).pdf

$(JOB).pdf: $(MAIN).tex preamble.tex paper.bib $(wildcard paper/*.tex) $(wildcard fig/*.tex)
	$(LATEX) $(FLAGS) -jobname=$(JOB) $(MAIN).tex
	biber $(JOB)
	$(LATEX) $(FLAGS) -jobname=$(JOB) $(MAIN).tex
	$(LATEX) $(FLAGS) -jobname=$(JOB) $(MAIN).tex

check:
	python3 cta_spec.py

# Print exactly what hal-source would archive, without creating it.
hal-source-list:
	@$(HAL_LIST)

# Self-contained source archive for a HAL preprint deposit.  Depends on the
# PDF so that the shipped .bbl is current.
hal-source: $(JOB).pdf
	@rm -rf $(HAL_DIR) $(HAL_ZIP)
	@mkdir -p $(HAL_DIR)
	@$(HAL_FILES) | while read -r f; do \
	    mkdir -p "$(HAL_DIR)/$$(dirname "$$f")"; \
	    cp "$$f" "$(HAL_DIR)/$$f"; \
	  done
	@$(HAL_ALIAS) | while read -r f; do cp $(JOB).bbl "$(HAL_DIR)/$$f"; done
	@cd $(DIST) && zip -q -r -X $(HAL_NAME).zip $(HAL_NAME)
	@rm -rf $(HAL_DIR)
	@echo "Created $(HAL_ZIP)"

clean:
	rm -f $(JOB).aux $(JOB).bbl $(JOB).bcf $(JOB).blg $(JOB).log $(JOB).out \
	      $(JOB).run.xml $(JOB).toc $(JOB).fls $(JOB).fdb_latexmk $(JOB)-blx.bib

distclean: clean
	rm -rf $(DIST)
	rm -f $(JOB).pdf

.PHONY: all check clean distclean hal-source hal-source-list
