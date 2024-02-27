########################################
# Makefile for thesis                  #
########################################

# Compilers
TEX = xelatex -shell-escape
BIB = bibtex
# PY = pythontex
GLOS = makeglossaries
PYTHON = python3

# Documents
FULL = thesis.pdf
CHAPTERS =  \
	chapters/1_introduction.pdf \
	chapters/2_background.pdf \
	chapters/3_mesh_data_layouts.pdf \
	chapters/4_dsl.pdf \
	chapters/5_parallel.pdf \
	chapters/6_rest.pdf \

# Figures
FIG_DIR = figures
FIG_TEX = $(wildcard $(FIG_DIR)/*.tex)
FIGURES = $(patsubst $(FIG_DIR)/%.tex,$(FIG_DIR)/%.pdf,$(FIG_TEX))

# Make everything
all: $(CHAPTERS) $(FULL)
	sed -i 's/{Contents}{v}/{Contents}{vii}/' thesis.toc
	$(TEX) thesis
	make tidy

# Make full document (DOES NOT CHECK FOR CHAPTERS)
full: $(FULL)

# Make all individual chapters
chapters: $(CHAPTERS)

# LaTeX file recipe
# %.pdf: $(FIGURES) %.clean %.tex
%.pdf: %.clean %.tex
	$(TEX) $*
	-$(BIB) $*
	# -$(PY) $*
	# $(GLOS) $*
	$(TEX) $*
	#-$(GLOS) $*
	$(TEX) $*

# Specific recipe for pdf figures
# $(FIGURES): %.pdf: %.tex
# 	$(PYTHON) figures/pstricks2pdf.py -o $@ $<

# Figures recipe
.PHONY: figures
figures: $(FIGURES)
	#@echo FIG_DIR = $(FIG_DIR)
	#@echo FIGTEX = $(FIG_TEX)
	#@echo FIGURES = $(FIGURES)

# Remove all auxilliary files
.PHONY: tidy
tidy:
	-rm ./*.aux ./*.ac? ./*.alg ./*.gl? ./*ist ./*.toc ./*.lo? ./*.out ./*.pytxcode ./*.bbl ./*.blg ./.*.swp 2>/dev/null || true
	-rm -r pythontex-files-* 2>/dev/null || true
	
# Remove pdfs
.PHONY: clean
clean: tidy
	-rm ./*.pdf 2>/dev/null || true

# Remove specific pdf
%.clean: ALWAYS
	-rm ./$*.pdf 2>/dev/null || true

# Target always gets executed
ALWAYS:
