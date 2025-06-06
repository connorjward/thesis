# Compilers
TEX = latexmk -lualatex="lualatex -shell-escape" -lualatex

# Documents
THESIS = thesis.pdf

CHAPTERS_DIR = chapters
CHAPTERS_SRC = $(wildcard $(CHAPTERS_DIR)/*.tex)
CHAPTERS = $(patsubst $(CHAPTERS_DIR)/%.tex,%.pdf,$(CHAPTERS_SRC))

FIGURES_DIR = figures
FIGURES_SRC = $(wildcard $(FIGURES_DIR)/*.tex)
FIGURES = $(patsubst $(FIGURES_DIR)/%.tex,%.pdf,$(FIGURES_SRC))

# Make everything
.PHONY: all
all: $(THESIS)

# Make all individual chapters
.PHONY: chapters
chapters: ALWAYS $(CHAPTERS)
	@echo $(CHAPTERS)

# Make figures
.PHONY: figures
figures: $(FIGURES) roofline.pgf

# Remove temporary files
# .PHONY: tidy
# tidy:
# 	-rm -rf ./build/

# Remove PDFs
.PHONY: clean
# clean: tidy
clean:
	-rm ./*.pdf 2>/dev/null || true

# Remove specific PDF
# %.clean: ALWAYS
# 	-rm ./$*.pdf 2>/dev/null || true

$(THESIS): %.pdf: $(CHAPTERS) %.tex titlepage.tex abstract.tex declaration.tex acknowledgements.tex dedication.tex
	$(TEX) $*

$(CHAPTERS): %.pdf: $(FIGURES) $(CHAPTERS_DIR)/%.tex
	$(TEX) $(CHAPTERS_DIR)/$*

$(FIGURES): %.pdf: $(FIGURES_DIR)/%.tex
	$(TEX) $(FIGURES_DIR)/$*

roofline.pgf:
	python scripts/plot_roofline.py -o roofline.pgf

# Use to always run a command even if the files are unchanged
.PHONY: ALWAYS
ALWAYS:
