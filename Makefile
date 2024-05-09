# Compilers
# TEX = xelatex -shell-escape
TEX = lualatex -shell-escape
BIB = biber
# PY = pythontex
GLOS = makeglossaries
PYTHON = python3

# Documents
THESIS = thesis.pdf

CHAPTERS_DIR = chapters
CHAPTERS_SRC = $(wildcard $(CHAPTERS_DIR)/*.tex)
CHAPTERS = $(patsubst $(CHAPTERS_DIR)/%.tex,$(CHAPTERS_DIR)/%.pdf,$(CHAPTERS_SRC))

FIGURES_DIR = figures
FIGURES_SRC = $(wildcard $(FIGURES_DIR)/*.tex)
FIGURES = $(patsubst $(FIGURES_DIR)/%.tex,$(FIGURES_DIR)/%.pdf,$(FIGURES_SRC))

# Make everything
.PHONY: all
all: $(THESIS)

# Make all individual chapters
.PHONY: chapters
chapters: ALWAYS $(CHAPTERS)
	@echo $(CHAPTERS)

# Make figures
.PHONY: figures
figures: $(FIGURES)

# Remove build artefacts
# .PHONY: tidy
# tidy:
# 	-rm ./*.aux ./**/*.aux \
# 	./*.bbl ./**/*.bbl \
# 	./*.bcf ./**/*.bcf \
# 	./*.blg ./**/*.blg \
# 	./*.log ./**/*.log \
# 	./*.pyg ./**/*.pyg \
# 	./*.run.xml ./**/*.run.xml \
# 	./*.toc 2>/dev/null || true

# Remove PDFs
.PHONY: clean
# clean: tidy
clean:
	-rm ./*.pdf 2>/dev/null || true

# Remove specific PDF
# %.clean: ALWAYS
# 	-rm ./$*.pdf 2>/dev/null || true

# $(THESIS): %.pdf: $(CHAPTERS) %.clean %.tex
# 	$(TEX) $*
# 	-$(BIB) $*
# 	# $(GLOS) $*  # TODO
# 	$(TEX) $*
# 	# -$(GLOS) $*  # TODO
# 	$(TEX) $*

$(CHAPTERS): %.pdf: $(FIGURES) %.tex
	latexmk -aux-directory=build -output-directory=chapters -lualatex="lualatex -shell-escape" -lualatex $*

$(FIGURES): %.pdf: %.tex
	latexmk -aux-directory=build -output-directory=figures -lualatex="lualatex -shell-escape" -lualatex $*
