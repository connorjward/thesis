# Compilers
TEX = latexmk -aux-directory=build -lualatex="lualatex -shell-escape" -lualatex

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

# Remove temporary files
.PHONY: tidy
tidy:
	-rm -rf ./build/

# Remove PDFs
.PHONY: clean
clean: tidy
	-rm ./*.pdf ./chapters/*.pdf ./figures/*.pdf 2>/dev/null || true

# Remove specific PDF
%.clean: ALWAYS
	-rm ./$*.pdf 2>/dev/null || true

$(THESIS): %.pdf: $(CHAPTERS) %.clean %.tex
	$(TEX) $*

$(CHAPTERS): %.pdf: $(FIGURES) %.tex
	$(TEX) -output-directory=chapters $*

$(FIGURES): %.pdf: %.tex
	$(TEX) -output-directory=figures $*

# Use to always run a command even if the files are unchanged
.PHONY: ALWAYS
ALWAYS:
