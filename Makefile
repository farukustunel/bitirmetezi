# Minimal makefile for Sphinx documentation
#

# You can set these variables from the command line.
SPHINXOPTS    =
SPHINXBUILD   = sphinx-build
SOURCEDIR     = .
BUILDDIR      = _build

# Put it first so that "make" without argument is like "make help".
help:
	@$(SPHINXBUILD) -M help "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

.PHONY: help Makefile

poster/%.md: src/%.rst
	pandoc $^ --atx --wrap=none -o $@
	
poster/poster.html: poster/poster.yaml poster/Introduction.md poster/Results.md
	pandoc --section-divs -t html4 --template=poster-template.html --css poster-local.css $^ -o $@


word/%.md: src/%.rst
	pandoc $^ --atx --wrap=none -o $@
	
report:
	cd src; pandoc Introduction.rst MaterialsandMethods.rst Results.rst --reference-doc=../template.docx -o report.docx --toc -M title="Assembly of NGS reads to determine DNA sequence of Multidrug Resistant Plasmid isolated in Cow and Human E.coli and *C.freundii*" -M autor="Faruk Üstünel"

# Catch-all target: route all unknown targets to Sphinx using the new
# "make mode" option.  $(O) is meant as a shortcut for $(SPHINXOPTS).
#%: Makefile
#	@$(SPHINXBUILD) -M $@ "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)
	
	
