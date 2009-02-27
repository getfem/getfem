# Makefile for Sphinx documentation
#

# You can set these variables from the command line.
PYTHON        = python
SVNROOT       = http://svn.python.org/projects
SPHINXOPTS    =
SPHINXBUILD   = $(PYTHON) tools/sphinx-build.py
PAPER         =
SOURCES       = source
DISTVERSION   = 3.1

# Internal variables.
PAPEROPT_a4     = -D latex_paper_size=a4
PAPEROPT_letter = -D latex_paper_size=letter
ALLSPHINXOPTS   = -b $(BUILDER) -d build/doctrees $(PAPEROPT_$(PAPER)) \
                  $(SPHINXOPTS) $(SOURCES) build/$(BUILDER)

.PHONY: help checkout update build html htmlhelp clean latex linkcheck dist

help:
	@echo "Please use \`make <target>' where <target> is one of"
	@echo "  html      to make standalone HTML files"
	@echo "  htmlhelp  to make HTML files and a HTML help project"
	@echo "  latex     to make LaTeX files, you can set PAPER=a4 or PAPER=letter"
	@echo "  alldoc    to make pdfLaTeX and HTML files, you can set PAPER=a4 or PAPER=letter"
	@echo "  changes   to make an overview over all changed/added/deprecated items"
	@echo "  linkcheck to check all external links for integrity"
	@echo "  dist      to create a \"dist\" directory with archived docs for download"

checkout:
	@if [ ! -d tools ]; then \
	  echo "Checking out tools..."; \
	  svn checkout $(SVNROOT)/python/trunk/Doc/tools tools; \
	fi
	@if [ ! -d tools/sphinx ]; then \
	  echo "Checking out Sphinx..."; \
	  svn checkout $(SVNROOT)/doctools/trunk/sphinx tools/sphinx; \
	fi
	@if [ ! -d tools/docutils ]; then \
	  echo "Checking out Docutils..."; \
	  svn checkout $(SVNROOT)/external/docutils-0.5/docutils tools/docutils; \
	fi
	@if [ ! -d tools/jinja2 ]; then \
	  echo "Checking out Jinja..."; \
	  svn checkout $(SVNROOT)/external/Jinja-2.1.1/jinja2 tools/jinja2; \
	fi
	@if [ ! -d tools/pygments ]; then \
	  echo "Checking out Pygments..."; \
	  svn checkout $(SVNROOT)/external/Pygments-0.11.1/pygments tools/pygments; \
	fi

update: checkout
	svn update tools/sphinx
	svn update tools/docutils
	svn update tools/jinja2
	svn update tools/pygments

build: checkout
	mkdir -p build/$(BUILDER) build/doctrees
	$(SPHINXBUILD) $(ALLSPHINXOPTS)
	@echo

html: BUILDER = html
html: build
	@echo "Build finished. The HTML pages are in build/html."

htmlhelp: BUILDER = htmlhelp
htmlhelp: build
	@echo "Build finished; now you can run HTML Help Workshop with the" \
	      ".hhp project file in build/htmlhelp."

latex: BUILDER = latex
latex: build
	@echo "Build finished; the LaTeX files are in build/latex."
	@echo "Run \`make all-pdf' or \`make all-ps' in that directory to" \
	      "run these through (pdf)latex."

linkcheck: BUILDER = linkcheck
linkcheck: build
	@echo "Link check complete; look for any errors in the above output " \
	      "or in build/$(BUILDER)/output.txt."

htmlview: html
	$(PYTHON) -c "import webbrowser; webbrowser.open('build/html/index.html')"

clean:
	-rm -rf build/
	-rm -rf html/
	-rm -rf tools/

dist: DIST = build/dist
dist: DESTINATION = build/html/dist/$(DISTVERSION)
dist:
	-rm -rf $(DIST)
	-rm -rf $(DESTINATION)
	mkdir -p $(DIST)
	mkdir -p $(DESTINATION)

	# archive the HTML
	make html
	cp -pPR build/html $(DIST)/getfem-$(DISTVERSION)-docs-html
	tar -C $(DIST) -cf $(DIST)/getfem-$(DISTVERSION)-docs-html.tar getfem-$(DISTVERSION)-docs-html
	bzip2 -9 -k $(DIST)/getfem-$(DISTVERSION)-docs-html.tar
	(cd $(DIST); zip -q -r -9 getfem-$(DISTVERSION)-docs-html.zip getfem-$(DISTVERSION)-docs-html)
	rm -r $(DIST)/getfem-$(DISTVERSION)-docs-html
	rm $(DIST)/getfem-$(DISTVERSION)-docs-html.tar

	# archive the A4 latex
	-rm -r build/latex
	make latex PAPER=a4
	(cd build/latex; make clean && make all-pdf && make FMT=pdf zip bz2)
	cp build/latex/docs-pdf.zip $(DIST)/getfem-$(DISTVERSION)-docs-pdf-a4.zip
	cp build/latex/docs-pdf.tar.bz2 $(DIST)/getfem-$(DISTVERSION)-docs-pdf-a4.tar.bz2

	# archive the letter latex
	rm -r build/latex
	make latex PAPER=letter
	(cd build/latex; make clean && make all-pdf && make FMT=pdf zip bz2)
	cp build/latex/docs-pdf.zip $(DIST)/getfem-$(DISTVERSION)-docs-pdf-letter.zip
	cp build/latex/docs-pdf.tar.bz2 $(DIST)/getfem-$(DISTVERSION)-docs-pdf-letter.tar.bz2

	# move to destination
	mv $(DIST)/* $(DESTINATION)
