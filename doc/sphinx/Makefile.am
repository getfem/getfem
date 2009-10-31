# Makefile for Sphinx documentation
#

# You can set these variables from the command line.
PYTHON        = python
SVNROOT       = http://svn.python.org/projects
SPHINXOPTS    =
SPHINXBUILD   = $(PYTHON) tools/sphinx-build.py
PAPER         =
SOURCES       = $(srcdir)/source
DISTVERSION   = @VERSION@


# Internal variables.
PAPEROPT_a4     = -D latex_paper_size=a4
PAPEROPT_letter = -D latex_paper_size=letter
ALLSPHINXOPTS   = -b $(BUILDER) -d build/doctrees $(PAPEROPT_$(PAPER)) \
                  $(SPHINXOPTS) $(SOURCES) build/$(BUILDER)

.PHONY: help checkout update build html htmlhelp clean latex linkcheck dist

help:
	@echo "Please use \`make <target>' where <target> is one of"
	@echo "  view      to view standalone HTML files"
	@echo "  html      to make standalone HTML files"
	@echo "  htmlhelp  to make HTML files and a HTML help project"
	@echo "  htmlview  to make and view standalone HTML files"
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
	svn update tools
	svn update tools/sphinx
	svn update tools/docutils
	svn update tools/jinja2
	svn update tools/pygments

images:
	-cd source/userdoc/images/; make png
	-cd source/project/images/; make png


build: checkout images
	mkdir -p build/$(BUILDER) build/doctrees
	$(SPHINXBUILD) $(ALLSPHINXOPTS)
	@echo

view:
	$(PYTHON) -c "import webbrowser; webbrowser.open('build/html/index.html')"

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

dist: DIST = build/dist
dist: DESTINATION = build/html/dist/
dist:
	-rm -rf $(DIST)
	-rm -rf $(DESTINATION)
	mkdir -p $(DIST)
	mkdir -p $(DESTINATION)

	make html # archive the HTML
	cp -pPR build/html $(DIST)/getfem-$(DISTVERSION)-docs-html
	tar -C $(DIST) -cf $(DIST)/getfem-$(DISTVERSION)-docs-html.tar getfem-$(DISTVERSION)-docs-html
	gzip -9 $(DIST)/getfem-$(DISTVERSION)-docs-html.tar
	rm -rf $(DIST)/getfem-$(DISTVERSION)-docs-html

	rm -rf build/latex # archive the letter latex
	make latex PAPER=letter
	(cd build/latex; make clean && make all-pdf && make FMT=pdf tar)
	gzip -9 build/latex/docs-pdf.tar
	cp build/latex/docs-pdf.tar.gz $(DIST)/getfem-$(DISTVERSION)-docs-pdf-letter.tar.gz

	mv $(DIST)/* $(DESTINATION) # move to destination

htmlupload: html
	-mv build/html sphinxdoc
	../../bin/upload_documentation sphinxdoc
