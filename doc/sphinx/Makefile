# Makefile for Sphinx documentation
#

# You can set these variables from the command line.
PYTHON        = python
SVNROOT       = http://svn.python.org/projects
SPHINXOPTS    =
SPHINXBUILD   = $(PYTHON) tools/sphinx-build.py
PAPER         =

# Internal variables.
PAPEROPT_a4     = -D latex_paper_size=a4
PAPEROPT_letter = -D latex_paper_size=letter
ALLSPHINXOPTS   = -d build/doctrees $(PAPEROPT_$(PAPER)) $(SPHINXOPTS) source

.PHONY: help checkout update clean html htmlhelp latex alldoc changes linkcheck

help:
	@echo "Please use \`make <target>' where <target> is one of"
	@echo "  html      to make standalone HTML files"
	@echo "  htmlhelp  to make HTML files and a HTML help project"
	@echo "  latex     to make LaTeX files, you can set PAPER=a4 or PAPER=letter"
	@echo "  alldoc    to make pdfLaTeX and HTML files, you can set PAPER=a4 or PAPER=letter"
	@echo "  changes   to make an overview over all changed/added/deprecated items"
	@echo "  linkcheck to check all external links for integrity"

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

clean:
	-rm -rf build/
	-rm -rf html/
	-rm -rf tools/

html: update
	mkdir -p build/html build/doctrees
	$(SPHINXBUILD) -b html $(ALLSPHINXOPTS) build/html
	@echo
	@echo "Build finished. The HTML pages are in build/html."

htmldoc: html
	rm -rf html/
	mv -f build/html ./

htmlhelp:
	mkdir -p build/htmlhelp build/doctrees
	$(SPHINXBUILD) -b htmlhelp $(ALLSPHINXOPTS) build/htmlhelp
	@echo
	@echo "Build finished; now you can run HTML Help Workshop with the" \
	      ".hhp project file in build/htmlhelp."

latex:
	mkdir -p build/latex build/doctrees
	$(SPHINXBUILD) -b latex $(ALLSPHINXOPTS) build/latex
	@echo
	@echo "Build finished; the LaTeX files are in build/latex."
	@echo "Run \`make all-pdf' or \`make all-ps' in that directory to" \
	      "run these through (pdf)latex."

alldoc: htmldoc latex
	cd build/latex && make all-pdf
	mv build/latex/PythonGetfemInterface.pdf html/

changes:
	mkdir -p build/changes build/doctrees
	$(SPHINXBUILD) -b changes $(ALLSPHINXOPTS) build/changes
	@echo
	@echo "The overview file is in build/changes."

linkcheck:
	mkdir -p build/linkcheck build/doctrees
	$(SPHINXBUILD) -b linkcheck $(ALLSPHINXOPTS) build/linkcheck
	@echo
	@echo "Link check complete; look for any errors in the above output " \
	      "or in build/linkcheck/output.txt."
