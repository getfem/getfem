# Makefile for Sphinx documentation
#

# You can set these variables from the command line.
PYTHON        = python
SPHINXROOT    = http://svn.python.org/projects
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

.PHONY: help checkout update images build view html htmlview htmlhelp latex pdf linkcheck clean upload


help:
	@echo "Please use \`make <target>' where <target> is one of"
	@echo "  view      to view standalone HTML files"
	@echo "  html      to make standalone HTML files"
	@echo "  htmlview  to make and view standalone HTML files"
	@echo "  htmlhelp  to make HTML files and a HTML help project"
	@echo "  latex     to make LaTeX files, you can set PAPER=a4 or PAPER=letter"
	@echo "  pdf       to make pdfLaTeX files, you can set PAPER=a4 or PAPER=letter"
	@echo "  linkcheck to check all external links for integrity"
	@echo "  upload    to make HTML and LATEX docs and upload it to gna.org"

checkout:
	@if [ ! -d tools ]; then \
	  echo "Checking out tools..."; \
	  svn checkout $(SPHINXROOT)/python/trunk/Doc/tools tools; \
	fi
	@if [ ! -d tools/sphinx ]; then \
	  echo "Checking out Sphinx..."; \
	  svn checkout $(SPHINXROOT)/external/Sphinx-0.6.5/sphinx tools/sphinx; \
	fi
	@if [ ! -d tools/docutils ]; then \
	  echo "Checking out Docutils..."; \
	  svn checkout $(SPHINXROOT)/external/docutils-0.6/docutils tools/docutils; \
	fi
	@if [ ! -d tools/jinja2 ]; then \
	  echo "Checking out Jinja..."; \
	  svn checkout $(SPHINXROOT)/external/Jinja-2.3.1/jinja2 tools/jinja2; \
	fi
	@if [ ! -d tools/pygments ]; then \
	  echo "Checking out Pygments..."; \
	  svn checkout $(SPHINXROOT)/external/Pygments-1.3.1/pygments tools/pygments; \
	fi

update: checkout
	svn update tools
	svn update tools/sphinx
	svn update tools/docutils
	svn update tools/jinja2
	svn update tools/pygments

images:
	-cd $(srcdir)/source/userdoc/images/; make png
	-cd $(srcdir)/source/project/images/; make png
	-cd $(srcdir)/source/matlab/images/; make png
	-cd $(srcdir)/source/python/images/; make png
	-cd $(srcdir)/source/scilab/images/; make png

build: $(srcdir)/source/matlab/cmdref.rst $(srcdir)/source/python/cmdref.rst $(srcdir)/source/scilab/cmdref.rst checkout images
	rm -fr build/$(BUILDER)/_images
	rm -fr build/$(BUILDER)/*.png
	mkdir -p build/$(BUILDER) build/doctrees
	PYTHONPATH=../../interface/src/python && $(SPHINXBUILD) $(ALLSPHINXOPTS)
	@echo

view:
	$(PYTHON) -c "import webbrowser; webbrowser.open('build/html/index.html')"

html: BUILDER = html
html: build
	@echo "Build finished. The HTML pages are in build/html."

htmlview: html
	$(PYTHON) -c "import webbrowser; webbrowser.open('build/html/index.html')"

htmlhelp: BUILDER = htmlhelp
htmlhelp: build
	@echo "Build finished; now you can run HTML Help Workshop with the" \
	      ".hhp project file in build/htmlhelp."

latex: BUILDER = latex
latex: build
	@echo "Build finished; the LaTeX files are in build/latex."
	@echo "Run \`make all-pdf' or \`make all-ps' in that directory to" \
	      "run these through (pdf)latex."

pdf: BUILDER = latex
pdf: build
	(cd build/latex && make all-pdf)
	@echo "Build finished; the pdfLaTeX files are in build/latex."

linkcheck: BUILDER = linkcheck
linkcheck: build
	@echo "Link check complete; look for any errors in the above output " \
	      "or in build/$(BUILDER)/output.txt."

clean:
	-rm -rf build/
	-rm -rf html/
	-cd $(srcdir)/source/userdoc/images/; make clean
	-cd $(srcdir)/source/project/images/; make clean
	-cd $(srcdir)/source/matlab/images/; make clean
	-cd $(srcdir)/source/python/images/; make clean
	-cd $(srcdir)/source/scilab/images/; make clean
	-rm $(srcdir)/source/matlab/cmdref.rst
	-rm $(srcdir)/source/scilab/cmdref.rst
	-rm $(srcdir)/source/python/cmdref.rst

$(srcdir)/source/matlab/cmdref.rst : $(top_srcdir)/interface/src/*.cc $(top_srcdir)/bin/extract_doc
	$(top_srcdir)/bin/extract_doc $(top_srcdir)/interface/src matlab-doc > $(srcdir)/source/matlab/cmdref.rst || (rm -f $(srcdir)/source/matlab/cmdref.rst; /bin/false )

$(srcdir)/source/scilab/cmdref.rst : $(top_srcdir)/interface/src/*.cc $(top_srcdir)/bin/extract_doc
	$(top_srcdir)/bin/extract_doc $(top_srcdir)/interface/src scilab-doc-rst > $(srcdir)/source/scilab/cmdref.rst || (rm -f $(srcdir)/source/scilab/cmdref.rst; /bin/false )

$(srcdir)/source/python/cmdref.rst : $(top_srcdir)/interface/src/*.cc $(top_srcdir)/bin/extract_doc
	$(top_srcdir)/bin/extract_doc $(top_srcdir)/interface/src python-doc > $(srcdir)/source/python/cmdref.rst || (rm -f $(srcdir)/source/python/cmdref.rst; /bin/false )

upload:
	make html
	-rm -fr homepage *.tar.gz
	-mv build/html homepage
	(export srcdir=$(srcdir) && $(top_srcdir)/bin/upload_html --delete homepage)
	-rm -fr getfem-$(DISTVERSION)-docs-html
	-mv homepage getfem-$(DISTVERSION)-docs-html
	tar -cf getfem-$(DISTVERSION)-docs-html.tar getfem-$(DISTVERSION)-docs-html
	gzip -9 getfem-$(DISTVERSION)-docs-html.tar
	(export srcdir=$(srcdir) && $(top_srcdir)/bin/upload_documentation --delete getfem-$(DISTVERSION)-docs-html.tar.gz)
	rm -fr getfem-$(DISTVERSION)-docs-html *.tar.gz
	make pdf
	-mv build/latex/*.pdf .
	(export srcdir=$(srcdir) && $(top_srcdir)/bin/upload_documentation --delete scilab_interface.pdf)
	(export srcdir=$(srcdir) && $(top_srcdir)/bin/upload_documentation --delete python_interface.pdf)
	(export srcdir=$(srcdir) && $(top_srcdir)/bin/upload_documentation --delete matlab_interface.pdf)
	(export srcdir=$(srcdir) && $(top_srcdir)/bin/upload_documentation --delete getfem_userdoc.pdf)
	(export srcdir=$(srcdir) && $(top_srcdir)/bin/upload_documentation --delete getfem_project.pdf)
	(export srcdir=$(srcdir) && $(top_srcdir)/bin/upload_documentation --delete gmm_userdoc.pdf)
	rm -fr *.pdf
	rm -fr build
