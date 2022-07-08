#  Copyright (C) 1999-2020 Yves Renard
#
#  This file is a part of GetFEM
#
#  GetFEM  is  free software;  you  can  redistribute  it  and/or modify it
#  under  the  terms  of the  GNU  Lesser General Public License as published
#  by  the  Free Software Foundation;  either version 3 of the License,  or
#  (at your option) any later version along with the GCC Runtime Library
#  Exception either version 3.1 or (at your option) any later version.
#  This program  is  distributed  in  the  hope  that it will be useful,  but
#  WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
#  or  FITNESS  FOR  A PARTICULAR PURPOSE.  See the GNU Lesser General Public
#  License and GCC Runtime Library Exception for more details.
#  You  should  have received a copy of the GNU Lesser General Public License
#  along  with  this program;  if not, write to the Free Software Foundation,
#  Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301, USA.

# Makefile for Sphinx documentation
#

# You can set these variables from the command line.
PYTHON        = python3
SPHINXROOT    = http://svn.python.org/projects
LANGUAGE      = en
SPHINXOPTS    = -D language=$(LANGUAGE)
# SPHINXBUILD   = $(PYTHON) tools/sphinx-build.py
# SPHINXBUILD   = sphinx-build
SPHINXBUILD   = python3 -msphinx
SPHINXINTL    = sphinx-intl
PAPER         =
SOURCES       = $(srcdir)/source
DISTVERSION   = @VERSION@


# Internal variables.
PAPEROPT_a4     = -D latex_paper_size=a4
PAPEROPT_letter = -D latex_paper_size=letter
ALLSPHINXOPTS   = -b $(BUILDER) -d build/doctrees $(PAPEROPT_$(PAPER)) \
                  $(SPHINXOPTS) $(SOURCES) build/$(BUILDER)
I18NSPHINXOPTS  = $(SPHINXOPTS) ./source

.PHONY: help checkout update images build view html htmlview htmlhelp latex pdf linkcheck clean upload gettext


help:
	@echo "Please use \`make <target>' where <target> is one of"
	@echo "  view      to view standalone HTML files"
	@echo "  html      to make standalone HTML files"
	@echo "  htmlview  to make and view standalone HTML files"
	@echo "  htmlhelp  to make HTML files and a HTML help project"
	@echo "  latex     to make LaTeX files, you can set PAPER=a4 or PAPER=letter"
	@echo "  pdf       to make pdfLaTeX files, you can set PAPER=a4 or PAPER=letter"
	@echo "  linkcheck to check all external links for integrity"
	@echo "  upload    to make HTML and LATEX docs and upload it to Savannah.org"

# checkout:
#	@if [ ! -d tools ]; then \
#	  echo "Checking out tools..."; \
#	  svn checkout $(SPHINXROOT)/python/trunk/Doc/tools tools; \
#	fi
#	@if [ ! -d tools/sphinx ]; then \
#	  echo "Checking out Sphinx..."; \
#	  svn checkout $(SPHINXROOT)/external/Sphinx-1.2/sphinx tools/sphinx; \
#	fi
#	@if [ ! -d tools/docutils ]; then \
#	  echo "Checking out Docutils..."; \
#	  svn checkout $(SPHINXROOT)/external/docutils-0.6/docutils tools/docutils; \
#	fi
#	@if [ ! -d tools/jinja2 ]; then \
#	  echo "Checking out Jinja..."; \
#	  svn checkout $(SPHINXROOT)/external/Jinja-2.3.1/jinja2 tools/jinja2; \
#	fi
#	@if [ ! -d tools/pygments ]; then \
#	  echo "Checking out Pygments..."; \
#	  svn checkout $(SPHINXROOT)/external/Pygments-1.6/pygments tools/pygments; \
#	fi
#
#update: checkout
#	svn update tools
#	svn update tools/sphinx
#	svn update tools/docutils
#	svn update tools/jinja2
#	svn update tools/pygments

images:
	-cd $(srcdir)/source/userdoc/images/; make png
	-cd $(srcdir)/source/project/images/; make png
	-cd $(srcdir)/source/matlab_octave/images/; make png
	-cd $(srcdir)/source/python/images/; make png
	-cd $(srcdir)/source/scilab/images/; make png

build: $(srcdir)/source/matlab_octave/cmdref.rst $(srcdir)/source/python/cmdref.rst $(srcdir)/source/scilab/cmdref.rst checkout images
	echo # $(SPHINXINTL) build # sphinx-intl produces a segmentation fault ...
	echo # rm -fr build/$(BUILDER)/_images
	echo # rm -fr build/$(BUILDER)/*.png
	mkdir -p build/$(BUILDER) build/doctrees
	echo # PYTHONPATH=../../interface/src/python LD_PRELOAD=../../src/.libs/libgetfem.so $(SPHINXBUILD) $(ALLSPHINXOPTS)
	PYTHONPATH=../../interface/src/python $(SPHINXBUILD) $(ALLSPHINXOPTS)
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
	@echo "Run \`make all-pdf' in that directory to" \
	      "run these through (pdf)latex."

pdf: BUILDER = latex
pdf: build
	(cd build/latex && make all-pdf)
	@echo "Build finished; the pdfLaTeX files are in build/latex."

linkcheck: BUILDER = linkcheck
linkcheck: build
	@echo "Link check complete; look for any errors in the above output " \
	      "or in build/$(BUILDER)/output.txt."
gettext:
	$(SPHINXBUILD) -b gettext $(I18NSPHINXOPTS) locale
	$(SPHINXINTL) update -p locale -l $(LANGUAGE)
	@echo
	@echo "Build finished. The message catalogs are in locale."

clean:
	-rm -rf build/
	-rm -rf html/
	-cd $(srcdir)/source/userdoc/images/ && make clean
	-cd $(srcdir)/source/project/images/ && make clean
	-cd $(srcdir)/source/matlab_octave/images/ && make clean
	-cd $(srcdir)/source/python/images/ && make clean
	-cd $(srcdir)/source/scilab/images/ && make clean
	-rm -f $(srcdir)/source/matlab_octave/cmdref.rst
	-rm -f $(srcdir)/source/scilab/cmdref.rst
	-rm -f $(srcdir)/source/python/cmdref.rst

$(srcdir)/source/matlab_octave/cmdref.rst : $(top_srcdir)/interface/src/*.cc $(top_srcdir)/bin/extract_doc $(top_srcdir)/bin/split_cmdref
	($(top_srcdir)/bin/extract_doc $(top_srcdir)/interface/src matlab-doc > $(srcdir)/source/matlab_octave/cmdref.rst &&  $(top_srcdir)/bin/split_cmdref $(srcdir)/source/matlab_octave) || (rm -f $(srcdir)/source/matlab_octave/cmdref.rst; /bin/false )

$(srcdir)/source/scilab/cmdref.rst : $(top_srcdir)/interface/src/*.cc $(top_srcdir)/bin/extract_doc $(top_srcdir)/bin/split_cmdref
	($(top_srcdir)/bin/extract_doc $(top_srcdir)/interface/src scilab-doc-rst > $(srcdir)/source/scilab/cmdref.rst  &&  $(top_srcdir)/bin/split_cmdref $(srcdir)/source/scilab) || (rm -f $(srcdir)/source/scilab/cmdref.rst; /bin/false )

$(srcdir)/source/python/cmdref.rst : $(top_srcdir)/interface/src/*.cc $(top_srcdir)/bin/extract_doc $(top_srcdir)/bin/split_cmdref
	($(top_srcdir)/bin/extract_doc $(top_srcdir)/interface/src python-doc > $(srcdir)/source/python/cmdref.rst  &&  $(top_srcdir)/bin/split_cmdref $(srcdir)/source/python) || (rm -f $(srcdir)/source/python/cmdref.rst; /bin/false )

upload:
	rm -fr build
	make html
	-rm -fr homepage *.tar.gz
	-mv build/html homepage
	-(export srcdir=$(srcdir) && $(top_srcdir)/bin/upload_html --delete homepage)
	-test ! -d "getfem_site" && git clone "git@plmlab.math.cnrs.fr:renard/getfem_site.git"
	-(cd homepage && rsync --delete  -rlptDv . ../getfem_site/public)
	-(cd getfem_site && git add . && git commit -m "update" && git push origin master)
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
	(export srcdir=$(srcdir) && $(top_srcdir)/bin/upload_documentation --delete matlab_octave_interface.pdf)
	(export srcdir=$(srcdir) && $(top_srcdir)/bin/upload_documentation --delete getfem_userdoc.pdf)
	(export srcdir=$(srcdir) && $(top_srcdir)/bin/upload_documentation --delete getfem_project.pdf)
	(export srcdir=$(srcdir) && $(top_srcdir)/bin/upload_documentation --delete gmm_userdoc.pdf)
	rm -fr *.pdf


EXTRA_DIST = \
	source/userdoc/images/Makefile                  \
	source/project/images/Makefile                  \
	source/matlab_octave/images/Makefile            \
	source/python/images/Makefile                   \
	source/scilab/images/Makefile
