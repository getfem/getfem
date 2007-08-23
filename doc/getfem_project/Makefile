all : getfem_project.pdf

FIGS=getfemuserelemf.fig diagram.fig getfemelemelem.fig getfemelemtrans.fig \
     getfemlistsymbols.fig getfemlistsegmentPk.fig getfemlisttriangleP1.fig \
     getfemlisttriangleP2.fig getfemlisttriangleP3.fig                      \
     getfemlisttriangleP6.fig getfemlisttetrahedronP1.fig                   \
     getfemlisttetrahedronP2.fig getfemlisttetrahedronP4.fig                \
     getfemlistquadQ1.fig getfemlistquadQ3.fig getfemlistcubeQ1.fig         \
     getfemlistcubeQ3.fig getfemlistprismP1.fig getfemlistprismP3.fig       \
     getfemlistprismP2P1.fig getfemlistquad8.fig getfemlistsegmenthier.fig  \
     getfemlisttriangleP1comp.fig getfemlisttriangleP1comphier.fig          \
     getfemlistRT0.fig getfemlistnedelec.fig getfemlistsegmenthermite.fig   \
     getfemlistsegmentbubble.fig getfemlisttriangleP1bubble.fig             \
     getfemlisttriangleP2bubble.fig


PDFFIGS=$(FIGS:.fig=.pdf)
PNGFIGS=$(PDFFIGS:.pdf=.png)

.SUFFIXES: .tex .dvi .ps .pdf .eps .fig .png

.fig.eps:
	fig2eps $(@:.eps=.fig)
#	fig2dev -L eps $(@:.eps=.fig) > $@

.eps.pdf:
	epstopdf $(@:.pdf=.eps) --outfile=$@

.pdf.png:
	convert $(@:.png=.pdf) $@

doxygenlinks.tex: updatedoxlinks.py
	python ./updatedoxlinks.py

getfemuserelemf.png: getfemuserelemf.pdf
	convert -resize 500x500 $(@:.png=.pdf) $@

getfemelemelem.png: getfemelemelem.pdf
	convert -resize 500x500 $(@:.png=.pdf) $@

diagram.png: diagram.pdf
	convert $(@:.png=.pdf) $@

TEXOPTS='-interaction=nonstopmode'
TEXMSGFILTER=grep 'LaTeX\|[Ww]arning\|^l\.\|^\!\|^<'

getfem_project.pdf: getfem_project.tex $(PDFFIGS) doxygenlinks.tex
	-pdflatex $(TEXOPTS) getfem_project.tex | $(TEXMSGFILTER) && if (grep Rerun getfem_project.log || grep 'undefined references' getfem_project.log) ; then echo 'RERUN!'; pdflatex $(TEXOPTS) getfem_project.tex | $(TEXMSGFILTER); fi;

html:	getfem_project.tex getfem_project.idx $(PNGFIGS)
	-rm -rf getfem_project/
	hyperlatex getfem_project.tex
	(cd getfem_project && ../cleanup_html_doc.pl)

pdfupload: getfem_project.pdf
	../../bin/upload_documentation getfem_project.pdf
#if [ -d ../../../getfem_html ]; then \
#          cp getfem_project.pdf ../../../getfem_html; \
#fi

htmlupload: html
	cp $(PNGFIGS) getfem_project/
	cp docstyle.css getfem_project/
	cp *.png getfem_project/
	cp next.gif up.gif previous.gif getfem_project/
	../../bin/upload_documentation getfem_project

#tar czvf html_getfem_project.tar.gz getfem_project
#if [ -d ../../../getfem_html ]; then \
#         cp html_getfem_project.tar.gz ../../../getfem_html; \
#fi

all: htmlupload pdfupload

clean:
	-rm -f *.dvi *.log *.toc *.bbl *.aux *.tmp *.ps.gz getfem_project.ps getfem_project.pdf getfem_project.blg getfem_project.out
	-find . -name '*~' -exec rm \{\} \;
	-find . -name '*.bak' -exec rm \{\} \;
