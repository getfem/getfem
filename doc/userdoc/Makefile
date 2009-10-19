all : getfemuser.pdf

FIGS=getfemuserelemf.fig getfemuserelem.fig getfemuserrefine.fig            \
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
     getfemlisttriangleP2bubble.fig getfemlisttriangleP1withP2face.fig      \
     getfemlisttriangleP1bubbleface.fig getfemlisttriangleP1linbubble.fig   \
     getfemlisttriangleP1nonconforming.fig getfemlisttrianglehermite.fig    \
     getfemlistmorley.fig getfemlistargyris.fig getfemlistHCT.fig           \
     getfemlistreducedHCT.fig getfemlistquadc1composite.fig                 \
     getfemlistreducedquadc1composite.fig getfemlisttetrahedronP1bubble.fig \
     getfemlisttetrahedronP2bubble.fig getfemlisttetrahedronP3bubble.fig    \
     getfemlisttetrahedronP1bubbleface.fig getfemlisttetrahedronhermite.fig \
     getfemlistintmethodtriangle1.fig getfemlistintmethodtriangle2.fig      \
     getfemlistintmethodtriangle3.fig getfemlistintmethodtriangle4.fig      \
     getfemlistintmethodtriangle5.fig getfemlistintmethodtriangle6.fig      \
     getfemlistintmethodtriangle7.fig getfemlistintmethodquad2.fig          \
     getfemlistintmethodquad3.fig getfemlistintmethodquad5.fig              \
     getfemlistintmethodtriangle2comp.fig getfemuserlinsysDir.fig           \
     getfemlistintmethodtetrahedron1.fig getfemuserlinearsys.fig            \
     getfemlistintmethodtetrahedron2.fig                                    \
     getfemlistintmethodtetrahedron3.fig                                    \
     getfemlistintmethodtetrahedron5.fig

PDFFIGS=$(FIGS:.fig=.pdf)
PNGFIGS=$(PDFFIGS:.pdf=.png)

.SUFFIXES: .tex .dvi .ps .pdf .eps .fig .png

.fig.eps:
	../../bin/fig2eps $(@:.eps=.fig)
#	fig2dev -L eps $(@:.eps=.fig) > $@

.eps.pdf:
	epstopdf $(@:.pdf=.eps) --outfile=$@

.pdf.png:
	convert $(@:.png=.pdf) $@

doxygenlinks.tex: updatedoxlinks.py
	python ./updatedoxlinks.py

getfemuserelemf.png: getfemuserelemf.pdf
	convert -density 80x80 $(@:.png=.pdf) $@

getfemuserelem.png: getfemuserelem.pdf
	convert -density 80x80 $(@:.png=.pdf) $@

getfemuserlinsysDir.png: getfemuserlinsysDir.pdf
	convert -density 100x100 $(@:.png=.pdf) $@

getfemuserlinearsys.png: getfemuserlinearsys.pdf
	convert -density 100x100 $(@:.png=.pdf) $@


TEXOPTS='-interaction=nonstopmode'
TEXMSGFILTER=grep 'LaTeX\|[Ww]arning\|^l\.\|^\!\|^<'

getfemuser.pdf: getfemuser.tex $(PDFFIGS) doxygenlinks.tex
	-pdflatex $(TEXOPTS) getfemuser.tex | $(TEXMSGFILTER) && if (grep Rerun getfemuser.log || grep 'undefined references' getfemuser.log) ; then echo 'RERUN!'; pdflatex $(TEXOPTS) getfemuser.tex | $(TEXMSGFILTER); fi;

#getfemuser.dvi : getfemuser.tex
#	latex getfemuser.tex; makeindex getfemuser.idx; latex getfemuser.tex; makeindex getfemuser.idx; latex getfemuser.tex

#getfemuser.ps : getfemuser.dvi
#	dvips getfemuser -z -Pamz -Pcmz -o
#	ps2pdf getfemuser.ps getfemuser.pdf; \
#	cp getfemuser.ps getfemuser.pdf ../../../getfem_html/;

html:	getfemuser.tex getfemuser.idx $(PNGFIGS)
	-rm -rf getfemuser/
	hyperlatex getfemuser.tex
	(cd getfemuser && ../cleanup_html_doc.pl)

pdfupload: getfemuser.pdf
	../../bin/upload_documentation getfemuser.pdf
#if [ -d ../../../getfem_html ]; then \
#          cp getfemuser.pdf ../../../getfem_html; \
#fi

htmlupload: html
	cp $(PNGFIGS) getfemuser/
	cp docstyle.css getfemuser/
	cp logogetfem.png logo_getfem_small.png getfemuser/
	cp getfemuserrefine.png getfemuser/
	cp next.gif up.gif previous.gif getfemuser/
	../../bin/upload_documentation getfemuser

#tar czvf html_getfemuser.tar.gz getfemuser
#if [ -d ../../../getfem_html ]; then \
#         cp html_getfemuser.tar.gz ../../../getfem_html; \
#fi

all: pdfupload htmlupload

clean:
	-rm -f *.dvi *.log *.toc *.bbl *.aux *.tmp *.ps.gz *.pdf getfemuser.ps getfemuser.blg getfemuser.out getfemuser.idx getfemuser.brf
	-find . -name '*~' -exec rm \{\} \;
	-find . -name '*.bak' -exec rm \{\} \;
