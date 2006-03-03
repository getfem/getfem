all : getfemuser.pdf

FIGS=getfemuserelemf.fig getfemuserelem.fig
PDFFIGS=$(FIGS:.fig=.pdf)
PNGFIGS=$(PDFFIGS:.pdf=.png)

.SUFFIXES: .tex .dvi .ps .pdf .eps .fig .png

.fig.eps:
	fig2dev -L eps $(@:.eps=.fig) > $@

.eps.pdf:
	epstopdf $(@:.pdf=.eps) --outfile=$@

getfemuserelemf.png: getfemuserelemf.pdf
	convert -resize 500x500 $(@:.png=.pdf) $@

getfemuserelem.png: getfemuserelem.pdf
	convert -resize 500x500 $(@:.png=.pdf) $@

TEXOPTS='-interaction=nonstopmode'
TEXMSGFILTER=grep 'LaTeX\|[Ww]arning\|^l\.\|^\!\|^<'

getfemuser.pdf: getfemuser.tex $(PDFFIGS)
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

pdf: getfemuser.pdf
	if [ -d ../../../getfem_html ]; then \
           cp getfemuser.pdf ../../../getfem_html; \
        fi

htmltar: html
	cp $(PNGFIGS) getfemuser/
	cp docstyle.css getfemuser/
	cp logogetfem.png logo_getfem_small.png getfemuser/
	cp getfemuserrefine.png getfemuser/
	cp next.gif up.gif previous.gif getfemuser/
	tar czvf html_getfemuser.tar.gz getfemuser
	if [ -d ../../../getfem_html ]; then \
           cp html_getfemuser.tar.gz ../../../getfem_html; \
        fi

all: htmltar pdf

clean:
	-rm -f *.dvi *.log *.toc *.bbl *.aux *.tmp *.ps.gz getfemuser.ps getfemuser.pdf getfemuser.blg getfemuser.out
	-find . -name '*~' -exec rm \{\} \;
	-find . -name '*.bak' -exec rm \{\} \;
