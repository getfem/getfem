all : gmmuser.pdf

# FIGS=gmmuserelemf.fig gmmuserelem.fig
PDFFIGS=$(FIGS:.fig=.pdf)
PNGFIGS=$(PDFFIGS:.pdf=.png) gmmlogo.png gmmlogo_small.png

.SUFFIXES: .tex .dvi .ps .pdf .eps .fig .png

.fig.eps:
	fig2dev -L eps $(@:.eps=.fig) > $@

.eps.pdf:
	epstopdf $(@:.pdf=.eps) --outfile=$@

gmmuser.idx : gmmuser.tex
	touch gmmuser.idx

gmmuser.ilg: gmmuser.idx
	makeindex gmmuser.idx

# gmmuserelemf.png: gmmuserelemf.pdf
#	convert -resize 500x500 $(@:.png=.pdf) $@
#
# gmmuserelem.png: gmmuserelem.pdf
#	convert -resize 500x500 $(@:.png=.pdf) $@

TEXOPTS='-interaction=nonstopmode'
TEXMSGFILTER=grep 'LaTeX\|[Ww]arning\|^l\.\|^\!\|^<'

gmmuser.pdf: gmmuser.tex $(PDFFIGS) gmmuser.ilg
	-pdflatex $(TEXOPTS) gmmuser.tex | $(TEXMSGFILTER) && if (grep Rerun gmmuser.log || grep 'undefined references' gmmuser.log) ; then echo 'RERUN!'; pdflatex $(TEXOPTS) gmmuser.tex | $(TEXMSGFILTER); fi;

html:	gmmuser.tex gmmuser.ilg $(PNGFIGS)
	-rm -rf gmmuser/
	hyperlatex gmmuser.tex
	(cd gmmuser && ../cleanup_html_doc.pl)

pdfupload: gmmuser.pdf
	../../bin/upload_documentation gmmuser.pdf

#if [ -d ../../../getfem_html ]; then \
#          cp gmmuser.pdf ../../../getfem_html; \
#fi

htmlupload: html
	cp $(PNGFIGS) gmmuser/
	cp docstyle.css gmmuser/
	cp next.gif up.gif previous.gif gmmuser/
	../../bin/upload_documentation gmmuser

#tar czvf html_gmmuser.tar.gz gmmuser
#if [ -d ../../../getfem_html ]; then \
#          cp html_gmmuser.tar.gz ../../../getfem_html; \
#       fi

all: htmlupload pdfupload

clean:
	-rm -f *.dvi *.log *.toc *.bbl *.aux *.tmp *.ps.gz gmmuser.ps gmmuser.pdf gmmuser.blg gmmuser.out
	-find . -name '*~' -exec rm \{\} \;
	-find . -name '*.bak' -exec rm \{\} \;
