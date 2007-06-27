OUTPUT=$@
TEXOPTS='-interaction=nonstopmode'
TEXMSGFILTER=grep 'LaTeX\|[Ww]arning\|^l\.\|^\!\|^<'
FIGS=hierarchy.fig



PDFFIGS=$(FIGS:.fig=.pdf) fempk51.pdf
PNGFIGS=$(PDFFIGS:.pdf=.png)

BASEPNGFIGS=donut_small.png tripodvonmiseswithmesh_small.png logogetfem.png logo_getfem_small.png

#png not in html doc
OTHERPNGFIGS=donut.png tripodvonmiseswithmesh.png logogetfemwhitebg.png

.SUFFIXES: .tex .dvi .ps .pdf .eps .fig .png

.fig.pdf:
	fig2dev -L eps $(@:.pdf=.fig) > $(@:.pdf=.eps)
	epstopdf $(@:.pdf=.eps) --outfile=$@

all : pdfupload htmlupload
	if [ -d ../../getfem_html ]; then \
           cp getfem_matlab.pdf ../../getfem_html; \
        fi

.eps.pdf:
	epstopdf $(@:.pdf=.eps) --outfile=$@

hierarchy.png: hierarchy.pdf
	convert -resize 400x400 $(@:.png=.pdf) $@

fempk51.png: fempk51.pdf
	convert -resize 300x200 $(@:.png=.pdf) $@

#donut.pdf:
#	convert -resize 500x500 $(@:.pdf=.png) $@

#tripodvonmiseswithmesh.pdf:
#	convert -resize 500x500 $(@:.pdf=.png) $@

#.eps2 ; mv $@.eps2 $@


#dvi::
#	latex $(TEXOPTS) gfm.tex | $(TEXMSGFILTER)

gfm.tex : getfemmatlab.tex
	perl ../bin/latexize_mcode.pl < getfemmatlab.tex > gfm.tex;


demolaplacian.tex : ../tests/matlab/demo_laplacian.m
	../bin/latexize_mfile.sh ../tests/matlab/demo_laplacian.m demolaplacian.tex

demotripod.tex: ../tests/matlab/demo_tripod.m
	../bin/latexize_mfile.sh ../tests/matlab/demo_tripod.m demotripod.tex

demorefine.tex: ../tests/matlab/demo_refine.m
	../bin/latexize_mfile.sh ../tests/matlab/demo_refine.m demorefine.tex

#gfm.dvi : gfm.tex $(PDFFIGS)
#	-latex $(TEXOPTS) gfm.tex | $(TEXMSGFILTER) && if (grep Rerun gfm.log || grep 'undefined references' gfm.log) ; then echo 'RERUN!'; latex $(TEXOPTS) gfm.tex | $(TEXMSGFILTER); fi;
#
#getfem_matlab.ps : gfm.dvi
#	dvips gfm -z -Pamz -Pcmz -o getfem_matlab.ps -p1

getfem_matlab.pdf : $(PDFFIGS) $(PNGFIGS) gfm.tex gfm.idx demolaplacian.tex demotripod.tex demorefine.tex
	-pdflatex $(TEXOPTS) gfm.tex | $(TEXMSGFILTER) && if (grep Rerun gfm.log || grep 'undefined references' gfm.log) ; then echo 'RERUN!'; pdflatex $(TEXOPTS) gfm.tex | $(TEXMSGFILTER); fi;
	mv gfm.pdf getfem_matlab.pdf

gfm.idx : gfm.tex
	touch -a gfm.idx
	makeindex gfm.idx

.PHONY : clean
clean:
	-rm -f *.dvi *.log *.toc *.bbl *.aux *.tmp *.ps.gz gfm.ps gfm.pdf gfm.blg gfm.out
	-find . -name '*~' -exec rm \{\} \;
	-find . -name '*.bak' -exec rm \{\} \;

html:	gfm.tex gfm.idx
	-rm -rf gfm/
	hyperlatex gfm.tex
	( cd gfm && ../cleanup_html_doc.pl ) && mv gfm getfem_matlab

htmlupload: html
	cp $(PNGFIGS) $(BASEPNGFIGS) getfem_matlab/
	cp docstyle.css getfem_matlab/
	cp next.gif up.gif previous.gif getfem_matlab/
	../../bin/upload_documentation getfem_matlab
	../../bin/upload_documentation getfem_python_reference.html

pdfupload: getfem_matlab.pdf
	../../bin/upload_documentation getfem_matlab.pdf

#tar czvf html_gfm.tar.gz gfm
#if [ -d ../../../getfem_html ]; then \
#          cp html_gfm.tar.gz ../../../getfem_html; \
#       fi
