all : getfemelem.ps getfemlist.ps

getfemlist_segment_Pk.eps : getfemlist_segment_Pk.fig
	fig2eps getfemlist_segment_Pk.fig

getfemelem.dvi : getfemelem.tex
	latex getfemelem.tex; latex getfemelem.tex

getfemelem.ps : getfemelem.dvi
	dvips getfemelem -z -Pamz -Pcmz -o -p1
	cp getfemelem.ps $(HOME)/public_html/getfem
	cd $(HOME)/public_html/getfem; make

getfemlist.dvi : getfemlist.tex getfemlist_segment_Pk.eps
	latex getfemlist.tex; latex getfemlist.tex

getfemlist.ps : getfemlist.dvi
	dvips getfemlist -z -Pamz -Pcmz -o -p1
	cp getfemlist.ps $(HOME)/public_html/getfem
	cd $(HOME)/public_html/getfem; make
