all : getfemelem.ps getfemlist.ps

getfemlist_segment_Pk.eps : getfemlist_segment_Pk.fig
	fig2eps getfemlist_segment_Pk.fig

getfemlist_extrans.eps : getfemlist_extrans.fig
	fig2eps getfemlist_extrans.fig

getfemelem.dvi : getfemelem.tex getfemelemint.tex getfemeleminta.tex getfemelemfem.tex
	latex getfemelem.tex; latex getfemelem.tex

getfemelem.ps : getfemelem.dvi
	dvips getfemelem -z -Pamz -Pcmz -o 
	ps2pdf getfemelem.ps getfemelem.pdf; \
	cp getfemelem.ps getfemelem.pdf ../../../getfem_html/;

# dvips getfemelem -z -Pamz -Pcmz -o -p1

getfemlist.dvi : getfemlist.tex getfemlist_segment_Pk.eps getfemlist_extrans.eps 
	latex getfemlist.tex; latex getfemlist.tex

getfemlist.ps : getfemlist.dvi
	dvips getfemlist -z -Pamz -Pcmz -o 
	ps2pdf getfemlist.ps getfemlist.pdf; \
	cp getfemlist.ps getfemlist.pdf ../../../getfem_html/;

#	dvips getfemlist -z -Pamz -Pcmz -o -p1
