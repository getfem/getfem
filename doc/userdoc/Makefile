all : getfemuser.ps

getfemuser.dvi : getfemuser.tex
	latex getfemuser.tex; latex getfemuser.tex

getfemuser.ps : getfemuser.dvi
	dvips getfemuser -z -Pamz -Pcmz -o -p1
	ps2pdf getfemuser.ps getfemuser.pdf; \
	cp getfemuser.ps getfemuser.pdf ../../../getfem_html/;


