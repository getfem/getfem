all : getfem_interpole.ps

getfem_interpole.dvi : getfem_interpole.tex
	latex getfem_interpole.tex; latex getfem_interpole.tex

getfem_interpole.ps : getfem_interpole.dvi
	dvips getfem_interpole -z -Pamz -Pcmz -o -p1
	ps2pdf getfem_interpole.ps getfem_interpole.pdf; \
	cp getfem_interpole.ps getfem_interpole.pdf ../../../getfem_html/

