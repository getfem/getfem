all : getfem_interpole.ps

getfem_interpole.dvi : getfem_interpole.tex
	latex getfem_interpole.tex; latex getfem_interpole.tex

getfem_interpole.ps : getfem_interpole.dvi
	dvips getfem_interpole -z -Pamz -Pcmz -o -p1
	cp getfem_interpole.ps $(HOME)/public_html/getfem
	cd $(HOME)/public_html/getfem; make

