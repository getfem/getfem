all : getfemuser.ps

getfemuser.dvi : getfemuser.tex
	latex getfemuser.tex; latex getfemuser.tex

getfemuser.ps : getfemuser.dvi
	dvips getfemuser -z -Pamz -Pcmz -o -p1
	cp getfemuser.ps $(HOME)/public_html/getfem
	cd $(HOME)/public_html/getfem; make

