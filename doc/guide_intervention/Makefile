

all : guide_intervention.ps



guide_intervention.ps : guide_intervention.tex files_structure.eps \
	 	getfemlab_orga.eps getfem_orga.eps guide_intervention_pg.eps
	latex guide_intervention; latex guide_intervention
	dvips -o guide_inter.ps -p1 guide_intervention
#	mpage -2 guide_inter.ps > guide_intervention.ps # marche pas ...
	cp guide_inter.ps guide_intervention.ps
	rm -f guide_inter.ps *.toc *.aux *.log *.bak *.dvi

clean :
	rm -f guide_intervention.ps *.toc *.aux *.dvi *.log *.bak *.dvi

