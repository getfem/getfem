

all : guide_intervention.ps



guide_intervention.ps : guide_intervention.tex files_structure.eps \
	 	getfemlab_orga.eps getfem_orga.eps guide_intervention_pg.eps
	latex guide_intervention; latex guide_intervention
	dvips -o guide_inter1.ps -p1 guide_intervention
	pstops "2:0L@.7(21cm,0)+1L@.7(21cm,14.85cm)" guide_inter1.ps guide_inter2.ps
	pstops "2:0,1LL@1(21cm,29.7cm)" guide_inter2.ps guide_intervention.ps
#	mpage -2 guide_inter.ps > guide_intervention.ps # marche pas ...
#	cp guide_inter.ps guide_intervention.ps
	rm -f guide_inter?.ps *.toc *.aux *.log *.bak *.dvi

clean :
	rm -f guide_intervention.ps *.toc *.aux *.dvi *.log *.bak *.dvi

