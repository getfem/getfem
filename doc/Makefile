


all : doc.ps guide_intervention/guide_intervention.ps

doc.ps : ../include/*.h
	./do_doc++
	rm -f *.toc *.aux *.log *.bak *.dvi doc.tex doc.dvi

guide_intervention/guide_intervention.ps : guide_intervention/guide_intervention.tex
	cd guide_intervention; make

