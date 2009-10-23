FIGS=diagram.fig \
     getfemelemelem.fig \
     getfemtransgeo.fig

PNGFIGS=$(FIGS:.fig=.png)

.SUFFIXES: .fig .eps .png

%.fig:
	-cp ../../../../getfem_project/$(@) ./

.fig.eps:
	../../../../../bin/fig2eps $(@:.eps=.fig)

.eps.png:
	convert -density 80x80 $(@:.png=.eps) $@

png: $(PNGFIGS)

clean:
	-rm -rf *.eps
	-rm -rf *.png
	-rm -rf *.aux
	-rm -rf *.dvi
	-rm -rf *.log
