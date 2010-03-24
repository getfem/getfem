FIGS=diagram.fig \
     getfemelemelem.fig \
     getfemtransgeo.fig

EPSFIGS=$(FIGS:.fig=.eps)
PNGFIGS=$(FIGS:.fig=.png)

.SUFFIXES: .fig .eps .png

.fig.eps:
	../../../../../bin/fig2eps $(@:.eps=.fig)

.eps.png:
	convert $(@:.png=.eps) $@

png: $(PNGFIGS)

clean:
	-rm -rf *.aux
	-rm -rf *.dvi
	-rm -rf *.log
	-rm -rf *.bak
	-rm -rf $(EPSFIGS)
	-rm -rf $(PNGFIGS)
