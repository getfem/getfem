FIGS=hierarchy.fig

EPSFIGS=$(FIGS:.fig=.eps)
PNGFIGS=$(FIGS:.fig=.png)

.SUFFIXES: .fig .eps .png

.fig.eps:
	../../../../../bin/fig2eps $(@:.eps=.fig)

.eps.png:
	convert -density 100x100 $(@:.png=.eps) $@

png: $(PNGFIGS)

clean:
	-rm -rf $(EPSFIGS)
	-rm -rf $(PNGFIGS)
