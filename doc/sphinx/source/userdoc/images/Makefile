FIGS=getfemlistsymbols.fig    \
     getfemlistsegmentPk.fig  \
     getfemuserelem.fig       \
     getfemuserelemf.fig      \
     getfemuserlinearsys.fig  \
     getfemuserlinsysDir.fig  \
     getfemuserrefine.fig     \
     getfemlisttriangleP1.fig \
     getfemlisttriangleP3.fig \
     getfemlisttriangleP2.fig \
     getfemlisttriangleP6.fig

PNGFIGS=$(FIGS:.fig=.png)

.SUFFIXES: .fig .eps .png

%.fig:
	-cp ../../../../userdoc/$(@) ./

.fig.eps:
	../../../../../bin/fig2eps $(@:.eps=.fig)

.eps.png:
	convert -density 100x100 $(@:.png=.eps) $@

getfemlisttriangleP1.png: getfemlisttriangleP1.eps
	convert -density 50x50 getfemlisttriangleP1.eps getfemlisttriangleP1.png
	rm -f getfemlisttriangleP1.eps
getfemlisttriangleP2.png: getfemlisttriangleP2.eps
	convert -density 50x50 getfemlisttriangleP2.eps getfemlisttriangleP2.png
	rm -f getfemlisttriangleP1.eps
getfemlisttriangleP3.png: getfemlisttriangleP3.eps
	convert -density 50x50 getfemlisttriangleP3.eps getfemlisttriangleP3.png
	rm -f getfemlisttriangleP1.eps
getfemlisttriangleP6.png: getfemlisttriangleP6.eps
	convert -density 50x50 getfemlisttriangleP6.eps getfemlisttriangleP6.png
	rm -f getfemlisttriangleP1.eps

png: $(PNGFIGS)

clean:
	-rm -rf *.fig
	-rm -rf *.eps
	-rm -rf *.png
	-rm -rf *.aux
	-rm -rf *.dvi
	-rm -rf *.log
