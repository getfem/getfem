png: eps
	convert -density 100x100 getfemlistsegmentPk.eps getfemlistsegmentPk.png
	convert -density 100x100 getfemuserrefine.eps getfemuserrefine.png
	convert -density 100x100 getfemlistsymbols.eps getfemlistsymbols.png
	convert -density 100x100 getfemuserlinsysDir.eps getfemuserlinsysDir.png
	convert -density 100x100 getfemuserlinearsys.eps getfemuserlinearsys.png
	convert -density 100x100 getfemuserelem.eps getfemuserelem.png
	convert -density 100x100 getfemuserelemf.eps getfemuserelemf.png

.fig.eps:
	echo $@

eps:
	../../../../../bin/fig2eps ../../../../userdoc/getfemlistsegmentPk.fig
	../../../../../bin/fig2eps ../../../../userdoc/getfemuserrefine.fig
	../../../../../bin/fig2eps ../../../../userdoc/getfemlistsymbols.fig
	../../../../../bin/fig2eps ../../../../userdoc/getfemuserlinsysDir.fig
	../../../../../bin/fig2eps ../../../../userdoc/getfemuserlinearsys.fig
	../../../../../bin/fig2eps ../../../../userdoc/getfemuserelem.fig
	../../../../../bin/fig2eps ../../../../userdoc/getfemuserelemf.fig

clean:
	-rm -rf *.png
	-rm -rf *.aux
	-rm -rf *.dvi
	-rm -rf *.log
	-rm -rf *.eps
