FIGS=getfemlistargyris.fig                   \
     getfemlistcubeQ1.fig                    \
     getfemlistcubeQ3.fig                    \
     getfemlistHCT.fig                       \
     getfemlistintmethodquad2.fig            \
     getfemlistintmethodquad3.fig            \
     getfemlistintmethodquad5.fig            \
     getfemlistintmethodtetrahedron1.fig     \
     getfemlistintmethodtetrahedron2.fig     \
     getfemlistintmethodtetrahedron3.fig     \
     getfemlistintmethodtetrahedron5.fig     \
     getfemlistintmethodtriangle1.fig        \
     getfemlistintmethodtriangle2comp.fig    \
     getfemlistintmethodtriangle2.fig        \
     getfemlistintmethodtriangle3.fig        \
     getfemlistintmethodtriangle4.fig        \
     getfemlistintmethodtriangle5.fig        \
     getfemlistintmethodtriangle6.fig        \
     getfemlistintmethodtriangle7.fig        \
     getfemlistmorley.fig                    \
     getfemlistnedelec.fig                   \
     getfemlistprismP1.fig                   \
     getfemlistprismP2P1.fig                 \
     getfemlistprismP3.fig                   \
     getfemlistquad8.fig                     \
     getfemlistquadc1composite.fig           \
     getfemlistquadQ1.fig                    \
     getfemlistquadQ3.fig                    \
     getfemlistreducedHCT.fig                \
     getfemlistreducedquadc1composite.fig    \
     getfemlistRT0.fig                       \
     getfemlistsegmentbubble.fig             \
     getfemlistsegmenthermite.fig            \
     getfemlistsegmenthier.fig               \
     getfemlistsegmentPk.fig                 \
     getfemlistsymbols.fig                   \
     getfemlisttetrahedronhermite.fig        \
     getfemlisttetrahedronP1bubbleface.fig   \
     getfemlisttetrahedronP1bubble.fig       \
     getfemlisttetrahedronP1.fig             \
     getfemlisttetrahedronP2bubble.fig       \
     getfemlisttetrahedronP2.fig             \
     getfemlisttetrahedronP3bubble.fig       \
     getfemlisttetrahedronP4.fig             \
     getfemlisttrianglehermite.fig           \
     getfemlisttriangleP1bubbleface.fig      \
     getfemlisttriangleP1bubble.fig          \
     getfemlisttriangleP1comp.fig            \
     getfemlisttriangleP1comphier.fig        \
     getfemlisttriangleP1.fig                \
     getfemlisttriangleP1linbubble.fig       \
     getfemlisttriangleP1nonconforming.fig   \
     getfemlisttriangleP1withP2face.fig      \
     getfemlisttriangleP2bubble.fig          \
     getfemlisttriangleP2.fig                \
     getfemlisttriangleP3.fig                \
     getfemlisttriangleP6.fig                \
     getfemuserelemf.fig                     \
     getfemuserelem.fig                      \
     getfemuserlinearsys.fig                 \
     getfemuserlinsysDir.fig                 \
     getfemuserrefine.fig



PNGFIGS=$(FIGS:.fig=.png)

.SUFFIXES: .fig .eps .png

.fig.eps:
	../../../../../bin/fig2eps $(@:.eps=.fig)

.eps.png:
	convert -density 100x100 $(@:.png=.eps) $@

getfemlisttriangleP1.png: getfemlisttriangleP1.eps
	convert getfemlisttriangleP1.eps getfemlisttriangleP1.png

getfemlisttriangleP2.png: getfemlisttriangleP2.eps
	convert getfemlisttriangleP2.eps getfemlisttriangleP2.png

getfemlisttriangleP3.png: getfemlisttriangleP3.eps
	convert getfemlisttriangleP3.eps getfemlisttriangleP3.png

getfemlisttriangleP6.png: getfemlisttriangleP6.eps
	convert getfemlisttriangleP6.eps getfemlisttriangleP6.png


png: $(PNGFIGS)

clean:
	-rm -rf *.eps
	-rm -rf *.png
	-rm -rf *.aux
	-rm -rf *.dvi
	-rm -rf *.log
