#  Copyright (C) 1999-2020 Yves Renard
#
#  This file is a part of GetFEM++
#
#  GetFEM++  is  free software;  you  can  redistribute  it  and/or modify it
#  under  the  terms  of the  GNU  Lesser General Public License as published
#  by  the  Free Software Foundation;  either version 3 of the License,  or
#  (at your option) any later version along with the GCC Runtime Library
#  Exception either version 3.1 or (at your option) any later version.
#  This program  is  distributed  in  the  hope  that it will be useful,  but
#  WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
#  or  FITNESS  FOR  A PARTICULAR PURPOSE.  See the GNU Lesser General Public
#  License and GCC Runtime Library Exception for more details.
#  You  should  have received a copy of the GNU Lesser General Public License
#  along  with  this program;  if not, write to the Free Software Foundation,
#  Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301, USA.

FIGS=getfemlistargyris.fig                   \
     getfemlistcubeQ1.fig                    \
     getfemlistcubeQ3.fig                    \
     getfemlistHCT.fig                       \
     getfemlistincomplete.fig                \
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
     getfemlistsymbols00.fig                 \
     getfemlistsymbols01.fig                 \
     getfemlistsymbols02.fig                 \
     getfemlistsymbols03.fig                 \
     getfemlistsymbols04.fig                 \
     getfemlistsymbols05.fig                 \
     getfemlistsymbols06.fig                 \
     getfemlistsymbols07.fig                 \
     getfemlistsymbols08.fig                 \
     getfemlistsymbols09.fig                 \
     getfemlistsymbols10.fig                 \
     getfemlistsymbols11.fig                 \
     getfemlistsymbols12.fig                 \
     getfemlistsymbols13.fig                 \
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
     getfemuserrefine.fig		     \
     getfemusermodelmasterslave.fig          \
     getfemusermodeldetectcontact.fig        \
     getfemusermodelfalsecontact1.fig        \
     getfemusermodelfalsecontact2.fig        \
     getfemlistpyramidP0.fig                 \
     getfemlistpyramidP1.fig                 \
     getfemlistpyramidP2.fig                 \
     ALE_rotating_body.fig		     \
     ALE_translation_body.fig		     \
     ALE_rotating_conf.fig

EPSFIGS=$(FIGS:.fig=.eps)
PNGFIGS=$(FIGS:.fig=.png)

.SUFFIXES: .fig .eps .png

.fig.eps:
	../../../../../bin/fig2eps $(@:.eps=.fig)

.eps.png:
	convert -density 200 $(@:.png=.eps) $@

png: $(PNGFIGS)

clean:
	-rm -rf *.aux
	-rm -rf *.dvi
	-rm -rf *.log
	-rm -rf *.bak
	-rm -rf $(EPSFIGS)
	-rm -rf $(PNGFIGS)
