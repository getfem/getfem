# Python GetFEM++ interface
#
# Copyright (C) 2004-2009 Yves Renard, Julien Pommier.
#                                                       
# This file is a part of GETFEM++                                         
#                                                                         
# GetFEM++  is  free software;  you  can  redistribute  it  and/or modify it
# under  the  terms  of the  GNU  Lesser General Public License as published
# by  the  Free Software Foundation;  either version 2.1 of the License,  or
# (at your option) any later version.
# This program  is  distributed  in  the  hope  that it will be useful,  but
# WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
# or  FITNESS  FOR  A PARTICULAR PURPOSE.  See the GNU Lesser General Public
# License for more details.
# You  should  have received a copy of the GNU Lesser General Public License
# along  with  this program;  if not, write to the Free Software Foundation,
# Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301, USA.
#

from getfem import *
from numpy import *

m0=Mesh('cartesian',[0,1,2,3],[0,1,2],[-3,-2])

m0.add_convex(GeoTrans('GT_QK(2,2)'),[[0,0,0,.4,.6,.5,1.2,1,1],
                                     [0,.3,1,0,.5,1,-.1,.5,1.1],
                                     [0,0,0,0,0,0,0,0,0]])

m0.add_convex(GeoTrans('GT_PK(2,2)'),[[2,2,2,2.6,2.5,3.2],
                                     [0,.3,1,0,.5,0], [0,0,0,0,0,0]])

m0.add_convex(GeoTrans('GT_PK(1,3)'),[[3.1, 2.8, 3.2, 3.7],
                                     [0, .3, .7, 1.3],
                                     [0,0,0,0]])

m0.add_convex(GeoTrans('GT_PRISM(3,1)'), [[0, 1, 0, 0, 1, 0],
                                         [0, 0, 1, 0, 0, 1],
                                         [0, 0, 0, 1, 1, 1]]);


m0.add_convex(GeoTrans('GT_PK(3,2)'), [array([0, .5, 1, 0, .5, 0, 0, .5, 0, 0])-1.5,
                                      array([0, 0, 0, .5, .5, 1, 0, 0, .5, 0])-1,
                                      array([0, 0, 0, 0, 0, 0, .5, .5, .5, 1])-1.1])

m0.add_convex(GeoTrans('GT_QK(3,2)'), [array([0, .5, 1, 0, .5, 1, 0, .5, 1, 0, .5, 1, 0, .5, 1, 0, .5, 1, 0, .5, 1, 0, .5, 1, 0, .5, 1])-1.5,
                                      array([0, 0, 0, .5, .5, .5, 1, 1, 1, 0, 0, 0, .5, .5, .5, 1, 1, 1, 0, 0, 0, .5, .5, .5, 1, 1, 1])-1,
                                      array([0, 0, 0, 0, 0, 0, 0, 0, 0, .5, .5, .5, .5, .5, .5, .5, .5, .5, 1, 1, 1, 1, 1, 1, 1, 1, 1])])

m1=Mesh('cartesian',[0,1,2,3],[0,1,2],[-3,-2])

mf0=MeshFem(m0); mf0.set_classical_fem(1); mf0.export_to_vtk('check_export0.vtk','ascii')
mf1=MeshFem(m1); mf1.set_classical_fem(1); mf1.export_to_vtk('check_export1.vtk','ascii')

m0.export_to_vtk('check_export2.vtk','quality')
m1.export_to_vtk('check_export3.vtk','quality')

try:
    m0.export_to_dx('check_export0.dx')
except RuntimeError:
    pass
m1.export_to_dx('check_export0.dx','ascii','edges')
m1.export_to_dx('check_export0.dx','ascii','append')

sl=Slice(('boundary',),m0,6)
sl.export_to_dx('check_export1.dx','ascii')
sl.export_to_dx('check_export1.dx','append','edges')

U=random.standard_normal(mf0.nbdof())
sl.export_to_dx('check_export1.dx','append','serie','rndserie',mf0,[U])
sl.export_to_dx('check_export1.dx','append','serie','rndserie',mf0,[U])
