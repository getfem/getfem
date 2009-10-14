#!/usr/bin/env python
# -*- coding: UTF8 -*-
# Python GetFEM++ interface
#
# Copyright (C) 2004-2009 Yves Renard, Julien Pommier.
#
# This file is a part of GetFEM++
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
############################################################################
"""  Plate problem test.

  This program is used to check that python-getfem is working. This is
  also a good example of use of GetFEM++.

  $Id$
"""
from getfem import *
from numpy import *

NX=10.0
mixed = True
thickness = 0.01;

m=Mesh('regular simplices', arange(0,1.01,1/NX), arange(0,1.01,1/NX))
mfut = MeshFem(m,2)
mfu3 = MeshFem(m,1)
mfth = MeshFem(m,2)
mfd  = MeshFem(m,1)

mfut.set_fem(Fem('FEM_PK(2,2)'))
mfu3.set_fem(Fem('FEM_PK(2,1)'))
mfth.set_fem(Fem('FEM_PK(2,2)'))
mfd.set_fem(Fem('FEM_PK(2,2)'))

mim = MeshIm(m, Integ('IM_TRIANGLE(5)'))


#get the list of faces whose normal is [-1,0]
flst = m.outer_faces();
fnor = m.normal_of_faces(flst);
fleft = compress(abs(fnor[1,:]+1) < 1e-14, flst, axis=1);
fright= compress(abs(fnor[1,:]-1) < 1e-14, flst, axis=1);
CLAMPED_BOUNDARY = 1;
m.set_region(CLAMPED_BOUNDARY, fleft);
SIMPLE_SUPPORT_BOUNDARY = 2
m.set_region(SIMPLE_SUPPORT_BOUNDARY, fright);

E=1e3
Nu=0.3
Lambda = E*Nu/((1+Nu)*(1-2*Nu))
Mu =E/(2*(1+Nu))


if not mixed:
    b0 = MdBrick('isotropic_linearized_plate',mim,mim,mfut,mfu3,mfth,thickness)
else:
    b0 = MdBrick('mixed_isotropic_linearized_plate',mim,mfut,mfu3,mfth,thickness)

b1 = MdBrick('plate_source_term', b0)
b1.set('param', 'M', mfd, mfd.eval('[0, x[1]*x[1]/1000]'))

b2 = MdBrick('plate clamped support', b1, CLAMPED_BOUNDARY, 'augmented');

b3 = MdBrick('plate simple support', b2, SIMPLE_SUPPORT_BOUNDARY, 'augmented');

b4 = b3
if mixed:
    b4 = MdBrick('plate closing', b3)

mds=MdState(b4)
print 'running solve...'
b4.solve(mds, 'noisy', 'lsolver','superlu')
print 'solve done!'


U=mds.state()
nut = mfut.nbdof()
nu3 = mfu3.nbdof()
nth = mfth.nbdof()
ut=U[0:nut]
u3=U[nut:(nut+nu3)]
th=U[(nut+nu3):(nut+nu3+nth)]


sl=Slice(('none',), mfu3, 4)
sl.export_to_vtk('plate.vtk', mfu3, u3, 'Displacement')
sl.export_to_pos('plate.pos', mfu3, u3,'Displacement')

print 'You can view the solution with (for example):'
print 'mayavi -d ./plate.vtk -f WarpScalar -m BandedSurfaceMap'
print 'or'
print 'mayavi2 -d plate.vtk -f WarpScalar -m Surface'
print 'or'
print 'gmsh plate.pos'
