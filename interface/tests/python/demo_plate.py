#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Python GetFEM interface
#
# Copyright (C) 2004-2020 Yves Renard, Julien Pommier.
#
# This file is a part of GetFEM
#
# GetFEM  is  free software;  you  can  redistribute  it  and/or modify it
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
  also a good example of use of GetFEM.

  $Id$
"""
import numpy as np

import getfem as gf

NX=10.0
thickness = 0.01;
f = -5.*pow(thickness,3.);

m=gf.Mesh('regular simplices', np.arange(0,1.01,1/NX), np.arange(0,1.01,1/NX))
mfu3 = gf.MeshFem(m,1)
mfth = gf.MeshFem(m,2)
mfd  = gf.MeshFem(m,1)

mfu3.set_fem(gf.Fem('FEM_PK(2,1)'))
mfth.set_fem(gf.Fem('FEM_PK(2,2)'))
mfd.set_fem(gf.Fem('FEM_PK(2,2)'))

mim = gf.MeshIm(m, gf.Integ('IM_TRIANGLE(5)'))


#get the list of faces whose normal is [-1,0]
flst = m.outer_faces();
fnor = m.normal_of_faces(flst);
fleft = np.compress(abs(fnor[1,:]+1) < 1e-14, flst, axis=1);
fright= np.compress(abs(fnor[1,:]-1) < 1e-14, flst, axis=1);
CLAMPED_BOUNDARY = 1;
m.set_region(CLAMPED_BOUNDARY, fleft);
SIMPLE_SUPPORT_BOUNDARY = 2
m.set_region(SIMPLE_SUPPORT_BOUNDARY, fright);

E=1e3
Nu=0.3


md = gf.Model('real')
md.add_fem_variable('u3', mfu3)
md.add_fem_variable('theta', mfth)
md.add_initialized_data('E', E)
md.add_initialized_data('nu', Nu)
md.add_initialized_data('epsilon', thickness)
md.add_initialized_data('kappa', 5./6.)
md.add_Mindlin_Reissner_plate_brick(mim, mim, 'u3', 'theta', 'E', 'nu', 'epsilon', 'kappa', 2)
md.add_initialized_data('VolumicData', f)
md.add_source_term_brick(mim, 'u3', 'VolumicData')

md.add_Dirichlet_condition_with_multipliers(mim, 'u3', mfu3, CLAMPED_BOUNDARY);
md.add_Dirichlet_condition_with_multipliers(mim, 'theta', mfth, CLAMPED_BOUNDARY);
md.add_Dirichlet_condition_with_multipliers(mim, 'u3', mfu3, SIMPLE_SUPPORT_BOUNDARY);
                                              
                                             


print('running solve...')
md.solve()
print('solve done!')


u3 = md.variable('u3')


sl=gf.Slice(('none',), mfu3, 4)
sl.export_to_vtk('plate.vtk', mfu3, u3, 'Displacement')
sl.export_to_pos('plate.pos', mfu3, u3, 'Displacement')

print('You can view the solution with (for example):')
print('mayavi2 -d plate.vtk -f WarpScalar -m Surface')
print('or')
print('gmsh plate.pos')
