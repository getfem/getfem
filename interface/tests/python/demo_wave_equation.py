#!/usr/bin/env python
# -*- coding: UTF8 -*-
# Python GetFEM++ interface
#
# Copyright (C) 2015-2015 FABRE Mathieu, SECK Mamadou, DALLERIT Valentin,
#                         Yves Renard.
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
"""  Simple demo of a wave equation solved with a Newmark scheme with the
     Getfem tool for time integration schemes

  This program is used to check that python-getfem is working. This is
  also a good example of use of GetFEM++.

"""
import numpy as np
import getfem as gf
import os

NX = 20
m = gf.Mesh('cartesian', np.arange(0., 1.+1./NX,1./NX),
            np.arange(0., 1.+1./NX,1./NX));


## create a mesh_fem of for a field of dimension 1 (i.e. a scalar field)
mf = gf.MeshFem(m,1)
mf.set_classical_fem(2)

## Integration which will be used
mim = gf.MeshIm(m, 4)


## Detect the border of the mesh
border = m.outer_faces()
m.set_region(1, border)

## Interpolate the initial data
U0 = mf.eval('y*(y-1.)*x*(x-1.)*x*x')
V0 = 0.*U0

md=gf.Model('real');
md.add_fem_variable('u', mf);
md.add_Laplacian_brick(mim, 'u');
md.add_Dirichlet_condition_with_multipliers(mim, 'u', mf, 1);
# md.add_Dirichlet_condition_with_penalization(mim, 'u', 1E9, 1);
# md.add_Dirichlet_condition_with_simplification('u', 1);

## Transient part.
T = 5.0;
dt = 0.025;
beta = 0.25;
gamma = 0.5;

md.add_Newmark_scheme('u', beta, gamma)
md.add_mass_brick(mim, 'Dot2_u')
md.set_time_step(dt)

## Initial data.
md.set_variable('Previous_u',  U0)
md.set_variable('Previous_Dot_u',  V0)

## Initialisation of the acceleration 'Previous_Dot2_u'
md.perform_init_time_derivative(dt/2.)
md.solve()

A0 = md.variable('Previous_Dot2_u')


os.system('mkdir results');
mf.export_to_vtk('results/displacement_0.vtk', U0)
mf.export_to_vtk('results/velocity_0.vtk', V0)
mf.export_to_vtk('results/acceleration_0.vtk', A0)

## Iterations
n = 1;
for t in np.arange(0.,T,dt):
  print ('Time step %g' % t)
  md.solve();
  U = md.variable('u')
  V = md.variable('Dot_u')
  A = md.variable('Dot2_u')
  
  mf.export_to_vtk('results/displacement_%d.vtk' % n, U)
  mf.export_to_vtk('results/velocity_%d.vtk' % n, V)
  mf.export_to_vtk('results/acceleration_%d.vtk' % n, A)

  n += 1
  md.shift_variables_for_time_integration()
  

print ('You can view solutions with for instance:\nmayavi2 -d results/displacement_0.vtk -f WarpScalar -m Surface')
