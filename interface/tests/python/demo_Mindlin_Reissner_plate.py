#!/usr/bin/env python
# -*- coding: UTF8 -*-
# Python GetFEM interface
#
# Copyright (C) 2015-2020 FABRE Mathieu, SECK Mamadou, DALLERIT Valentin,
#                         Yves Renard.
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
"""  Simple supported Mindlin-Reissner plate demonstration.

  This program is used to check that python-getfem is working. This is
  also a good example of use of GetFEM.

"""
import numpy as np

import getfem as gf

## Parameters
Emodulus = 1.          # Young Modulus
nu       = 0.5         # Poisson Coefficient
epsilon  = 0.001       # Plate thickness
kappa     = 5./6.      # Shear correction factor
f = -5.*pow(epsilon,3.) # Prescribed force on the top of the plate

variant = 0            # 0 : not reduced, 1 : with reduced integration,
                       # 2 : MITC reduction
quadrangles = True     # Locking free only on quadrangle for the moment
K = 1                  # Degree of the finite element method
dirichlet_version = 1  # 0 = simplification, 1 = with multipliers,
                       # 2 = penalization
r = 1.E8               # Penalization parameter.
NX = 80                # Number of element per direction


if (quadrangles):
  m = gf.Mesh('cartesian', np.arange(0., 1.+1./NX, 1./NX),np.arange(0., 1.+1./NX, 1./NX))
else:
  m = gf.Mesh('import','structured',
              'GT="GT_PK(2,1)";SIZES=[1,1];NOISED=0;NSUBDIV=[%d,%d];'
              % (NX, NX))

## Create a mesh_fem for a 2D vector field
mftheta = gf.MeshFem(m, 2)
mfu = gf.MeshFem(m, 1);
mftheta.set_classical_fem(K)
mfu.set_classical_fem(K)
mim = gf.MeshIm(m, 6)
mim_reduced = gf.MeshIm(m, 1)

## Detect the border of the mesh and assign it the boundary number 1
border = m.outer_faces()
m.set_region(1, border)

## Build the model
md=gf.Model('real')
md.add_fem_variable('u', mfu)
md.add_fem_variable('theta', mftheta)
md.add_initialized_data('E', Emodulus)
md.add_initialized_data('nu', nu)
md.add_initialized_data('epsilon', epsilon)
md.add_initialized_data('kappa', kappa)

md.add_Mindlin_Reissner_plate_brick(mim, mim_reduced, 'u', 'theta', 'E', 'nu', 'epsilon', 'kappa', variant)

md.add_initialized_data('VolumicData', f)
md.add_source_term_brick(mim, 'u', 'VolumicData');
md.add_initialized_data('DirichletData', 0);

if (dirichlet_version == 0):
    md.add_Dirichlet_condition_with_simplification('u', 1, 'DirichletData');   
elif (dirichlet_version == 1):
    md.add_Dirichlet_condition_with_multipliers(mim, 'u', mfu, 1, 'DirichletData');
elif (dirichlet_version == 2):
    md.add_Dirichlet_condition_with_penalization(mim, 'u', r, 1, 'DirichletData');

print ('Number of Dof: %d' % md.nbdof()) 

md.solve()
U = md.variable('u');

mfu.export_to_vtk('Deflection.vtk', U)
print ('You can view solutions with for instance:\nmayavi2 -d Deflection.vtk -f WarpScalar -m Surface')
