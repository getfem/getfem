#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Python GetFEM interface
#
# Copyright (C) 2019-2020 Yves Renard.
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
"""  Discretization of a Poisson problem using HHO methods.

  This program is used to check that python-getfem is working. This is
  also a good example of use of GetFEM.

  $Id$
"""
import numpy as np

# Import basic modules
import getfem as gf

## Parameters
NX = 20                            # Mesh parameter.
N = 2
Dirichlet_with_multipliers = True  # Dirichlet condition with multipliers
                                   # or penalization
dirichlet_coefficient = 1e10       # Penalization coefficient
using_HHO = True                   # Use HHO method or standard Lagrange FEM

# Create a simple cartesian mesh
I = np.arange(0,1+1./NX,1./NX)
if (N == 2):
  m = gf.Mesh('regular_simplices', I, I)
elif (N == 3):
  m = gf.Mesh('regular_simplices', I, I, I)
  
  

# Meshfems
mfu   = gf.MeshFem(m, 1)
mfgu  = gf.MeshFem(m, N)
mfur  = gf.MeshFem(m, 1)
mfrhs = gf.MeshFem(m, 1)

if (using_HHO):
  mfu.set_fem(gf.Fem('FEM_HHO(FEM_SIMPLEX_IPK(%d,2),FEM_SIMPLEX_CIPK(%d,2))' % (N,N-1)))
  mfur.set_fem(gf.Fem('FEM_PK(%d,3)' % N))
else:
  mfu.set_fem(gf.Fem('FEM_PK(%d,2)' % N))
  mfur.set_fem(gf.Fem('FEM_PK(%d,2)' % N))

mfgu.set_fem(gf.Fem('FEM_PK(%d,2)' % N))
mfrhs.set_fem(gf.Fem('FEM_PK(%d,2)' % N))

print('nbdof : %d' % mfu.nbdof());

#  Integration method used
if (N == 2):
  mim = gf.MeshIm(m, gf.Integ('IM_TRIANGLE(4)'))
else:
  mim = gf.MeshIm(m, gf.Integ('IM_TETRAHEDRON(5)'))

# Boundary selection
flst  = m.outer_faces()
fnor  = m.normal_of_faces(flst)
tleft = abs(fnor[1,:]+1) < 1e-14
ttop  = abs(fnor[0,:]-1) < 1e-14
fleft = np.compress(tleft, flst, axis=1)
ftop  = np.compress(ttop, flst, axis=1)
fneum = np.compress(np.logical_not(ttop + tleft), flst, axis=1)

# Mark it as boundary
DIRICHLET_BOUNDARY_NUM1 = 1
DIRICHLET_BOUNDARY_NUM2 = 2
NEUMANN_BOUNDARY_NUM = 3
m.set_region(DIRICHLET_BOUNDARY_NUM1, fleft)
m.set_region(DIRICHLET_BOUNDARY_NUM2, ftop)
m.set_region(NEUMANN_BOUNDARY_NUM, fneum)

# Faces for stabilization term
all_faces = m.all_faces()
ALL_FACES = 4
m.set_region(ALL_FACES, all_faces)

# Interpolate the exact solution (Assuming mfu is a Lagrange fem)
Ue = mfur.eval('y*(y-1)*x*(x-1)+x*x*x*x*x')

# Interpolate the source term
F1 = mfrhs.eval('-(2*(x*x+y*y)-2*x-2*y+20*x*x*x)')
if (N == 2):
  F2 = mfrhs.eval('[y*(y-1)*(2*x-1) + 5*x*x*x*x, x*(x-1)*(2*y-1)]')
else:
  F2 = mfrhs.eval('[y*(y-1)*(2*x-1) + 5*x*x*x*x, x*(x-1)*(2*y-1), 0]')

# Model
md = gf.Model('real')



# Main unknown
md.add_fem_variable('u', mfu)
md.add_fem_data('Gu', mfgu)
md.add_fem_data('ur', mfur)

# Needed reconstuction and stabilization operators
if (using_HHO):
  md.add_HHO_reconstructed_gradient('HHO_Grad')
  md.add_HHO_reconstructed_value('HHO_Val')
  md.add_HHO_stabilization('HHO_Stab')
  md.add_macro('HHO_Val_u', 'Elementary_transformation(u, HHO_Val, ur)')
  md.add_macro('HHO_Grad_u', 'Elementary_transformation(u, HHO_Grad, Gu)')
  md.add_macro('HHO_Grad_Test_u',
               'Elementary_transformation(Test_u, HHO_Grad, Gu)')
  md.add_macro('HHO_Stab_u', 'Elementary_transformation(u, HHO_Stab)')
  md.add_macro('HHO_Stab_Test_u',
               'Elementary_transformation(Test_u, HHO_Stab)')


# Laplacian term on u
if (using_HHO):
  # Laplacian term
  md.add_linear_term(mim, 'HHO_Grad_u.HHO_Grad_Test_u')
  # Stabilization term
  md.add_linear_term(mim, '10*HHO_Stab_u.HHO_Stab_Test_u', ALL_FACES)
else:
  md.add_Laplacian_brick(mim, 'u')

# Volumic source term
md.add_initialized_fem_data('VolumicData', mfrhs, F1)
md.add_source_term_brick(mim, 'u', 'VolumicData')

# Neumann condition.
md.add_initialized_fem_data('NeumannData', mfrhs, F2)
md.add_normal_source_term_brick(mim, 'u', 'NeumannData', NEUMANN_BOUNDARY_NUM)

# Dirichlet condition on the left.
md.add_initialized_fem_data("Ue", mfur, Ue)

if (Dirichlet_with_multipliers):
  md.add_Dirichlet_condition_with_multipliers(mim, 'u', mfu,
                                              DIRICHLET_BOUNDARY_NUM1, 'Ue')
else:
  md.add_Dirichlet_condition_with_penalization(mim, 'u', dirichlet_coefficient,
                                               DIRICHLET_BOUNDARY_NUM1, 'Ue')

# Dirichlet condition on the top.
# Two Dirichlet brick in order to test the multiplier
# selection in the intersection.
if (Dirichlet_with_multipliers):
  md.add_Dirichlet_condition_with_multipliers(mim, 'u', mfu,
                                              DIRICHLET_BOUNDARY_NUM2, 'Ue')
else:
  md.add_Dirichlet_condition_with_penalization(mim, 'u', dirichlet_coefficient,
                                               DIRICHLET_BOUNDARY_NUM2, 'Ue')

# Assembly of the linear system and solve.
md.solve()

# Error computation
U = md.variable('u')
L2error = gf.asm('generic', mim, 0, 'sqr(u-Ue)', -1, md)
H1error = gf.asm('generic', mim, 0, 'Norm_sqr(Grad_u-Grad_Ue)', -1, md)
H1error = np.sqrt(L2error + H1error); L2error = np.sqrt(L2error)
print('Error in L2 norm (without reconstruction): %g' % L2error)
print('Error in H1 norm (without reconstruction): %g' % H1error)
if (using_HHO):
  L2error = gf.asm('generic', mim, 0, 'sqr(HHO_Val_u-Ue)', -1, md)
  H1error = gf.asm('generic', mim, 0, 'Norm_sqr(HHO_Grad_u-Grad_Ue)', -1, md)
  H1error = np.sqrt(L2error + H1error); L2error = np.sqrt(L2error)
  print('Error in L2 norm (with reconstruction): %g' % L2error)
  print('Error in H1 norm (with reconstruction): %g' % H1error)


# Export data
# mfur.export_to_pos('laplacian_e.pos', Ue, 'Exact solution')
mfu.export_to_pos('laplacian.pos', U, 'Computed solution')
print('You can view the solution with (for example):')
print('gmsh laplacian.pos')


if (H1error > 3e-5):
    print('Error too large !')
    exit(1)
