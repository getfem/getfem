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
"""  2D elasticity problem using HHO methods.

  This program is used to check that python-getfem is working. This is
  also a good example of use of GetFEM.

  $Id$
"""
import numpy as np

# Import basic modules
import getfem as gf

## Parameters
NX = 40                           # Mesh parameter.
Dirichlet_with_multipliers = True # Dirichlet condition with multipliers
                                  # or penalization
dirichlet_coefficient = 1e10      # Penalization coefficient
using_HHO = True                  # Use HHO method or standard Lagrange FEM
using_symmetric_gradient = True   # Use symmetric gradient reconstruction or not

E = 1                             # Young's modulus
nu = 0.3                          # Poisson ratio

cmu = E/(2*(1+nu))                # Lame coefficient
clambda = 2*cmu*nu/(1-2*nu)       # Lame coefficient
use_quad = False                  # Quadrilaterals or triangles

# Create a simple cartesian mesh
I = np.arange(0,1+1./NX,1./NX)
if (use_quad):
  m = gf.Mesh('cartesian', I, I)
  # m=gf.Mesh('import','structured','GT="GT_QK(2,1)";SIZES=[1,1];NOISED=1;NSUBDIV=[1,1];')
else:
  m = gf.Mesh('regular_simplices', I, I)
  # m=gf.Mesh('import','structured','GT="GT_PK(2,1)";SIZES=[1,1];NOISED=1;NSUBDIV=[%d,%d];' % (NX, NX))
  

N = m.dim();

# Meshfems
mfu   = gf.MeshFem(m, N)
mfgu  = gf.MeshFem(m, N, N)
mfur  = gf.MeshFem(m, N)
mfrhs = gf.MeshFem(m, 1)

if (using_HHO):
  if (use_quad):
    mfu.set_fem(gf.Fem('FEM_HHO(FEM_QUAD_IPK(2,2),FEM_SIMPLEX_CIPK(1,2))'))
    # mfu.set_fem(gf.Fem('FEM_HHO(FEM_QK_DISCONTINUOUS(2,2,0.1),FEM_SIMPLEX_CIPK(1,2))'))
    # mfu.set_fem(gf.Fem('FEM_HHO(FEM_QK(2,2),FEM_PK(1,2))'))
    mfur.set_fem(gf.Fem('FEM_QUAD_IPK(2,3)'))
    # mfur.set_fem(gf.Fem('FEM_QK(2,3)'))
  else:
    mfu.set_fem(gf.Fem('FEM_HHO(FEM_SIMPLEX_IPK(2,2),FEM_SIMPLEX_CIPK(1,2))'))
    mfur.set_fem(gf.Fem('FEM_PK(2,3)'))
else:
  if (use_quad):
    mfu.set_fem(gf.Fem('FEM_QK(2,2)'))
    mfur.set_fem(gf.Fem('FEM_QK(2,2)'))
  else:
    mfu.set_fem(gf.Fem('FEM_PK(2,2)'))
    mfur.set_fem(gf.Fem('FEM_PK(2,2)'))
    
if (use_quad):
  mfgu.set_fem(gf.Fem('FEM_QUAD_IPK(2,2)'))
  # mfgu.set_fem(gf.Fem('FEM_QK(2,2)'))
  mfrhs.set_fem(gf.Fem('FEM_QK(2,3)'))
else:
  mfgu.set_fem(gf.Fem('FEM_PK(2,2)'))
  mfrhs.set_fem(gf.Fem('FEM_PK(2,3)'))

print('nbdof : %d' % mfu.nbdof());

#  Integration method used
if (use_quad):
  mim = gf.MeshIm(m, gf.Integ('IM_GAUSS_PARALLELEPIPED(2,6)'))
else:
  mim = gf.MeshIm(m, gf.Integ('IM_TRIANGLE(6)'))

# Boundary selection
flst  = m.outer_faces()
GAMMAD = 1
m.set_region(GAMMAD, flst)

# Faces for stabilization term
all_faces = m.all_faces()
ALL_FACES = 4
m.set_region(ALL_FACES, all_faces)

# Interpolate the exact solution (Assuming mfu is a Lagrange fem)
a = 8.
Ue = mfrhs.eval('[np.sin((%g)*x), -(%g)*y*np.cos((%g)*x)]' % (a,a,a), globals(), locals())

# Interpolate the source term
F1 = mfrhs.eval('(%g)*(pow(%g,2.))*np.sin((%g)*x), -(%g)*(pow(%g,3.))*y*np.cos((%g)*x)' % (cmu, a, a, cmu, a, a), globals(), locals())



# Model
md = gf.Model('real')



# Main unknown
md.add_fem_variable('u', mfu)
md.add_fem_data('Gu', mfgu)
md.add_fem_data('ur', mfur)

# Needed reconstuction and stabilization operators
if (using_HHO):
  if (using_symmetric_gradient):
    md.add_HHO_reconstructed_symmetrized_gradient('HHO_Grad')
    md.add_HHO_reconstructed_gradient('HHO_vGrad')
    md.add_HHO_reconstructed_symmetrized_value('HHO_Val')
    md.add_HHO_symmetrized_stabilization('HHO_Stab')
    md.add_macro('HHO_vGrad_u', 'Elementary_transformation(u, HHO_vGrad, Gu)')
  else :
    md.add_HHO_reconstructed_gradient('HHO_Grad')
    md.add_HHO_reconstructed_value('HHO_Val')
    md.add_HHO_stabilization('HHO_Stab')
    md.add_macro('HHO_vGrad_u', 'Elementary_transformation(u, HHO_Grad, Gu)')

    
  md.add_macro('HHO_Val_u', 'Elementary_transformation(u, HHO_Val, ur)')
  md.add_macro('HHO_Grad_u', 'Elementary_transformation(u, HHO_Grad, Gu)')
  md.add_macro('HHO_Grad_Test_u',
               'Elementary_transformation(Test_u, HHO_Grad, Gu)')
  md.add_macro('HHO_Stab_u', 'Elementary_transformation(u, HHO_Stab)')
  md.add_macro('HHO_Stab_Test_u',
               'Elementary_transformation(Test_u, HHO_Stab)')


# Elasticity term on u
md.add_initialized_data('cmu', [cmu])
md.add_initialized_data('clambda', [clambda])
if (using_HHO):
  # Elasticity term
  md.add_linear_term(mim, 'clambda*Trace(HHO_Grad_u)*Trace(HHO_Grad_Test_u)'
                     +    '+ 2*cmu*Sym(HHO_Grad_u):Sym(HHO_Grad_Test_u)')
  # Stabilization term
  md.add_linear_term(mim, '20*cmu*HHO_Stab_u.HHO_Stab_Test_u', ALL_FACES)
else:
  md.add_isotropic_linearized_elasticity_brick(mim, 'u', 'clambda', 'cmu')

# Volumic source term
md.add_initialized_fem_data('VolumicData', mfrhs, F1)
md.add_source_term_brick(mim, 'u', 'VolumicData')

# Dirichlet condition
md.add_initialized_fem_data("Ue", mfrhs, Ue)

if (Dirichlet_with_multipliers):
  md.add_Dirichlet_condition_with_multipliers(mim, 'u', mfu, GAMMAD, 'Ue')
else:
  md.add_Dirichlet_condition_with_penalization(mim, 'u', dirichlet_coefficient,
                                               GAMMAD, 'Ue')

# Assembly of the linear system and solve.
md.solve()

# Error computation
U = md.variable('u')
L2error = gf.asm('generic', mim, 0, 'Norm_sqr(u-Ue)', -1, md)
H1error = gf.asm('generic', mim, 0, 'Norm_sqr(Grad_u-Grad_Ue)', -1, md)
H1error = np.sqrt(L2error + H1error); L2error = np.sqrt(L2error)
print('Error in L2 norm (without reconstruction): %g' % L2error)
print('Error in H1 norm (without reconstruction): %g' % H1error)

if (using_HHO):
  L2error = gf.asm('generic', mim, 0, 'Norm_sqr(HHO_Val_u-Ue)', -1, md)
  H1error = gf.asm('generic', mim, 0, 'Norm_sqr(HHO_vGrad_u-Grad_Ue)', -1, md)
  H1error = np.sqrt(L2error + H1error); L2error = np.sqrt(L2error)
  print('Error in L2 norm (with reconstruction): %g' % L2error)
  print('Error in H1 norm (with reconstruction): %g' % H1error)


# Export data
# mfur.export_to_pos('elasticity_e.pos', Ue, 'Exact solution')
mfu.export_to_pos('elasticity.pos', U, 'Computed solution')
print('You can view the solution with (for example):')
print('gmsh elasticity.pos')


if (H1error > 0.011):
    print('Error too large !')
    exit(1)
