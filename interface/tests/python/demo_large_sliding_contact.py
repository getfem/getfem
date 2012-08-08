#!/usr/bin/env python
# -*- coding: UTF8 -*-
# Python GetFEM++ interface
#
# Copyright (C) 2012-2012 Yves Renard.
#
# This file is a part of GetFEM++
#
# Getfem++  is  free software;  you  can  redistribute  it  and/or modify it
# under  the  terms  of the  GNU  Lesser General Public License as published
# by  the  Free Software Foundation;  either version 3 of the License,  or
# (at your option) any later version along with the GCC Runtime Library
# Exception either version 3.1 or (at your option) any later version.
# This program  is  distributed  in  the  hope  that it will be useful,  but
# WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
# or  FITNESS  FOR  A PARTICULAR PURPOSE.  See the GNU Lesser General Public
# License and GCC Runtime Library Exception for more details.
# You  should  have received a copy of the GNU Lesser General Public License
# along  with  this program;  if not, write to the Free Software Foundation,
# Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301, USA.
#
############################################################################

import getfem as gf
import numpy as np

clambda = 1.    # Elasticity parameters
cmu = 1.
r = 1.0         # Augmentation parameter
f_coeff = 1.    # Friction coefficient
vf = 0.01       # Vertical force
penalty_parameter = 0.1

mesh1 = gf.Mesh('load', '../../../tests/meshes/disc_with_a_hole.mesh')
mesh2 = gf.Mesh('import', 'structured', 'GT="GT_PK(2,1)";ORG=[-0.5,0];SIZES=[1,0.1];NSUBDIV=[20,2]')

N = mesh1.dim()

mfu1 = gf.MeshFem(mesh1, N)
mfu1.set_classical_fem(2)

pre_mflambda1 = gf.MeshFem(mesh1, N)
pre_mflambda1.set_classical_fem(1)

mfvm1 = gf.MeshFem(mesh1)
mfvm1.set_classical_discontinuous_fem(1)

fb1 = mesh1.outer_faces()
CONTACT_BOUNDARY1 = 1
mesh1.set_region(CONTACT_BOUNDARY1, fb1)

dol1 = pre_mflambda1.basic_dof_on_region(CONTACT_BOUNDARY1)
mflambda1 = gf.MeshFem('partial', pre_mflambda1, dol1)

mim1 = gf.MeshIm(mesh1, 4)


mfu2 = gf.MeshFem(mesh2, N)
mfu2.set_classical_fem(1)

pre_mflambda2 = gf.MeshFem(mesh2, N)
pre_mflambda2.set_classical_fem(1)

mfvm2 = gf.MeshFem(mesh2)
mfvm2.set_classical_discontinuous_fem(1)

fb2 = mesh2.outer_faces()
CONTACT_BOUNDARY2 = 2
mesh2.set_region(CONTACT_BOUNDARY2, fb2)

dol2 = pre_mflambda2.basic_dof_on_region(CONTACT_BOUNDARY2)
mflambda2 = gf.MeshFem('partial',  pre_mflambda2, dol2)

mim2 = gf.MeshIm(mesh2, 8)


two_bodies = 1

md = gf.Model('real')
md.add_initialized_data('lambda', clambda)
md.add_initialized_data('mu', cmu)

if two_bodies:
   md.add_fem_variable('u1', mfu1)
   md.add_fem_variable('lambda1', mflambda1)
   md.add_isotropic_linearized_elasticity_brick(mim1, 'u1', 'lambda', 'mu')
#   md.add_initialized_data('cpoints1', [0 0.5 0 1.5 0 0.5 0 1.5])
#   md.add_initialized_data('cunitv1', [1 0 1 0 0 1 0 1])
#   md.add_initialized_data('cdata', [0 0 -0.01 -0.01])
#   md.add_pointwise_constraints_with_multipliers('u1', 'cpoints1', 'cunitv1', 'cdata')
   md.add_initialized_data('cpoints1', [0,0.5,0,1.5])
   md.add_initialized_data('cunitv1', [1,0,1,0])
   md.add_initialized_data('cdata', [0,0])
   md.add_pointwise_constraints_with_multipliers('u1', 'cpoints1', 'cunitv1', 'cdata')
   md.add_initialized_data('data1', [0,-vf])
   md.add_source_term_brick(mim1, 'u1', 'data1')
   md.add_initialized_data('penalty_param1', [penalty_parameter])
   md.add_mass_brick(mim1, 'u1', 'penalty_param1')

md.add_fem_variable('u2', mfu2)
md.add_fem_variable('lambda2', mflambda2)
md.add_isotropic_linearized_elasticity_brick(mim2, 'u2', 'lambda', 'mu')
md.add_initialized_data('cpoints2', [0,0])
md.add_initialized_data('cunitv2', [1,0])
md.add_pointwise_constraints_with_multipliers('u2', 'cpoints2', 'cunitv2');
md.add_initialized_data('data2', [0,-vf])
md.add_source_term_brick(mim2, 'u2', 'data2')
md.add_initialized_data('penalty_param2', [penalty_parameter])
md.add_mass_brick(mim2, 'u2', 'penalty_param2')

md.add_initialized_data('r', r)
md.add_initialized_data('f', f_coeff)

indb = md.add_integral_large_sliding_contact_brick(mim2, 'u2', 'lambda2', 'r', 'f', CONTACT_BOUNDARY2)

if two_bodies:
   md.add_boundary_to_large_sliding_contact_brick(indb, mim1, 'u1', 'lambda1', CONTACT_BOUNDARY1)

md.add_rigid_obstacle_to_large_sliding_contact_brick(indb, 'y')

# md.test_tangent_matrix(1E-6, 10, 0.00001)



for i in range(100):
   md.solve('noisy', 'max_iter', 50, 'max_res', 1e-8); # , 'lsearch', 'simplest')

   U2 = md.variable('u2')
   VM2 = md.compute_isotropic_linearized_Von_Mises_or_Tresca('u2', 'lambda', 'mu', mfvm2)

   # gf_plot(mfvm2,VM2,'mesh','off', 'deformation',U2,'deformation_mf',mfu2,'deformation_scale', 1, 'refine', 8); colorbar;
   mfvm2.export_to_vtk('large_sliding_contact_2_%i.vtk' % i, 'ascii', mfvm2,  VM2,
                       'Von Mises Stresses 2', mfu2, U2, 'Displacements 2')

   if two_bodies:
      # hold on
      U1 = md.variable('u1')
      VM1 = md.compute_isotropic_linearized_Von_Mises_or_Tresca('u1', 'lambda', 'mu', mfvm1)
      # gf_plot(mfvm1,VM1,'mesh','off', 'deformation',U1,'deformation_mf',mfu1,'deformation_scale', 1, 'refine', 8); colorbar;
      # hold off
      mfvm1.export_to_vtk('large_sliding_contact_1_%i.vtk' % i, 'ascii', mfvm1,  VM1,
                          'Von Mises Stresses 1', mfu1, U1, 'Displacements 1')

   # axis([-2, 2, -0.2, 3])
   # pause(1)

   vf = vf + 0.01;
   md.set_variable('data1', [0,-vf]);
   md.set_variable('data2', [0,-vf]);

