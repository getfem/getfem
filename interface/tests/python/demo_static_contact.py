#!/usr/bin/env python
# -*- coding: UTF8 -*-
# Python GetFEM++ interface
#
# Copyright (C) 2011 Yves Renard.
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

""" Static equilibrium of an elastic solid in contact with a rigid foundation

  This program is used to check that matlab-getfem is working. This is also
  a good example of use of GetFEM++.
"""

import getfem as gf
import numpy as np

# Import the mesh
m = gf.Mesh('load', '../../../tests/meshes/disc_P2_h2.mesh')
# m = gf.Mesh('load', '../../../tests/meshes/disc_P2_h1.mesh')
# m = gf.Mesh('load', '../../../tests/meshes/disc_P2_h0.5.mesh')
# m = gf.Mesh('load', '../../../tests/meshes/disc_P2_h0.25.mesh')
# m = gf.Mesh('load', '../../../tests/meshes/disc_P2_h0.15.mesh')
d = m.dim() # Mesh dimension

# Parameters of the model
Lambda = 1  # Lame coefficient
Mu = 1      # Lame coefficient
friction_coeff = 0.2 # coefficient of friction
r = 0.0001  # Augmentation parameter
version = 9   # 1 : frictionless contact and the basic contact brick
              # 2 : contact with 'static' Coulomb friction and basic contact
              #     brick
              # 3 : frictionless contact and the contact with a
              #     rigid obstacle brick
              # 4 : contact with 'static' Coulomb friction and the contact
              #     with a rigid obstacle brick
              # 5 : frictionless contact and the continuous brick
              #     Newton and Alart-Curnier augmented lagrangian,
              #     unsymmetric version
              # 6 : frictionless contact and the continuous brick
              #     Newton and Alart-Curnier augmented lagrangian, symmetric
              #     version.
              # 7 : frictionless contact and the continuous brick
              #     Newton and Alart-Curnier augmented lagrangian,
              #     unsymmetric version with an additional augmentation.
              # 8 : frictionless contact and the continuous brick
              #     New unsymmetric version.
              # 9 : frictionless contact and the continuous brick : Uzawa
              #     (not very adapted because it is a semi-coercive case)
              # 10 : contact with 'static' Coulomb friction and the continuous
              #     brick. Newton and Alart-Curnier augmented lagrangian,
              #     unsymmetric version.
              # 11 : contact with 'static' Coulomb friction and the continuous
              #     brick. Newton and Alart-Curnier augmented lagrangian,
              #     nearly symmetric version.
              # 12 : contact with 'static' Coulomb friction and the continuous
              #     brick. Newton and Alart-Curnier augmented lagrangian,
              #     unsymmetric version with an additional augmentation.
              # 13 : contact with 'static' Coulomb friction and the continuous
              #     brick. New unsymmetric version.

penalty_parameter = 1E-8    # For rigid motions.
uzawa_r = penalty_parameter # Descent coefficient for Uzawa method.
niter = 50  # Maximum number of iterations for Newton's algorithm.
# Signed distance representing the obstacle
if d == 2:
   obstacle = 'y'
else:
   obstacle = 'z'

# Selection of the contact boundary
border = m.outer_faces()
normals = m.normal_of_faces(border)
contact_boundary = border[:,np.nonzero(normals[d-1] < 0)[0]]
GAMMAC = 1
m.set_region(GAMMAC, contact_boundary)

# Finite element methods
mfu = gf.MeshFem(m, d)
mfu.set_classical_fem(1)

mfd = gf.MeshFem(m, 1)
mfd.set_classical_fem(1)

mfvm = gf.MeshFem(m, 1)
mfvm.set_classical_discontinuous_fem(1)

mflambda = gf.MeshFem(m, 1) # used only by version 5, 6, 7, 8 and 9
mflambda.set_classical_fem(1)

# Integration method
mim = gf.MeshIm(m, 4)
mim_friction = gf.MeshIm(m, gf.Integ('IM_STRUCTURED_COMPOSITE(IM_TRIANGLE(4),4)'))

# Volumic density of force
nbdofd = mfd.nbdof()
nbdofu = mfu.nbdof()
F = np.zeros(nbdofd*d)
F[d-1:nbdofd*d:d] = -0.02;

# Elasticity model
md = gf.Model('real')
md.add_fem_variable('u', mfu)
md.add_initialized_data('mu', [Mu])
md.add_initialized_data('lambda', [Lambda])
md.add_isotropic_linearized_elasticity_brick(mim, 'u', 'lambda', 'mu')
md.add_initialized_fem_data('volumicload', mfd, F)
md.add_source_term_brick(mim, 'u', 'volumicload')

# Small penalty term to avoid rigid motion (should be replaced by an
# explicit treatment of the rigid motion with a constraint matrix)
md.add_initialized_data('penalty_param', [penalty_parameter])
md.add_mass_brick(mim, 'u', 'penalty_param')


# The contact condition

cdof = mfu.dof_on_region(GAMMAC)
nbc = cdof.shape[0] / d

if version == 1 or version == 2: # defining the matrices BN and BT by hand
   contact_dof = cdof[d-1:nbc*d:d]
   contact_nodes = mfu.basic_dof_nodes(contact_dof)
   BN = gf.Spmat('empty', nbc, nbdofu)
   for i in range(nbc):
      BN[i, contact_dof[i]] = -1.
      gap[i] = contact_nodes[d-1, i]

   if version == 2:
      BT = gf.Spmat('empty', nbc*(d-1), nbdofu)
      for i in range(nbc):
         BT[i*(d-1)-1, contact_dof[i]-d+1] = 1.0
         if d > 2:
            BT[i*(d-1), contact_dof[i]-d+2] = 1.0

   md.add_variable('lambda_n', nbc)
   md.add_initialized_data('r', [r])
   if version == 2:
      md.add_variable('lambda_t', nbc * (d-1))
      md.add_initialized_data('friction_coeff', [friction_coeff])

   md.add_initialized_data('gap', gap)
   md.add_initialized_data('alpha', np.ones(nbc))
   if version == 1:
      md.add_basic_contact_brick('u', 'lambda_n', 'r', BN, 'gap', 'alpha', 0)
   else:
      md.add_basic_contact_brick('u', 'lambda_n', 'lambda_t', 'r', BN, BT, 'friction_coeff', 'gap', 'alpha', 0);

elif version == 3 or version == 4: # BN and BT defined by the contact brick

   md.add_variable('lambda_n', nbc)
   md.add_initialized_data('r', [r])
   if version == 3:
      md.add_contact_with_rigid_obstacle_brick(mim, 'u', 'lambda_n', 'r', GAMMAC, obstacle, 0);
   else:
      md.add_variable('lambda_t', nbc*(d-1))
      md.add_initialized_data('friction_coeff', [friction_coeff])
      md.add_contact_with_rigid_obstacle_brick(mim, 'u', 'lambda_n', 'lambda_t', 'r',
                                              'friction_coeff', GAMMAC, obstacle, 0)

elif version >= 5 and version <= 8: # The continuous version, Newton

   ldof = mflambda.dof_on_region(GAMMAC)
   mflambda_partial = gf.MeshFem('partial', mflambda, ldof)
   md.add_fem_variable('lambda_n', mflambda_partial)
   md.add_initialized_data('r', [r])
   OBS = mfd.eval(obstacle)
   md.add_initialized_fem_data('obstacle', mfd, OBS)
   md.add_continuous_contact_with_rigid_obstacle_brick(mim_friction, 'u', 'lambda_n',
                                                      'obstacle', 'r', GAMMAC, version-4);

elif version == 9: # The continuous version, Uzawa

   ldof = mflambda.dof_on_region(GAMMAC)
   mflambda_partial = gf.MeshFem('partial', mflambda, ldof)
   nbc = mflambda_partial.nbdof()
   lambda_n = np.zeros(nbc)
   lzeros = lambda_n
   # Not correct : the normal to the obstacle is assumed to be vertical.
   W = -gf.asm_boundary(GAMMAC, 'a=data(#2);V(#1)+=comp(vBase(#1).Base(#2))(:,2,i).a(i)',
                        mim, mfu, mflambda_partial, lambda_n)
   indb = md.add_explicit_rhs('u', W)
   OBS = np.array([mfd.eval(obstacle)])
   M = gf.asm_mass_matrix(mim, mflambda_partial, mflambda_partial, GAMMAC)
   M.display()
   for ii in range(100000):
      print 'iteration %d' % (ii+1)
      md.solve('max_res', 1E-12, 'very noisy')
      U = md.get('variable', 'u')
      lambda_n_old = lambda_n
      print lambda_n
      lambda_n = np.array([np.linalg.solve(M, gf.asm_contact_Uzawa_projection(GAMMAC, mim_friction, mfu, U, mflambda_partial, lambda_n, mfd, OBS, uzawa_r))])
      W = -gf.asm_boundary(GAMMAC, 'a=data(#2);V(#1)+=comp(vBase(#1).Base(#2))(:,2,i).a(i)', mim, mfu, mflambda_partial, lambda_n)
      md.set_private_rhs(indb, W)
      difff = max(abs(lambda_n-lambda_n_old));
      print 'diff : %g' % difff
      if difff < penalty_parameter:
         break

elif version >= 10 and version <= 13: # Continuous version with friction, Newton

   mflambda.set_qdim(2);
   ldof = mflambda.dof_on_region(GAMMAC)
   mflambda_partial = gf.MeshFem('partial', mflambda, ldof)
   md.add_fem_variable('lambda_n', mflambda_partial)
   md.add_initialized_data('r', [r])
   md.add_initialized_data('friction_coeff', [friction_coeff])
   OBS = mfd.eval(obstacle)
   md.add_initialized_fem_data('obstacle', mfd, OBS)
   md.add_continuous_contact_with_friction_with_rigid_obstacle_brick(mim_friction, 'u', 'lambda_n', 'obstacle', 'r', 'friction_coeff', GAMMAC, version-9);

else:
   print 'Unexistent version'

# Solve the problem
if version != 9:
   md.solve('max_res', 1E-9, 'very noisy', 'max_iter', niter, 'lsearch', 'default') #, 'with pseudo potential')

U = md.get('variable', 'u')
# lambda_n = md.get('variable', 'lambda_n')
VM = md.compute_isotropic_linearized_Von_Mises_or_Tresca('u', 'lambda', 'mu', mfvm)

mfd.export_to_vtk('static_contact.vtk', 'ascii', mfvm,  VM, 'Von Mises Stress', mfu, U, 'Displacement')

