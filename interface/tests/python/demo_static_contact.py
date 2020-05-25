#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Python GetFEM interface
#
# Copyright (C) 2011-2020 Yves Renard.
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

""" Static equilibrium of an elastic solid in contact with a rigid foundation

  This program is used to check that python-getfem is working. This is also
  a good example of use of GetFEM.
"""

import numpy as np

import getfem as gf

# Import the mesh : disc
# m = gf.Mesh('load', '../../../tests/meshes/disc_P2_h4.mesh')
#m = gf.Mesh('load', '../../../tests/meshes/disc_P2_h2.mesh')
# m = gf.Mesh('load', '../../../tests/meshes/disc_P2_h1.mesh')
# m = gf.Mesh('load', '../../../tests/meshes/disc_P2_h0_5.mesh')
# m = gf.Mesh('load', '../../../tests/meshes/disc_P2_h0_3.mesh')

# Import the mesh : sphere
# m = gf.Mesh('load', '../../../tests/meshes/sphere_with_quadratic_tetra_8_elts.mesh')
# m = gf.Mesh('load', '../../../tests/meshes/sphere_with_quadratic_tetra_80_elts.mesh')
m = gf.Mesh('load', '../../../tests/meshes/sphere_with_quadratic_tetra_400_elts.mesh')
# m = gf.Mesh('load', '../../../tests/meshes/sphere_with_quadratic_tetra_2000_elts.mesh')
# m = gf.Mesh('load', '../../../tests/meshes/sphere_with_quadratic_tetra_16000_elts.mesh')

d = m.dim() # Mesh dimension

# Parameters of the model
clambda = 1.          # Lame coefficient
cmu = 1.              # Lame coefficient
friction_coeff = 0.4  # coefficient of friction
vertical_force = 0.05 # Volumic load in the vertical direction
r = 10.               # Augmentation parameter
condition_type = 0 # 0 = Explicitely kill horizontal rigid displacements
                   # 1 = Kill rigid displacements using a global penalization
                   # 2 = Add a Dirichlet condition on the top of the structure
penalty_parameter = 1E-6    # Penalization coefficient for the global penalization

if d == 2:
   cpoints = [0, 0]   # constrained points for 2d
   cunitv  = [1, 0]   # corresponding constrained directions for 2d
else:
   cpoints = [0, 0, 0,   0, 0, 0,   5, 0, 5]  # constrained points for 3d
   cunitv  = [1, 0, 0,   0, 1, 0,   0, 1, 0]  # corresponding constrained directions for 3d

niter = 100  # Maximum number of iterations for Newton's algorithm.
version = 13  # 1 : frictionless contact and the basic contact brick
              # 2 : contact with 'static' Coulomb friction and basic contact brick
              # 3 : frictionless contact and the contact with a
              #     rigid obstacle brick
              # 4 : contact with 'static' Coulomb friction and the contact with a
              #     rigid obstacle brick
              # 5 : frictionless contact and the integral brick
              #     Newton and Alart-Curnier augmented lagrangian,
              #     unsymmetric version
              # 6 : frictionless contact and the integral brick
              #     Newton and Alart-Curnier augmented lagrangian, symmetric
              #     version.
              # 7 : frictionless contact and the integral brick
              #     Newton and Alart-Curnier augmented lagrangian,
              #     unsymmetric version with an additional augmentation.
              # 8 : frictionless contact and the integral brick
              #     New unsymmetric method.
              # 9 : frictionless contact and the integral brick : Uzawa
              #     on the Lagrangian augmented by the penalization term.
              # 10 : contact with 'static' Coulomb friction and the integral brick
              #     Newton and Alart-Curnier augmented lagrangian,
              #     unsymmetric version.
              # 11 : contact with 'static' Coulomb friction and the integral brick
              #     Newton and Alart-Curnier augmented lagrangian,
              #     nearly symmetric version.
              # 12 : contact with 'static' Coulomb friction and the integral brick
              #     Newton and Alart-Curnier augmented lagrangian,
              #     unsymmetric version with an additional augmentation.
              # 13 : contact with 'static' Coulomb friction and the integral brick
              #     New unsymmetric method.
              # 14 : contact with 'static' Coulomb friction and the integral brick : Uzawa
              #     on the Lagrangian augmented by the penalization term.
              # 15 : penalized contact with 'static' Coulomb friction (r is the penalization
              #     coefficient).

# Signed distance representing the obstacle
if d == 2:
   obstacle = 'y'
else:
   obstacle = 'z'

# Selection of the contact and Dirichlet boundaries
GAMMAC = 1
GAMMAD = 2

border = m.outer_faces()
normals = m.normal_of_faces(border)
contact_boundary = border[:,np.nonzero(normals[d-1] < -0.01)[0]]
m.set_region(GAMMAC, contact_boundary)
contact_boundary = border[:,np.nonzero(normals[d-1] > 0.01)[0]]
m.set_region(GAMMAD, contact_boundary)

# Finite element methods
u_degree = 2
lambda_degree = 2

mfu = gf.MeshFem(m, d)
mfu.set_classical_fem(u_degree)

mfd = gf.MeshFem(m, 1)
mfd.set_classical_fem(u_degree)

mflambda = gf.MeshFem(m, 1) # used only by version 5 to 13
mflambda.set_classical_fem(lambda_degree)

mfvm = gf.MeshFem(m, 1)
mfvm.set_classical_discontinuous_fem(u_degree-1)

# Integration method
mim = gf.MeshIm(m, 4)
if d == 2:
   mim_friction = gf.MeshIm(m,
     gf.Integ('IM_STRUCTURED_COMPOSITE(IM_TRIANGLE(4),4)'))
else:
   mim_friction = gf.MeshIm(m,
     gf.Integ('IM_STRUCTURED_COMPOSITE(IM_TETRAHEDRON(5),4)'))

# Volumic density of force
nbdofd = mfd.nbdof()
nbdofu = mfu.nbdof()
F = np.zeros(nbdofd*d)
F[d-1:nbdofd*d:d] = -vertical_force;

# Elasticity model
md = gf.Model('real')
md.add_fem_variable('u', mfu)
md.add_initialized_data('cmu', [cmu])
md.add_initialized_data('clambda', [clambda])
md.add_isotropic_linearized_elasticity_brick(mim, 'u', 'clambda', 'cmu')
md.add_initialized_fem_data('volumicload', mfd, F)
md.add_source_term_brick(mim, 'u', 'volumicload')

if condition_type == 2:
   Ddata = np.zeros(d)
   Ddata[d-1] = -5
   md.add_initialized_data('Ddata', Ddata)
   md.add_Dirichlet_condition_with_multipliers(mim, 'u', u_degree, GAMMAD, 'Ddata')
elif condition_type == 0:
   md.add_initialized_data('cpoints', cpoints)
   md.add_initialized_data('cunitv', cunitv)
   md.add_pointwise_constraints_with_multipliers('u', 'cpoints', 'cunitv')
elif condition_type == 1:
   # Small penalty term to avoid rigid motion (should be replaced by an
   # explicit treatment of the rigid motion with a constraint matrix)
   md.add_initialized_data('penalty_param', [penalty_parameter])
   md.add_mass_brick(mim, 'u', 'penalty_param')

# The contact condition

cdof = mfu.dof_on_region(GAMMAC)
nbc = cdof.shape[0] / d

solved = False
if version == 1 or version == 2: # defining the matrices BN and BT by hand
   contact_dof = cdof[d-1:nbc*d:d]
   contact_nodes = mfu.basic_dof_nodes(contact_dof)
   BN = gf.Spmat('empty', nbc, nbdofu)
   ngap = np.zeros(nbc)
   for i in range(nbc):
      BN[i, contact_dof[i]] = -1.
      ngap[i] = contact_nodes[d-1, i]

   if version == 2:
      BT = gf.Spmat('empty', nbc*(d-1), nbdofu)
      for i in range(nbc):
         for j in range(d-1):
            BT[j+i*(d-1), contact_dof[i]-d+j+1] = 1.0

   md.add_variable('lambda_n', nbc)
   md.add_initialized_data('r', [r])
   if version == 2:
      md.add_variable('lambda_t', nbc * (d-1))
      md.add_initialized_data('friction_coeff', [friction_coeff])

   md.add_initialized_data('ngap', ngap)
   md.add_initialized_data('alpha', np.ones(nbc))
   if version == 1:
      md.add_basic_contact_brick('u', 'lambda_n', 'r', BN, 'ngap', 'alpha', 1)
   else:
      md.add_basic_contact_brick('u', 'lambda_n', 'lambda_t', 'r', BN, BT, 'friction_coeff', 'ngap', 'alpha', 1);

elif version == 3 or version == 4: # BN and BT defined by the contact brick

   md.add_variable('lambda_n', nbc)
   md.add_initialized_data('r', [r])
   if version == 3:
      md.add_nodal_contact_with_rigid_obstacle_brick(mim, 'u', 'lambda_n', 'r', GAMMAC, obstacle, 1);
   else:
      md.add_variable('lambda_t', nbc*(d-1))
      md.add_initialized_data('friction_coeff', [friction_coeff])
      md.add_nodal_contact_with_rigid_obstacle_brick(mim, 'u', 'lambda_n', 'lambda_t', 'r',
                                                     'friction_coeff', GAMMAC, obstacle, 1)

elif version >= 5 and version <= 8: # The integral version, Newton

   ldof = mflambda.dof_on_region(GAMMAC)
   mflambda_partial = gf.MeshFem('partial', mflambda, ldof)
   md.add_fem_variable('lambda_n', mflambda_partial)
   md.add_initialized_data('r', [r])
   OBS = mfd.eval(obstacle)
   md.add_initialized_fem_data('obstacle', mfd, OBS)
   md.add_integral_contact_with_rigid_obstacle_brick(mim_friction, 'u', 'lambda_n',
                                                      'obstacle', 'r', GAMMAC, version-4);

elif version == 9: # The integral version, Uzawa on the augmented Lagrangian

   ldof = mflambda.dof_on_region(GAMMAC)
   mflambda_partial = gf.MeshFem('partial', mflambda, ldof)
   nbc = mflambda_partial.nbdof()
   M = gf.asm_mass_matrix(mim, mflambda_partial, mflambda_partial, GAMMAC)
   lambda_n = np.zeros(nbc)
   md.add_initialized_fem_data('lambda_n', mflambda_partial, lambda_n)
   md.add_initialized_data('r', [r])
   OBS = mfd.eval(obstacle) # np.array([mfd.eval(obstacle)])
   md.add_initialized_fem_data('obstacle', mfd, OBS)
   md.add_penalized_contact_with_rigid_obstacle_brick \
     (mim_friction, 'u', 'obstacle', 'r', GAMMAC, 2, 'lambda_n')

   for ii in range(100):
      print('iteration %d' % (ii+1))
      md.solve('max_res', 1E-9, 'max_iter', niter)
      U = md.get('variable', 'u')
      lambda_n_old = lambda_n
      sol = gf.linsolve_superlu(M, gf.asm_integral_contact_Uzawa_projection(GAMMAC, mim_friction, mfu, U, mflambda_partial, lambda_n, mfd, OBS, r))
      lambda_n = sol[0].transpose()
      md.set_variable('lambda_n', lambda_n)
      difff = max(abs(lambda_n-lambda_n_old))[0]/max(abs(lambda_n))[0]
      print('diff : %g' % difff)
      if difff < penalty_parameter:
         break

   solved = True

elif version >= 10 and version <= 13: # The integral version with friction, Newton

   mflambda.set_qdim(d);
   ldof = mflambda.dof_on_region(GAMMAC)
   mflambda_partial = gf.MeshFem('partial', mflambda, ldof)
   md.add_fem_variable('lambda', mflambda_partial)
   md.add_initialized_data('r', [r])
   md.add_initialized_data('friction_coeff', [friction_coeff])
   OBS = mfd.eval(obstacle)
   md.add_initialized_fem_data('obstacle', mfd, OBS)
   md.add_integral_contact_with_rigid_obstacle_brick \
     (mim_friction, 'u', 'lambda', 'obstacle', 'r', 'friction_coeff', GAMMAC, version-9)

elif version == 14: # The integral version, Uzawa on the augmented Lagrangian with friction
  
   mflambda.set_qdim(d)
   ldof = mflambda.dof_on_region(GAMMAC)
   mflambda_partial = gf.MeshFem('partial', mflambda, ldof)
   nbc = mflambda_partial.nbdof()
   md.add_initialized_data('friction_coeff', [friction_coeff])
   M = gf.asm_mass_matrix(mim, mflambda_partial, mflambda_partial, GAMMAC)
   lambda_nt = np.zeros(nbc)
   md.add_initialized_fem_data('lambda', mflambda_partial, lambda_nt)
   md.add_initialized_data('r', [r])
   OBS = mfd.eval(obstacle)
   md.add_initialized_fem_data('obstacle', mfd, OBS)
   md.add_penalized_contact_with_rigid_obstacle_brick \
     (mim_friction, 'u', 'obstacle', 'r', 'friction_coeff', GAMMAC, 2, 'lambda')

   for ii in range(100):
      print('iteration %d' % (ii+1))
      md.solve('max_res', 1E-9, 'max_iter', niter)
      U = md.get('variable', 'u')
      lambda_nt_old = lambda_nt
      sol = gf.linsolve_superlu(M,
        gf.asm_integral_contact_Uzawa_projection(
        GAMMAC, mim_friction, mfu, U, mflambda_partial, lambda_nt, mfd, OBS, r, friction_coeff))
      lambda_nt = sol[0].transpose()
      md.set_variable('lambda', lambda_nt)
      difff = max(abs(lambda_nt-lambda_nt_old))[0]/max(abs(lambda_nt))[0]
      print('diff : %g' % difff)
      if difff < penalty_parameter:
         break
 
   solved = True

elif version == 15:
 
   md.add_initialized_data('r', [r])
   md.add_initialized_data('friction_coeff', [friction_coeff])
   OBS = mfd.eval(obstacle)
   md.add_initialized_fem_data('obstacle', mfd, OBS);
   md.add_penalized_contact_with_rigid_obstacle_brick \
     (mim_friction, 'u', 'obstacle', 'r', 'friction_coeff', GAMMAC)

else:
   print('Inexistent version')

# Solve the problem
if not solved:
   md.solve('max_res', 1E-9, 'very noisy', 'max_iter', niter, 'lsearch', 'default') #, 'with pseudo potential')

U = md.get('variable', 'u')
# LAMBDA = md.get('variable', 'lambda_n')
VM = md.compute_isotropic_linearized_Von_Mises_or_Tresca('u', 'clambda', 'cmu', mfvm)

mfd.export_to_vtk('static_contact.vtk', 'ascii', mfvm,  VM, 'Von Mises Stress', mfu, U, 'Displacement')
