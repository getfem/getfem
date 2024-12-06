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
"""  2D-3D bilaplacian problem test.

  This program is used to check that python-getfem is working. This is
  also a good example of use of GetFEM.

  $Id$
"""

import numpy as np
import getfem as gf


## Parameters
N = 3                                   # Dimension of the problem
NN = 10
NX = NN; NY = NN; NZ = NN;               # Mesh parameters.

if (N == 2):
  m = gf.Mesh('regular_simplices', np.arange(0, 1.+1./NX, 1./NX),
              np.arange(0, 1.+1./NY, 1./NY))
else:
  m = gf.Mesh('regular_simplices', np.arange(0, 1.+1./NX, 1./NX),
              np.arange(0, 1.+1./NY, 1./NY), np.arange(0,1.+1./NZ,1./NZ))

useKL=1 # use the Kirchhoff-Love plate model, or just a pure
        # bilaplacian problem

D=1.    # Flexion modulus
NU=0.3  # Poisson ratio (0 <= NU <= 1) for KL version

if (N == 2):
  mim = gf.MeshIm(m, gf.Integ('IM_TRIANGLE(10)'))
  mfu = gf.MeshFem(m, 1)
  mfd = gf.MeshFem(m, 1)
  # mfu.set_fem(gf.Fem('FEM_ARGYRIS'))
  # mfu.set_fem(gf.Fem('FEM_HERMITE(2)'))
  mfu.set_fem(gf.Fem('FEM_MORLEY(2)'))
  mfd.set_fem(gf.Fem('FEM_PK(2,5)'))
else:
  mim = gf.MeshIm(m, 6)
  mfu = gf.MeshFem(m, 1)
  mfd = gf.MeshFem(m, 1)
  mfu.set_fem(gf.Fem('FEM_MORLEY(3)'))
  # mfu.set_fem(gf.Fem('FEM_HERMITE(3)'))
  mfd.set_fem(gf.Fem('FEM_PK(3,3)'))

ftot    = m.outer_faces()

# FORCE_BOUNDARY_NUM=1;
# MOMENTUM_BOUNDARY_NUM=2;
SIMPLE_SUPPORT_BOUNDARY_NUM=3;
CLAMPED_BOUNDARY_NUM=4;

m.set_region(SIMPLE_SUPPORT_BOUNDARY_NUM, ftot)
m.set_region(CLAMPED_BOUNDARY_NUM, ftot)

if (N == 2):
  sol_u=mfd.eval('((x**2)*((1-x)**2) + (y**2)*((1-y)**2))', globals(),locals());
  F = 48
else:
  sol_u=mfd.eval('((x**2)*((1-x)**2) + (y**2)*((1-y)**2) + (z**2)*((1-z)**2))', globals(),locals());
  F = 72

md=gf.Model('real');
md.add_fem_variable('u', mfu);

if useKL:
  md.add_initialized_data('D', D);
  md.add_initialized_data('nu', NU);
  md.add_Kirchhoff_Love_plate_brick(mim, 'u', 'D', 'nu');
  # M = np.zeros((N, N, mfd.nbdof()));
else :
  md.add_initialized_data('D', D);
  md.add_bilaplacian_brick(mim, 'u', 'D');
  # M = np.zeros((1, mfd.nbdof()));


md.add_initialized_fem_data('VolumicData', mfd, mfd.eval('%g' % F));
md.add_source_term_brick(mim, 'u', 'VolumicData');

#md.add_initialized_fem_data('M', mfd, M);
#md.add_normal_derivative_source_term_brick(mim, 'u', 'M', MOMENTUM_BOUNDARY_NUM);

#if (useKL): 
#  H = np.zeros((N, N, mfd.nbdof()));
#  FF = np.zeros((N, mfd.nbdof()));
#  md.add_initialized_fem_data('H', mfd, H);
#  md.add_initialized_fem_data('F', mfd, FF);
#  md.add_Kirchhoff_Love_Neumann_term_brick(mim, 'u', 'H', 'F', FORCE_BOUNDARY_NUM);
#else:
#  FF = np.zeros((1, N, mfd.nbdof()));
#  md.add_initialized_fem_data('F', mfd, FF);
#  md.add_normal_source_term_brick(mim, 'u', 'F', FORCE_BOUNDARY_NUM);

md.add_initialized_fem_data('SOL_U', mfd, sol_u)
md.add_normal_derivative_Dirichlet_condition_with_penalization(mim, 'u', 1e10, CLAMPED_BOUNDARY_NUM);
md.add_Dirichlet_condition_with_penalization(mim, 'u', 1e10, SIMPLE_SUPPORT_BOUNDARY_NUM, 'SOL_U');

md.solve('noisy');
U = md.variable('u');

Ud = gf.compute(mfu,U,'interpolate on',mfd)
err= Ud - sol_u;
print(gf.compute_L2_norm(mfd, err, mim));

mfd.export_to_pos('laplacian.pos', Ud, 'computed solution', err, 'error', sol_u, 'exact solution')















