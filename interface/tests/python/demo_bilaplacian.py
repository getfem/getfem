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

# Import basic modules
import getfem as gf

N = 3;

## Parameters
NX = 3; NY = 3; NZ = 3;                  # Mesh parameter.
if (N == 2):
  m = gf.Mesh('regular_simplices', np.arange(0,0.4+0.4/NX,0.4/NX),
              np.arange(0,1.2+1.2/NY,1.2/NY))
else:
  m = gf.Mesh('regular_simplices', np.arange(0,0.4+0.4/NX,0.4/NX),
              np.arange(0,1.2+1.2/NY,1.2/NY), np.arange(0,1.+1./NZ,1./NZ))

useKL=0; # use the Kirchhoff-Love plate model, or just a pure
         # bilaplacian problem

D=1.0;   # Flexion modulus

if (useKL): NU=0.3 # poisson ratio (0 <= NU <= 1)

if (N == 2):
  mim = gf.MeshIm(m, gf.Integ('IM_TRIANGLE(13)'))
  mfu = gf.MeshFem(m, 1)
  mfd = gf.MeshFem(m, 1)
  mfu.set_fem(gf.Fem('FEM_ARGYRIS'))
  mfd.set_fem(gf.Fem('FEM_PK(2,5)'))
  ftop    = m.outer_faces_with_direction([1, 0], 0.1);
  fbottom = m.outer_faces_with_direction([-1, 0], 0.1);
  fleft   = m.outer_faces_with_direction([0, -1], 0.1);
  fright  = m.outer_faces_with_direction([0,  1], 0.1);

else:
  mim = gf.MeshIm(m, 4)
  mfu = gf.MeshFem(m, 1)
  mfd = gf.MeshFem(m, 1)
  mfu.set_fem(gf.Fem('FEM_MORLEY(3)'))
  mfd.set_fem(gf.Fem('FEM_PK(3,2)'))
  ftop    = m.outer_faces_with_direction([1, 0, 0], 0.1);
  fbottom = m.outer_faces_with_direction([-1, 0, 0], 0.1);
  fleft   = m.outer_faces_with_direction([0, -1, 0], 0.1);
  fright  = m.outer_faces_with_direction([0,  1, 0], 0.1);

FORCE_BOUNDARY_NUM=1;
MOMENTUM_BOUNDARY_NUM=2;
SIMPLE_SUPPORT_BOUNDARY_NUM=3;
CLAMPED_BOUNDARY_NUM=4;

m.set_region(FORCE_BOUNDARY_NUM, fright);
m.set_region(SIMPLE_SUPPORT_BOUNDARY_NUM,
             np.concatenate((fleft, ftop, fbottom), axis=1));
m.set_region(CLAMPED_BOUNDARY_NUM,
             np.concatenate((fleft, ftop, fbottom), axis=1));
m.set_region(MOMENTUM_BOUNDARY_NUM, np.concatenate((ftop, fbottom), axis=1));


FT=2.;
sol_u=mfd.eval('np.sin(%g*(x+y))' % FT, globals(),locals());
sol_f=sol_u*FT*FT*FT*FT*N*N*D;
sol_lapl_u=-FT*FT*sol_u*N;

md=gf.Model('real');
md.add_fem_variable('u', mfu);

if useKL:
  md.add_initialized_data('D', D);
  md.add_initialized_data('nu', NU);
  md.add_Kirchhoff_Love_plate_brick(mim, 'u', 'D', 'nu');
  M = np.zeros((N,N, mfd.nbdof()));
else :
  md.add_initialized_data('D', D);
  md.add_bilaplacian_brick(mim, 'u', 'D');
  M = np.zeros((1, mfd.nbdof()));

md.add_initialized_fem_data('VolumicData', mfd,
                            mfd.eval('1-(x-y)**2'));

md.add_source_term_brick(mim, 'u', 'VolumicData');

md.add_initialized_fem_data('M', mfd, M);
md.add_normal_derivative_source_term_brick(mim, 'u', 'M', MOMENTUM_BOUNDARY_NUM);

if (useKL): 
  H = np.zeros((N, N, mfd.nbdof()));
  F = np.zeros((N, mfd.nbdof()));
  md.add_initialized_fem_data('H', mfd, H);
  md.add_initialized_fem_data('F', mfd, F);
  md.add_Kirchhoff_Love_Neumann_term_brick(mim, 'u', 'H', 'F', FORCE_BOUNDARY_NUM);
else:
  F = np.zeros((1, N, mfd.nbdof()));
  md.add_initialized_fem_data('F', mfd, F);
  md.add_normal_source_term_brick(mim, 'u', 'F', FORCE_BOUNDARY_NUM);


md.add_normal_derivative_Dirichlet_condition_with_penalization(mim, 'u', 1e10, CLAMPED_BOUNDARY_NUM);
 
md.add_Dirichlet_condition_with_penalization(mim, 'u', 1e10, SIMPLE_SUPPORT_BOUNDARY_NUM);

md.solve('noisy');
U = md.variable('u');

err=gf.compute(mfu,U,'interpolate on',mfd) - sol_u;
print(gf.compute_L2_norm(mfd, err, mim));
print(err)

















