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
"""  2D Poisson problem test.

  This program is used to check that python-getfem is working in parallel.
  This is also a good example of use of GetFEM.

  Run this script by invoking
  mpiexec -n 4 python demo_parallel_laplacian.py

  $Id: demo_parallel_laplacian.py 3809 2011-09-26 20:38:56Z logari81 $
"""
import time

import numpy as np

import getfem as gf
# import basic modules
import mpi4py.MPI as mpi

rank = mpi.COMM_WORLD.rank
if (rank == 0):
  print('Running Parallel Getfem with python interface')

print('Hello from thread ', rank)

## Parameters
NX = 100                           # Mesh parameter.
Dirichlet_with_multipliers = True  # Dirichlet condition with multipliers
                                   # or penalization
dirichlet_coefficient = 1e10       # Penalization coefficient

t = time.process_time()

# creation of a simple cartesian mesh
m = gf.Mesh('regular_simplices', np.arange(0,1+1./NX,1./NX), np.arange(0,1+1./NX,1./NX))
if (rank == 0):
  print('Time for building mesh', time.process_time()-t)
t = time.process_time()

# create a MeshFem for u and rhs fields of dimension 1 (i.e. a scalar field)
mfu   = gf.MeshFem(m, 1)
mfrhs = gf.MeshFem(m, 1)
# assign the P2 fem to all convexes of the both MeshFem
mfu.set_fem(gf.Fem('FEM_PK(2,2)'))
mfrhs.set_fem(gf.Fem('FEM_PK(2,2)'))

# an exact integration will be used
mim = gf.MeshIm(m, gf.Integ('IM_TRIANGLE(4)'))

# boundary selection
flst  = m.outer_faces()
fnor  = m.normal_of_faces(flst)
tleft = abs(fnor[1,:]+1) < 1e-14
ttop  = abs(fnor[0,:]-1) < 1e-14
fleft = np.compress(tleft, flst, axis=1)
ftop  = np.compress(ttop, flst, axis=1)
fneum = np.compress(np.logical_not(ttop + tleft), flst, axis=1)

# mark it as boundary
DIRICHLET_BOUNDARY_NUM1 = 1
DIRICHLET_BOUNDARY_NUM2 = 2
NEUMANN_BOUNDARY_NUM = 3
m.set_region(DIRICHLET_BOUNDARY_NUM1, fleft)
m.set_region(DIRICHLET_BOUNDARY_NUM2, ftop)
m.set_region(NEUMANN_BOUNDARY_NUM, fneum)

if (rank == 0):
  print('Time for building fem and im', time.process_time()-t)
t = time.process_time()

nb_dof = mfu.nbdof()
if (rank == 0):
  print('Nb dof for the main unknown: ', nb_dof)

if (rank == 0):
  print('Time for dof numbering', time.process_time()-t)
t = time.process_time()
  

# interpolate the exact solution (Assuming mfu is a Lagrange fem)
Ue = mfu.eval('y*(y-1)*x*(x-1)+x*x*x*x*x')

# interpolate the source term
F1 = mfrhs.eval('-(2*(x*x+y*y)-2*x-2*y+20*x*x*x)')
F2 = mfrhs.eval('[y*(y-1)*(2*x-1) + 5*x*x*x*x, x*(x-1)*(2*y-1)]')

if (rank == 0):
  print('Time for python interpolation', time.process_time()-t)
t = time.process_time()

# model
md = gf.Model('real')

# main unknown
md.add_fem_variable('u', mfu)

# laplacian term on u
md.add_Laplacian_brick(mim, 'u')

# volumic source term
md.add_initialized_fem_data('VolumicData', mfrhs, F1)
md.add_source_term_brick(mim, 'u', 'VolumicData')

# Neumann condition.
md.add_initialized_fem_data('NeumannData', mfrhs, F2)
md.add_normal_source_term_brick(mim, 'u', 'NeumannData',
                                NEUMANN_BOUNDARY_NUM)

# Dirichlet condition on the left.
md.add_initialized_fem_data("DirichletData", mfu, Ue)

if (Dirichlet_with_multipliers):
  md.add_Dirichlet_condition_with_multipliers(mim, 'u', mfu,
                                              DIRICHLET_BOUNDARY_NUM1,
                                              'DirichletData')
else:
  md.add_Dirichlet_condition_with_penalization(mim, 'u', dirichlet_coefficient,
                                               DIRICHLET_BOUNDARY_NUM1,
                                               'DirichletData')

# Dirichlet condition on the top.
# Two Dirichlet brick in order to test the multiplier
# selection in the intersection.
if (Dirichlet_with_multipliers):
  md.add_Dirichlet_condition_with_multipliers(mim, 'u', mfu,
                                              DIRICHLET_BOUNDARY_NUM2,
                                              'DirichletData')
else:
  md.add_Dirichlet_condition_with_penalization(mim, 'u', dirichlet_coefficient,
                                               DIRICHLET_BOUNDARY_NUM2,
                                               'DirichletData')

if (rank == 0):
  print('Time for model building', time.process_time()-t)
t = time.process_time()

md.nbdof
nb_dof = md.nbdof()
if (rank == 0):
  print('Nb dof for the model: ', nb_dof)

if (rank == 0):
  print('Time for model actualize sizes', time.process_time()-t)
t = time.process_time()

# assembly of the linear system and solve.
md.solve()

if (rank == 0):
  print('Time for model solve', time.process_time()-t)
t = time.process_time()

# main unknown
U = md.variable('u')
L2error = gf.compute(mfu, U-Ue, 'L2 norm', mim)
H1error = gf.compute(mfu, U-Ue, 'H1 norm', mim)

if (rank == 0):
  print('Error in L2 norm : ', L2error)
  print('Error in H1 norm : ', H1error)

if (rank == 0):
  print('Time for error computation', time.process_time()-t)
t = time.process_time()


# export data
# if (rank == 0):
#   mfu.export_to_pos('laplacian.pos', Ue,'Exact solution',
#                     U,'Computed solution')
#   print('You can view the solution with (for example):')
#   print('gmsh laplacian.pos')
