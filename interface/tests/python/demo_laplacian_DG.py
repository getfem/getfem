#!/usr/bin/env python
# Python GetFEM interface
#
# Copyright (C) 2015-2020 Yves Renard
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

  This program is used to check that python-getfem is working. This is
  also a good example of use of GetFEM.

  Poisson problem solved with a Discontinuous Galerkin method (or
  interior penalty method). See for instance
  "Unified analysis of discontinuous Galerkin methods for elliptic
  problems", D.N. Arnold, F. Brezzi, B. Cockburn, L.D. Marini, SIAM J.
  Numer. Anal. vol. 39:5, pp 1749-1779, 2002.

  $Id: demo_laplacian_DG.py 4429 2013-10-01 13:15:15Z renard $
"""
import numpy as np

# Import basic modules
import getfem as gf

## Parameters
NX = 20                            # Mesh parameter.
Dirichlet_with_multipliers = True  # Dirichlet condition with multipliers
                                   # or penalization
dirichlet_coefficient = 1e10       # Penalization coefficient
interior_penalty_factor = 1e2*NX   # Parameter of the interior penalty term
verify_neighbor_computation = True;


# Create a simple cartesian mesh
m = gf.Mesh('regular_simplices', np.arange(0,1+1./NX,1./NX), np.arange(0,1+1./NX,1./NX))

# Create a MeshFem for u and rhs fields of dimension 1 (i.e. a scalar field)
mfu   = gf.MeshFem(m, 1)
mfrhs = gf.MeshFem(m, 1)
# Assign the discontinuous P2 fem to all convexes of the both MeshFem
mfu.set_fem(gf.Fem('FEM_PK_DISCONTINUOUS(2,2)'))
mfrhs.set_fem(gf.Fem('FEM_PK(2,2)'))

# Integration method used
mim = gf.MeshIm(m, gf.Integ('IM_TRIANGLE(4)'))

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

# Inner edges for the interior penalty terms
in_faces = m.inner_faces()
INNER_FACES=4
m.set_region(INNER_FACES, in_faces)

if (verify_neighbor_computation):
  TEST_FACES=5
  adjf = m.adjacent_face(42, 0);
  if (len(adjf) != 2):
    print('No adjacent edge found, change the element number')
    exit(1)
  m.set_region(TEST_FACES, np.array([[42,adjf[0][0]], [0,adjf[1][0]]]));
  

# Interpolate the exact solution (Assuming mfu is a Lagrange fem)
Ue = mfu.eval('y*(y-1)*x*(x-1)+x*x*x*x*x')

# Interpolate the source term
F1 = mfrhs.eval('-(2*(x*x+y*y)-2*x-2*y+20*x*x*x)')
F2 = mfrhs.eval('[y*(y-1)*(2*x-1) + 5*x*x*x*x, x*(x-1)*(2*y-1)]')

# Model
md = gf.Model('real')

# Main unknown
md.add_fem_variable('u', mfu)

# Laplacian term on u
md.add_Laplacian_brick(mim, 'u')

# Volumic source term
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

# Interior penalty terms
md.add_initialized_data('alpha', [interior_penalty_factor])
jump = "((u-Interpolate(u,neighbor_element))*Normal)"
test_jump = "((Test_u-Interpolate(Test_u,neighbor_element))*Normal)"
grad_mean = "((Grad_u+Interpolate(Grad_u,neighbor_element))*0.5)"
grad_test_mean = "((Grad_Test_u+Interpolate(Grad_Test_u,neighbor_element))*0.5)"
md.add_linear_term(mim, "-(({F}).({G}))-(({H}).({I}))+alpha*(({J}).({K}))".format(F=grad_mean, G=test_jump, H=jump, I=grad_test_mean, J=jump, K=test_jump), INNER_FACES);

gf.memstats()
# md.listvar()
# md.listbricks()

# Assembly of the linear system and solve.
md.solve()

# Main unknown
U = md.variable('u')
L2error = gf.compute(mfu, U-Ue, 'L2 norm', mim)
H1error = gf.compute(mfu, U-Ue, 'H1 norm', mim)
print('Error in L2 norm : ', L2error)
print('Error in H1 norm : ', H1error)

# Export data
mfu.export_to_pos('laplacian.pos', Ue,'Exact solution',
                                    U,'Computed solution')
print('You can view the solution with (for example):')
print('gmsh laplacian.pos')

if (verify_neighbor_computation):
  A=gf.asm('generic', mim, 1, 'u*Test_u*(Normal.Normal)', TEST_FACES, md)
  B=gf.asm('generic', mim, 1, '-Interpolate(u,neighbor_element)*Interpolate(Test_u,neighbor_element)*(Interpolate(Normal,neighbor_element).Normal)', TEST_FACES, md)
  err_v = np.linalg.norm(A-B)
  A=gf.asm('generic', mim, 1, '(Grad_u.Normal)*(Grad_Test_u.Normal)', TEST_FACES, md)
  B=gf.asm('generic', mim, 1, '(Interpolate(Grad_u,neighbor_element).Normal)*(Interpolate(Grad_Test_u,neighbor_element).Normal)', TEST_FACES, md)
  err_v = err_v + np.linalg.norm(A-B)
  if (err_v > 1E-13):
    print('Test on neighbor element computation: error to big: ', err_v)
    exit(1)
  
if (H1error > 1e-3):
    print('Error too large !')
    exit(1)
