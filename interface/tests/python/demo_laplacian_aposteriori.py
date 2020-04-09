#!/usr/bin/env python
# -*- coding: utf-8 -*-
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

  Poisson problem solved a refined with an a posteriori error estimate
  Example of integration on edges of normal derivative jump.

  $Id: demo_laplacian_aposteriori.py 4429 2013-10-01 13:15:15Z renard $
"""
import numpy as np

# Import basic modules
import getfem as gf

## Parameters
h = 4.                             # Mesh parameter.
Dirichlet_with_multipliers = True  # Dirichlet condition with multipliers
                                   # or penalization
dirichlet_coefficient = 1e10       # Penalization coefficient
export_mesh = True                 # Draw the mesh after mesh generation or not


# Create a mesh
mo1 = gf.MesherObject('rectangle', [0., 50.], [100., 100.])
mo2 = gf.MesherObject('rectangle', [50., 0.], [100., 100.])
mo3 = gf.MesherObject('union', mo1, mo2)
mo4 = gf.MesherObject('ball', [25., 75], 8.)
mo5 = gf.MesherObject('ball', [75., 25.], 8.)
mo6 = gf.MesherObject('ball', [75., 75.], 8.)
mo7 = gf.MesherObject('union', mo4, mo5, mo6)
mo  = gf.MesherObject('set minus', mo3, mo7)

gf.util('trace level', 2)   # No trace for mesh generation
mesh = gf.Mesh('generate', mo, h, 3)

#
# Boundary selection
#
fb1 = mesh.outer_faces_with_direction([-1., 0.], 0.01) # Left   (Dirichlet)
fb2 = mesh.outer_faces_with_direction([0., -1.], 0.01) # Bottom (Neumann)
fb3 = mesh.outer_faces_in_box([-1., 10.], [101., 101.]) 
fb4 = mesh.outer_faces_in_box([10., -1.], [101., 101.]) 
LEFT_BOUND=1; BOTTOM_BOUND=2; AUX_BOUND1 = 3; AUX_BOUND2 = 4; 
mesh.set_region(  LEFT_BOUND, fb1)
mesh.set_region(BOTTOM_BOUND, fb2)
mesh.set_region(  AUX_BOUND1, fb3)
mesh.set_region(  AUX_BOUND2, fb4)
mesh.region_subtract(  LEFT_BOUND, AUX_BOUND2)
mesh.region_subtract(BOTTOM_BOUND, AUX_BOUND1)

# Create a MeshFem for u and rhs fields of dimension 1 (i.e. a scalar field)
mfu  = gf.MeshFem(mesh, 1)
mfP0 = gf.MeshFem(mesh, 1)
# Assign the discontinuous P2 fem to all convexes of the both MeshFem
mfu.set_fem(gf.Fem('FEM_PK(2,2)'))
mfP0.set_fem(gf.Fem('FEM_PK(2,0)'))

# Integration method used
mim = gf.MeshIm(mesh, gf.Integ('IM_TRIANGLE(4)'))

# Inner edges for the computation of the normal derivative jump
in_faces = mesh.inner_faces()
INNER_FACES=18
mesh.set_region(INNER_FACES, in_faces)

# Model
md = gf.Model('real')

# Main unknown
md.add_fem_variable('u', mfu)

# Laplacian term on u
md.add_Laplacian_brick(mim, 'u')

# Volumic source term
# md.add_initialized_fem_data('VolumicData', mfrhs, F1)
# md.add_source_term_brick(mim, 'u', 'VolumicData')

# Neumann condition.
md.add_initialized_data('NeumannData', [0.001])
md.add_source_term_brick(mim, 'u', 'NeumannData', BOTTOM_BOUND)

# Dirichlet condition on the left.
if (Dirichlet_with_multipliers):
  md.add_Dirichlet_condition_with_multipliers(mim, 'u', mfu, LEFT_BOUND)
else:
  md.add_Dirichlet_condition_with_penalization(mim, 'u', dirichlet_coefficient,
                                               LEFT_BOUND)

for refiter in range(5):
    
    if (export_mesh):
      mesh.export_to_vtk('mesh%d.vtk' % refiter)
    
    print('\nYou can view the mesh for instance with')
    print('mayavi2 -d mesh%d.vtk -f ExtractEdges -m Surface \n'%refiter)
    
    # Assembly of the linear system and solve.
    md.solve()
    
    # Main unknown
    U = md.variable('u')
    
    # Residual a posteriori estimator
    
    grad_jump = '( (Grad_u-Interpolate(Grad_u,neighbor_element)).Normal )'

    bulkresidual = 'sqr(element_size*Trace(Hess_u))*Test_psi'
    edgeresidual = '0.25*element_size*sqr(%s)*(Test_psi + Interpolate(Test_psi,neighbor_element))'%grad_jump

    ETA1tmp = gf.asm('generic',mim,1,bulkresidual,-1,md,'psi',1,mfP0,np.zeros(mfP0.nbdof()))
    ETA1 = ETA1tmp [ ETA1tmp.size - mfP0.nbdof() : ETA1tmp.size ]
    ETA2tmp = gf.asm('generic',mim,1,edgeresidual,INNER_FACES,md,'psi',1,mfP0,np.zeros(mfP0.nbdof()))
    ETA2 = ETA2tmp [ ETA2tmp.size - mfP0.nbdof() : ETA2tmp.size ]
    ETA  = np.sqrt ( ETA1 + ETA2 )
    
    # Export data
    mfu.export_to_pos('laplacian%d.pos'%refiter, U, 'Computed solution')
    print('You can view the solution with (for example):')
    print('gmsh laplacian%d.pos'%refiter)
    
    mfP0.export_to_pos('eta1_%d.pos'%refiter, ETA1, 'Bulk residual')
    print('You can view eta1 with (for example):')
    print('gmsh eta1_%d.pos'%refiter)
    
    mfP0.export_to_pos('eta2_%d.pos'%refiter, ETA2, 'Edge residual')
    print('You can view eta2 with (for example):')
    print('gmsh eta2_%d.pos'%refiter)
    
    mfu.export_to_vtk('laplacian%d.vtk'%refiter, mfu, U, 'u', mfP0, ETA1, 'eta1', mfP0, ETA2, 'eta2')
    print('mayavi2 -d laplacian%d.vtk -f WarpScalar -m Surface'%refiter)
    
    # Refine the mesh
    dd=mfP0.basic_dof_from_cvid()
    
    ETAElt = np.zeros(dd[0].size)
    for i in range(dd[0].size):
        ETAElt[i] = ETA[ dd[0][i] ] 
    
    mesh.refine(np.where( ETAElt > 0.6*np.max(ETA) ))
    mesh.optimize_structure()
