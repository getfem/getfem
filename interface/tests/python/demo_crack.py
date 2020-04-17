#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Python GetFEM interface
#
# Copyright (C) 2009-2020 Luis Saavedra, Yves Renard.
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
"""  Linear Elastostatic problem with a crack.

  This program is used to check that python-getfem is working. This is
  also a good example of use of GetFEM.

  $Id$
"""
try:
  import getfem as gf
except ImportError:
  import sys
  sys.path.append('../../src/python/')
  import getfem as gf
else:
  print("module getfem not found!")

import numpy as np

variant = 4
# variant : 1 : a single crack with cutoff enrichement
#           2 : a single crack with a fixed size area Xfem enrichment
#           3 : a supplementary crossing  crack with a fixed size area
#               Xfem enrichment
#           4 : variant 3 with the second crack closed by a penalisation of
#               the jump (example of use of xfem_plus and xfem_minus).
#           5 : variant 3 with the first crack closed by a penalisation of
#               the jump (example of use of xfem_plus and xfem_minus).
#           6 : variant 3 with the two cracks closed by a penalisation of
#               the jump (example of use of xfem_plus and xfem_minus).

# Parameters:  ######################
nx = 20                             #
                                    #
DIRICHLET  = 101                    #
                                    #
Lambda = 1.25E10 # Lame coefficient #
Mu = 1.875E10    # Lame coefficient #
#####################################

# Mesh definition:
m = gf.Mesh('regular_simplices', -0.5+np.arange(nx+1)/float(nx),
                                 -0.5+np.arange(nx+1)/float(nx))
#m = gf.Mesh('import','gmsh','quad.msh')

# Boundary set:
m.set_region(DIRICHLET, m.outer_faces())

# Global functions for asymptotic enrichment:
ck0 = gf.GlobalFunction('crack',0)
ck1 = gf.GlobalFunction('crack',1)
ck2 = gf.GlobalFunction('crack',2)
ck3 = gf.GlobalFunction('crack',3)
if variant == 1: # Cutoff enrichement 
   coff = gf.GlobalFunction('cutoff',2,0.4,0.01,0.4)
   ckoff0 = ck0*coff # gf.GlobalFunction('product',ck0,coff)
   ckoff1 = ck1*coff
   ckoff2 = ck2*coff
   ckoff3 = ck3*coff

# Levelset definition:
ls = gf.LevelSet(m,1,'y','x')

mls = gf.MeshLevelSet(m)
mls.add(ls)
if variant > 2:
   ls2 =  gf.LevelSet(m,1,'x+0.125','abs(y)-0.375')
   mls.add(ls2);
mls.adapt()

# Basic mesh_fem without enrichment:
mf_pre_u = gf.MeshFem(m)
mf_pre_u.set_fem(gf.Fem('FEM_PK(2,1)'))

# Definition of the enriched finite element method (MeshFemLevelSet):
mfls_u = gf.MeshFem('levelset',mls,mf_pre_u)

if variant == 1: # Cutoff enrichement 
   # MeshFemGlobalFunction:
   mf_sing_u = gf.MeshFem('global function',m,ls,[ckoff0,ckoff1,ckoff2,ckoff3],1)
   # MeshFemDirectSum:
   mf_u      = gf.MeshFem('sum',mf_sing_u,mfls_u)
else:
   mf_sing_u = gf.MeshFem('global function',m,ls,[ck0,ck1,ck2,ck3],1)
   mf_part_unity = gf.MeshFem(m)
   mf_part_unity.set_classical_fem(1)
   DOFpts = mf_part_unity.basic_dof_nodes()
   # Search the dofs to be enriched with the asymptotic displacement.
   Idofs_center = np.nonzero(np.square(DOFpts[0,:]) +
                             np.square(DOFpts[1,:]) <= 0.1**2)[0]
   mf_xfem_sing = gf.MeshFem('product', mf_part_unity, mf_sing_u)
   mf_xfem_sing.set_enriched_dofs(Idofs_center)
   if variant > 2:
      Idofs_up = np.nonzero(np.square(DOFpts[0,:]+0.125) +
                            np.square(DOFpts[1,:]-0.375) <= 0.1**2)[0]
      Idofs_down = np.nonzero(np.square(DOFpts[0,:]+0.125) +
                              np.square(DOFpts[1,:]+0.375) <= 0.1**2)[0]
      mf_sing_u2 = gf.MeshFem('global function',m,ls2,[ck0,ck1,ck2,ck3],1)
      mf_xfem_sing2 = gf.MeshFem('product', mf_part_unity, mf_sing_u2)
      mf_xfem_sing2.set_enriched_dofs(np.union1d(Idofs_up, Idofs_down))
   if variant == 2:
      mf_u = gf.MeshFem('sum', mf_xfem_sing, mfls_u)
   else:
      mf_u = gf.MeshFem('sum', mf_xfem_sing, mf_xfem_sing2, mfls_u)

mf_u.set_qdim(2)

# MeshIm definition (MeshImLevelSet):
mim = gf.MeshIm('levelset', mls, 'all',
        gf.Integ('IM_STRUCTURED_COMPOSITE(IM_TRIANGLE(6),3)'),
        gf.Integ('IM_STRUCTURED_COMPOSITE(IM_GAUSS_PARALLELEPIPED(2,6),9)'),
        gf.Integ('IM_STRUCTURED_COMPOSITE(IM_TRIANGLE(6),5)'))

# Exact solution for a single crack:
mf_ue = gf.MeshFem('global function',m,ls,[ck0,ck1,ck2,ck3])
A = 2+2*Mu/(Lambda+2*Mu); B=-2*(Lambda+Mu)/(Lambda+2*Mu)
Ue = np.zeros([2,4])
Ue[0,0] =   0; Ue[1,0] = A-B # sin(theta/2)
Ue[0,1] = A+B; Ue[1,1] = 0   # cos(theta/2)
Ue[0,2] =  -B; Ue[1,2] = 0   # sin(theta/2)*sin(theta)
Ue[0,3] =   0; Ue[1,3] = B   # cos(theta/2)*cos(theta)
Ue /= 2*np.pi
Ue = Ue.T.reshape(1,8)

# Model definition:
md = gf.Model('real')
md.add_fem_variable('u', mf_u)
# data
md.add_initialized_data('lambda', [Lambda])
md.add_initialized_data('mu', [Mu])
md.add_isotropic_linearized_elasticity_brick(mim,'u','lambda','mu')
md.add_initialized_fem_data("DirichletData", mf_ue, Ue)
md.add_Dirichlet_condition_with_penalization(mim,'u', 1e12, DIRICHLET, 'DirichletData')

if variant == 5 or variant == 6: # Penalisation of the jump over the first crack
   mim_bound1 = gf.MeshIm('levelset', mls, 'boundary(a)',
                          gf.Integ('IM_STRUCTURED_COMPOSITE(IM_TRIANGLE(6),3)'))
   #gf.asm_generic(mim_bound1, 0, '1', -1) # length of the crack
   md.add_linear_term\
   (mim_bound1, '1e17*(Xfem_plus(u)-Xfem_minus(u)).(Xfem_plus(Test_u)-Xfem_minus(Test_u))')

if variant == 4 or variant == 6: # Penalisation of the jump over the second crack
   mim_bound2 = gf.MeshIm('levelset', mls, 'boundary(b)',
                          gf.Integ('IM_STRUCTURED_COMPOSITE(IM_TRIANGLE(6),3)'))
   md.add_linear_term\
   (mim_bound2, '1e17*(Xfem_plus(u)-Xfem_minus(u)).(Xfem_plus(Test_u)-Xfem_minus(Test_u))')

# Assembly of the linear system and solve:
md.solve()
U = md.variable('u')

# Interpolation of the solution on a cut mesh for the drawing purpose
cut_mesh = mls.cut_mesh()
mfv = gf.MeshFem(cut_mesh, 2)
mfv.set_classical_discontinuous_fem(2, 0.001)
mf_ue.set_qdim(2)

V  = gf.compute_interpolate_on(mf_u, U, mfv)
Ve = gf.compute_interpolate_on(mf_ue, Ue, mfv)

# Computation of the Von Mises stress
mfvm = gf.MeshFem(cut_mesh)
mfvm.set_classical_discontinuous_fem(2, 0.001)
md.add_initialized_fem_data('u_cut', mfv, V)
VM = md.compute_isotropic_linearized_Von_Mises_or_Tresca('u_cut', 'lambda', 'mu', mfvm)

mfv.export_to_pos('crack.pos', V, 'V', Ve, 'Ve', mfvm, VM, 'Von Mises')

print('You can view the solution with (for example):')
print('gmsh crack.pos')
