#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Python GetFEM++ interface
#
# Copyright (C) 2009-2010 Luis Saavedra, Yves Renard.
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
"""  Linear Elastostatic problem with a crack.

  This program is used to check that python-getfem is working. This is
  also a good example of use of GetFEM++.

  $Id$
"""
try:
  import getfem as gf
except ImportError:
  import sys
  sys.path.append('../../src/python/')
  import getfem as gf
else:
  print "module getfem not found!"

import numpy as np
import math

# Parameters:  ######################
nx = 20                             #
                                    #
DIRICHLET  = 101                    #
                                    #
Lambda = 1.25E10 # Lamé coefficient #
Mu = 1.875E10    # Lamé coefficient #
#####################################

# Global Functions: ###############################
ck0 = gf.GlobalFunction('crack',0)                #
ck1 = gf.GlobalFunction('crack',1)                #
ck2 = gf.GlobalFunction('crack',2)                #
ck3 = gf.GlobalFunction('crack',3)                #
                                                  #
coff = gf.GlobalFunction('cutoff',2,0.4,0.01,0.4) #
                                                  #
ckoff0 = ck0*coff                                 #
ckoff1 = ck1*coff                                 #
ckoff2 = ck2*coff                                 #
ckoff3 = ck3*coff                                 #
###################################################

# Mesh in action:
m = gf.Mesh('regular_simplices', np.arange(-0.5,.5+1./nx,1./nx), np.arange(-.5,.5+1./nx,1./nx))
#m = gf.Mesh('import','gmsh','quad.msh')

# boundary set:
m.set_region(DIRICHLET, m.outer_faces())

# MeshFem in action:
mf_pre_u = gf.MeshFem(m)
mf_pre_u.set_fem(gf.Fem('FEM_PK(2,1)'))

# Levelset in action:
ls = gf.LevelSet(m,1,'y','x')

mls = gf.MeshLevelSet(m)
mls.add(ls)
mls.adapt()

# MeshFemLevelSet:
mfls_u = gf.MeshFem('levelset',mls,mf_pre_u)

# MeshFemGlobalFunction:
mf_sing_u = gf.MeshFem('global function',m,ls,[ckoff0,ckoff1,ckoff2,ckoff3])

# MeshImLevelSet:
mim = gf.MeshIm('levelset', mls, 'all',
        gf.Integ('IM_STRUCTURED_COMPOSITE(IM_TRIANGLE(6),3)'),
        gf.Integ('IM_STRUCTURED_COMPOSITE(IM_GAUSS_PARALLELEPIPED(2,6),9)'),
        gf.Integ('IM_STRUCTURED_COMPOSITE(IM_TRIANGLE(6),5)'))

# MeshFemDirectSum:
mf_u = gf.MeshFem('sum',mf_sing_u,mfls_u)
mf_u.set_qdim(2)

# exact solution:
mf_ue = gf.MeshFem('global function',m,ls,[ck0,ck1,ck2,ck3])

A = 2+2*Mu/(Lambda+2*Mu); B=-2*(Lambda+Mu)/(Lambda+2*Mu)
Ue = np.zeros([2,4])
Ue[0,0] =   0; Ue[1,0] = A-B # sin(theta/2)
Ue[0,1] = A+B; Ue[1,1] = 0   # cos(theta/2)
Ue[0,2] =  -B; Ue[1,2] = 0   # sin(theta/2)*sin(theta)
Ue[0,3] =   0; Ue[1,3] = B   # cos(theta/2)*cos(theta)
Ue /= 2*np.pi
Ue = Ue.T.reshape(1,8)

# Model in action:
md = gf.Model('real')
md.add_fem_variable('u', mf_u)

# data
md.add_initialized_data('lambda', [Lambda])
md.add_initialized_data('mu', [Mu])
md.add_isotropic_linearized_elasticity_brick(mim,'u','lambda','mu')

# il fault!!!
#md.add_variable('mult_spec',6)
#BB = gf.Spmat('empty',6,mf_u.nbdof())
#md.add_constraint_with_multipliers('u','mult_spec',BB,plain_vector(6))

md.add_initialized_fem_data("DirichletData", mf_ue, Ue)
md.add_Dirichlet_condition_with_penalization(mim,'u', 1e12, DIRICHLET, 'DirichletData')

# assembly of the linear system and solve:
md.solve()

U = md.variable('u')



# export to pos
cut_mesh = mls.cut_mesh();
mfv = gf.MeshFem(cut_mesh, 2)
mfv.set_classical_discontinuous_fem(2, 0.001)
mf_ue.set_qdim(2)

V  = gf.compute_interpolate_on(mf_u, U, mfv)
Ve = gf.compute_interpolate_on(mf_ue, Ue, mfv)

mfvm = gf.MeshFem(cut_mesh);
mfvm.set_classical_discontinuous_fem(2, 0.001);
md.add_initialized_fem_data('u_cut', mfv, V);
VM = md.compute_isotropic_linearized_Von_Mises_or_Tresca('u_cut', 'lambda', 'mu', mfvm);

mfv.export_to_pos('crack.pos', V, 'V', Ve, 'Ve', mfvm, VM, 'Von Mises')

print 'You can view the solution with (for example):'
print 'gmsh crack.pos'
