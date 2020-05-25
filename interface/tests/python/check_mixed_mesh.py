#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Python GetFEM interface
#
# Copyright (C) 2018-2020 Huu Phuc Bui
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
"""  Test of computation of inter-element terms on a mixed 3D mesh with
     hexaedrons, prims, pyramids and tetrahedrons.

  This program is used to check that Python-GetFEM interface, and more
  generally GetFEM are working.

  $Id$
"""

import os

import numpy as np

import getfem as gf

# parameters
E = 1.0e3
Nu = 0.25
Lambda = E*Nu/((1.0+Nu)*(1.0-2.0*Nu))
Mu = E/(2*(1+Nu))


# define constants
NEUMANN_BOUNDARY = 1
NEUMANN_BOUNDARY_NO_LOAD = 2
DIRICHLET_BOUNDARY = 3
OMEGA = 4
INNER_FACES=6

degree = 1

make_check=('srcdir' in os.environ);
filename='../meshes/mixed_mesh.gmf'
if (make_check):
    filename=os.environ['srcdir']+'/'+filename

m = gf.Mesh('load', filename)

#--------------- boundary condtions
# detect some boundary of the mesh

face_left = m.outer_faces_with_direction([0, 0., -1.0], 0.01)
face_right = m.outer_faces_with_direction([0., 0., 1.0], 0.01)

face_top = m.outer_faces_with_direction([0., 1., 0], 0.01)
face_bottom = m.outer_faces_with_direction([0., -1., 0], 0.01)
face_top1 = m.outer_faces_with_direction([1., 0., 0], 0.01)
face_bottom1 = m.outer_faces_with_direction([-1., 0., 0], 0.01)
face_top_bottom =  np.append(np.append(np.append(face_top,face_bottom,axis=1),face_top1 , axis=1), face_bottom1, axis=1)

# create boundary regions
m.set_region(NEUMANN_BOUNDARY,face_right)
m.set_region(NEUMANN_BOUNDARY_NO_LOAD,face_top_bottom)
m.set_region(DIRICHLET_BOUNDARY,face_left)

in_faces = m.inner_faces()
m.set_region(INNER_FACES, in_faces)


listTetra = [];
listHexa = [];
listPrism = [];
listPyramid = [];

for i in range(m.nbcvs()):
    gt = m.geotrans(i) # id of element
    if str(gt[0])=='GT_PK(3,1)':
        #print 'index of TETRA: ', i
        listTetra.append(i)
    elif str(gt[0])=='GT_PYRAMID(1)':
        #print 'index of PYRAMID: ', i
        listPyramid.append(i)
    elif str(gt[0])=='GT_PRISM(3,1)':
        #print 'index of PRISM: ', i
        listPrism.append(i)
    elif str(gt[0])=='GT_QK(3,1)':
        #print 'index of HEXA: ', i
        listHexa.append(i)
    else:
        print('Geometric transformation: ', gt[0])

print('num Tetra: ', len(listTetra))
print('num Hexa: ', len(listHexa))
print('num Hexa: ', len(listHexa))
print('num Hexa: ', len(listHexa))


mfu = gf.MeshFem(m, 3)

mfu.set_fem(gf.Fem('FEM_QK(3,{d})'.format(d=degree)),listHexa)
mfu.set_fem(gf.Fem('FEM_PYRAMID_LAGRANGE({d})'.format(d=degree)),listPyramid)
mfu.set_fem(gf.Fem('FEM_PK_PRISM(3,{d})'.format(d=degree)),listPrism)
mfu.set_fem(gf.Fem('FEM_PK(3,{d})'.format(d=degree)),listTetra)

mim = gf.MeshIm(m, 3)


# Model
md = gf.Model('real')
md.add_fem_variable('u',mfu)
md.add_initialized_data('mu_para', Mu)
md.add_initialized_data('lambda_para', Lambda)
md.add_linear_generic_assembly_brick(mim,"lambda_para*Div_u*Div_Test_u + 2*mu_para*Sym(Grad_u):Grad_Test_u")
md.add_initialized_data('Fdata',[0.0,-1.0, 0.0])
md.add_source_term_brick(mim, 'u', 'Fdata', NEUMANN_BOUNDARY)
md.add_initialized_data('DirichletData', [0, 0, 0])
md.add_Dirichlet_condition_with_simplification('u', DIRICHLET_BOUNDARY,'DirichletData')

md.solve('max_res', 1E-9, 'max_iter', 100, 'noisy', 'lsearch', 'simplest',  'alpha min', 0.8)
U = md.variable('u');


mfer = gf.MeshFem(m,1)
mfer.set_fem(gf.Fem('FEM_PK(3,{d})'.format(d=0)),listTetra)
mfer.set_fem(gf.Fem('FEM_QK(3,{d})'.format(d=0)),listHexa)
mfer.set_fem(gf.Fem('FEM_PYRAMID_LAGRANGE({d})'.format(d=0)),listPyramid)
mfer.set_fem(gf.Fem('FEM_PK_PRISM(3,{d})'.format(d=0)),listPrism)

divsigma = '(lambda_para+ mu_para)*(Hess_u(1,1,:) + Hess_u(2,2,:) + Hess_u(3,3,:) ) + mu_para*(Hess_u:Id( qdim(u) ))' 
# 1a) interior residual
bulkresidual = 'sqr(element_size)*Norm_sqr({divsigma})*Test_psi'.format(divsigma=divsigma)

ETA1tmp = gf.asm_generic(mim,1,bulkresidual,-1
                        ,md
                        ,'psi',True,mfer,np.zeros(mfer.nbdof())) 
ETA1 = ETA1tmp [ ETA1tmp.size - mfer.nbdof() : ETA1tmp.size ]


# 1b) jump at inner faces    
sig_u = "(lambda_para*Trace(Grad_u)*Id(qdim(u)) + mu_para*(Grad_u + Grad_u'))"
grad_u_neighbor = "Interpolate(Grad_u,neighbor_element)"
sig_u_neighbor = "(lambda_para*Trace({Grad_u})*Id(qdim(u)) + mu_para*(({Grad_u}) + ({Grad_u})'))".format(Grad_u=grad_u_neighbor)

stress_jump_inner = "((({sig_u}) - ({sig_u_neighbor}))*Normal )".format(sig_u=sig_u,sig_u_neighbor=sig_u_neighbor)
edgeresidual = "0.5*(element_size*Norm_sqr({stress_jump_inner})*2*0.5*(Test_psi + Interpolate(Test_psi,neighbor_element)))".format(stress_jump_inner=stress_jump_inner)

ETA2tmp = gf.asm_generic(mim,1,edgeresidual,INNER_FACES
                        ,md
                        ,'psi',True,mfer,np.zeros(mfer.nbdof()))
ETA2 = ETA2tmp [ ETA2tmp.size - mfer.nbdof() : ETA2tmp.size ]


# 1c) jump at NEUMANN_BOUNDARY  

g = "[0; -1.0; 0]"
stress_jump_at_Neumann = "(({g}) - ({Lambda}*Trace(Grad_u)*Id(qdim(u)) + {Mu}*(Grad_u+Grad_u'))*Normal )".format(g = g, Lambda= Lambda, Mu = Mu)

edgeresidual_Neumann = "(element_size*Norm_sqr({stress_jump_at_Neumann})*Test_psi)".format(stress_jump_at_Neumann=stress_jump_at_Neumann)

ETA3tmp = gf.asm_generic(mim,1,edgeresidual_Neumann,NEUMANN_BOUNDARY
                        ,md
                        ,'psi',True,mfer,np.zeros(mfer.nbdof()))

ETA3 = ETA3tmp [ ETA3tmp.size - mfer.nbdof() : ETA3tmp.size ]

print('sum(ETA3): ', sum(ETA3))

# 1d) jump at NEUMANN_BOUNDARY_NO_LOAD
g = "[0; 0; 0]"
stress_jump_at_Neumann = "(({g}) - ({Lambda}*Trace(Grad_u)*Id(qdim(u)) + {Mu}*(Grad_u+Grad_u'))*Normal )".format(g = g, Lambda= Lambda, Mu = Mu)

edgeresidual_Neumann = "(element_size*Norm_sqr({stress_jump_at_Neumann})*Test_psi)".format(stress_jump_at_Neumann=stress_jump_at_Neumann)

ETA4tmp = gf.asm_generic(mim,1,edgeresidual_Neumann,NEUMANN_BOUNDARY_NO_LOAD
                        ,md ,'psi',True,mfer,np.zeros(mfer.nbdof()))

ETA4 = ETA4tmp [ ETA4tmp.size - mfer.nbdof() : ETA4tmp.size ]

print('sum(ETA4): ', sum(ETA4))


ETA_square = ETA1 + ETA2 + ETA3 + ETA4 # element wise

ETA  = np.sqrt( ETA_square )# element wise
