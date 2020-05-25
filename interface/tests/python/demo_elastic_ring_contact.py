#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Python GetFEM interface
#
# Copyright (C) 2013-2020 Konstantinos Poulios.
#
# This file is a part of GetFEM
#
# GetFEM  is  free software;  you  can  redistribute  it  and/or modify it
# under  the  terms  of the  GNU  Lesser General Public License as published
# by  the  Free Software Foundation;  either version 3 of the License,  or
# (at your option) any later version along with the GCC Runtime Library
# Exception either version 3.1 or (at your option) any later version.
# This program  is  distributed  in  the  hope  that it will be useful,  but
# WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
# or  FITNESS  FOR  A PARTICULAR PURPOSE.  See the GNU Lesser General Public
# License and GCC Runtime Library Exception for more details.
# You  should  have received a copy of the GNU Lesser General Public License
# along  with  this program;  if not, write to the Free Software Foundation,
# Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301, USA.
#
############################################################################

import time

import numpy as np

import getfem as gf

gf.util_trace_level(1)

# Input data
ri = 90.   # ring inner diameter
t1 = 5.    # inner layer thickness
t2 = 5.    # outer layer thickness
Ncirc = 64 # number of elements in the ring circumferential direction
Nt1 = 1    # number of elements in first layer thickness direction
Nt2 = 1    # number of elements in second layer thickness direction
lx = 260.  # length of obstacle block
ly = 50.   # height of obstacle block
Nx = 52    # number of elements in block length direction
Ny = 10    # number of elements in block height direction

E1 = 1e5   # Young's modulus
nu1 = 0.3  # Poisson's ratio
E2 = 1e3   # Young's modulus
nu2 = 0.3  # Poisson's ratio
E = 1e8    # Young's modulus
nu = 0.    # Poisson's ratio

g0 = 20.   # initial gap between ring and block
dg = -0.5  # vertical displacement per load step
steps = 80 # number of load steps

#------------------------------------
geotrans_R = 'GT_QK(2,2)'  # geometric transformation for the ring mesh
geotrans_B = 'GT_QK(2,2)'  # geometric transformation for the block mesh

fem_disp_order_R = 2  # displacements finite element order for the ring
fem_disp_order_B = 2  # displacements finite element order for the block
fem_mult_order_R = 1  # multiplier finite element order for the ring
fem_mult_order_B = 1  # multiplier finite element order for the block

integration_degree_R = 4
integration_degree_B = 4
integration_contact_degree_R = 4
integration_contact_degree_B = 4

is_master_R = False
is_master_B = True

r_aug = 300.             # Augmentation parameter
alpha = 1.               # Alpha coefficient for "sliding velocity"
f_coeff = 0.             # Friction coefficient
release_dist = 5.


#------------------------------------
clambda1 = E1*nu1 / ((1+nu1)*(1-2*nu1))
cmu1 = E1 / (2*(1+nu1))
clambda2 = E2*nu2 / ((1+nu2)*(1-2*nu2))
cmu2 = E2 / (2*(1+nu2))
clambda = E*nu / ((1+nu)*(1-2*nu))
cmu = E / (2*(1+nu))


mesh_R = gf.Mesh('import', 'structured',
                 'GT="%s";ORG=[-1,-1];SIZES=[2,1];NSUBDIV=[%i,%i]'
                 % (geotrans_R, Ncirc, Nt1+Nt2))
mesh_B = gf.Mesh('import', 'structured',
                 'GT="%s";ORG=[%f,%f];SIZES=[%f,%f];NSUBDIV=[%i,%i]'
                 % (geotrans_B, -lx/2, -ly, lx, ly, Nx, Ny))

N = mesh_R.dim()

CONTACT_BOUNDARY_R = 1
DIRICHLET_BOUNDARY_R = 3
CONTACT_BOUNDARY_B = 5
DIRICHLET_BOUNDARY_B = 7
RING1 = 11
RING2 = 12

outer_R = mesh_R.outer_faces()
normals_R = mesh_R.normal_of_faces(outer_R)
contact_boundary_R = outer_R[:,np.nonzero(normals_R[1] < -0.95)[0]]
mesh_R.set_region(CONTACT_BOUNDARY_R, contact_boundary_R)

dirichlet_boundary_R = outer_R[:,np.nonzero(np.absolute(normals_R[0]) > 0.99)[0]]
mesh_R.set_region(DIRICHLET_BOUNDARY_R, dirichlet_boundary_R)

pts_R = mesh_R.pts()
y_interf = -float(Nt1)/float(Nt1+Nt2)

is_in_ring1 = pts_R[1,:] > y_interf-0.001
is_in_ring2 = pts_R[1,:] < y_interf+0.001

cvids = mesh_R.cvid()
(pid,idx) = mesh_R.pid_from_cvid(cvids)
cvs_ring1 = []
cvs_ring2 = []
for i in range(idx.size-1):
   cv = cvids[i]
   if all(is_in_ring1[pid[idx[i]:idx[i+1]]]):
      cvs_ring1.append(cv)
   elif all(is_in_ring2[pid[idx[i]:idx[i+1]]]):
      cvs_ring2.append(cv)
mesh_R.set_region(RING1, np.array(cvs_ring1, ndmin=2))
mesh_R.set_region(RING2, np.array(cvs_ring2, ndmin=2))

for ip in range(pts_R.shape[1]):
   x = pts_R[0,ip]
   y = pts_R[1,ip]
   if y >= y_interf: # ring 1
      y *= t1/(-y_interf);
      r = ri
   else:
      y -= y_interf
      y *= t2/(1+y_interf);
      r = ri + t1
   pts_R[0,ip] = (r-y) * np.sin(np.pi*x/2.)
   pts_R[1,ip] = g0+ri+t1+t2 - (r-y) * np.cos(np.pi*x/2.)
mesh_R.set_pts(pts_R)

outer_B = mesh_B.outer_faces()
normals_B = mesh_B.normal_of_faces(outer_B)
contact_boundary_B = outer_B[:,np.nonzero(normals_B[1] > 0.95)[0]]
dirichlet_boundary_B = outer_B[:,np.nonzero(normals_B[1] < -0.95)[0]]
mesh_B.set_region(CONTACT_BOUNDARY_B, contact_boundary_B)
mesh_B.set_region(DIRICHLET_BOUNDARY_B, dirichlet_boundary_B)

#pts_B = mesh_B.pts()
#for ip in range(pts_B.shape[1]):
#   x = pts_B[0,ip]
#   y = pts_B[1,ip]
#   pts_B[1,ip] = y + 0.02*x**2
#mesh_B.set_pts(pts_B)

#mesh_R.export_to_vtk('/tmp/mesh_R.vtk')
#mesh_B.export_to_vtk('/tmp/mesh_B.vtk')

# Ring
mfu_R = gf.MeshFem(mesh_R, N)
mfu_R.set_classical_fem(fem_disp_order_R)

pre_mflambda_R = gf.MeshFem(mesh_R, N)
pre_mflambda_R.set_classical_fem(fem_mult_order_R)

mfvm_R = gf.MeshFem(mesh_R)
mfvm_R.set_classical_discontinuous_fem(fem_disp_order_R-1)

mim_R = gf.MeshIm(mesh_R, integration_degree_R)
mim_R_contact = gf.MeshIm(mesh_R, integration_contact_degree_R)

# Block
mfu_B = gf.MeshFem(mesh_B, N)
mfu_B.set_classical_fem(fem_disp_order_B)

pre_mflambda_B = gf.MeshFem(mesh_B, N)
pre_mflambda_B.set_classical_fem(fem_mult_order_B)

mfvm_B = gf.MeshFem(mesh_B)
mfvm_B.set_classical_discontinuous_fem(fem_disp_order_B-1)

mim_B = gf.MeshIm(mesh_B, integration_degree_B)
mim_B_contact = gf.MeshIm(mesh_B, integration_contact_degree_B)

# Model
md = gf.Model('real')

md.add_fem_variable('uR', mfu_R)
if is_master_B:
   md.add_filtered_fem_variable('lambda_ring', pre_mflambda_R, CONTACT_BOUNDARY_R)
if f_coeff > 1e-10:
   md.add_fem_data('wR', mfu_R)

#lawname = 'neo Hookean'
#params_R1 = [cmu1/2., clambda1/2+cmu1/3]
#params_R2 = [cmu2/2., clambda2/2+cmu2/3]
#params_B = [cmu/2., clambda/2+cmu/3]

lawname = 'neo Hookean Ciarlet'
params_R1 = [clambda1, cmu1]
params_R2 = [clambda2, cmu2]
params_B = [clambda, cmu]

#lawname = 'Ciarlet Geymonat'
#params_R1 = [clambda1, cmu1, cmu1/2-clambda1/8]
#params_R2 = [clambda2, cmu2, cmu2/2-clambda2/8]
#params_B = [clambda, cmu, cmu/2-clambda/8]

md.add_initialized_data('params_ring1', params_R1)
md.add_initialized_data('params_ring2', params_R2)
md.add_nonlinear_elasticity_brick(mim_R, 'uR', lawname, 'params_ring1', RING1)
md.add_nonlinear_elasticity_brick(mim_R, 'uR', lawname, 'params_ring2', RING2)

md.add_fem_variable('uB', mfu_B)
if is_master_R:
   md.add_filtered_fem_variable('lambda_block', pre_mflambda_B, CONTACT_BOUNDARY_B)
if f_coeff > 1e-10:
   md.add_fem_data('wB', mfu_B)

md.add_initialized_data('params_block', params_B)
md.add_nonlinear_elasticity_brick(mim_B, 'uB', lawname, 'params_block')

md.add_initialized_data('dirichlet_ring', np.zeros(N))
md.add_Dirichlet_condition_with_multipliers(mim_R, 'uR', mfu_R, DIRICHLET_BOUNDARY_R, 'dirichlet_ring')

md.add_initialized_data('dirichlet_block', np.zeros(N))
md.add_Dirichlet_condition_with_multipliers(mim_B, 'uB', mfu_B, DIRICHLET_BOUNDARY_B, 'dirichlet_block')

md.add_initialized_data('r', r_aug)
md.add_initialized_data('alpha', alpha)
md.add_initialized_data('f', f_coeff)
ibc = md.add_integral_large_sliding_contact_brick_raytracing('r', release_dist, 'f', 'alpha', 0)

wR_str = ''
wB_str = ''
if f_coeff > 1e-10:
   wR_str = 'wR'
   wB_str = 'wB'

if not is_master_R:
   md.add_slave_contact_boundary_to_large_sliding_contact_brick\
   (ibc, mim_R_contact, CONTACT_BOUNDARY_R, 'uR', 'lambda_ring', wR_str)
elif is_master_R and is_master_B:
   md.add_master_slave_contact_boundary_to_large_sliding_contact_brick\
   (ibc, mim_R_contact, CONTACT_BOUNDARY_R, 'uR', 'lambda_ring', wR_str)
else:
   md.add_master_contact_boundary_to_large_sliding_contact_brick\
   (ibc, mim_R_contact, CONTACT_BOUNDARY_R, 'uR', wR_str)

if not is_master_B:
   md.add_slave_contact_boundary_to_large_sliding_contact_brick\
   (ibc, mim_B_contact, CONTACT_BOUNDARY_B, 'uB', 'lambda_block', wB_str)
elif is_master_B and is_master_R:
   md.add_master_slave_contact_boundary_to_large_sliding_contact_brick\
   (ibc, mim_B_contact, CONTACT_BOUNDARY_B, 'uB', 'lambda_block', wB_str)
else:
   md.add_master_contact_boundary_to_large_sliding_contact_brick\
   (ibc, mim_B_contact, CONTACT_BOUNDARY_B, 'uB', wB_str)

dirichlet_R = np.zeros(N)
for nit in range(steps+1):

   if nit == 0:
     dirichlet_R[N-1] -= g0
   else:
     dirichlet_R[N-1] += dg
   md.set_variable('dirichlet_ring', dirichlet_R)

   if f_coeff > 1e-10:
      md.set_variable('wR', md.variable('uR'))
      md.set_variable('wB', md.variable('uB'))

   starttime = time.process_time()
   md.solve('noisy', 'max_iter', 40, 'max_res', 1e-8, #)[0]
            'lsearch', 'simplest', 'alpha max ratio', 1.5, 'alpha min', 0.2, 'alpha mult', 0.6)[0]
   print('solution time for iteration %i is %f sec' % (nit, time.process_time()-starttime))

   U_R = md.variable('uR')
   VM_R = md.compute_Von_Mises_or_Tresca('uR', lawname, 'params_ring1', mfvm_R)
   mfvm_R.export_to_vtk('lsc_R_%i.vtk' % nit, mfvm_R,  VM_R,
                        'Von Mises Stresses', mfu_R, U_R, 'Displacements')

   lambda_R = md.variable('lambda_ring')
   mf_lambda_R = md.mesh_fem_of_variable('lambda_ring')
   sl = gf.Slice(('boundary',), mf_lambda_R, CONTACT_BOUNDARY_R)
   sl.export_to_vtk('lsc_R_boundary_%i.vtk' % nit,
                    mfu_R, U_R, 'BDisplacements',
                    mf_lambda_R, lambda_R, 'BMultiplier')

   U_B = md.variable('uB')
   VM_B = md.compute_Von_Mises_or_Tresca('uB', lawname, 'params_block', mfvm_B)
   mfvm_B.export_to_vtk('lsc_B_%i.vtk' % nit, mfvm_B,  VM_B,
                        'Von Mises Stresses', mfu_B, U_B, 'Displacements')

   sl = gf.Slice(('boundary',), mfu_B, CONTACT_BOUNDARY_B)
   sl.export_to_vtk('lsc_B_boundary_%i.vtk' % nit,
                    mfu_B, U_B, 'BDisplacements')
