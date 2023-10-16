#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Python GetFEM++ interface
#
# Copyright (C) 2010-2021 Konstantinos Poulios.
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
"""  This example computes a planetary gear model incorporating different
     contact mechanisms like contact with a rigid obstacle and contact
     between elastic bodies of non matching meshes.

     This program is used to check that python-getfem is working. This is
     also a good example of use of GetFEM++.
"""
import getfem as gf
from math import sin,cos,pi
gf.util_trace_level(1)
gf.util_warning_level(1)

z_1 = 20
z_2 = -64
z_p = 22

a = 99.
R_i = 31.

#rot_angle = 2e-2
torsion = 1000.e3
steps = 10
anglestep = 2*pi/abs(z_2)/float(steps)

Lambda = 1.18e5
Mu = 0.83e5

f_coeff = 0.

disp_fem_order = 2
mult_fem_order = 1 # but discontinuous

integration_degree = 4 # 9 gauss points per quad
integration_contact_degree = 10 # 6 gauss points per face

# mesh import
m_1 = gf.Mesh('import', 'gmsh', './static_contact_planetary_1.msh')
m_2 = gf.Mesh('import', 'gmsh', './static_contact_planetary_2.msh')
m_p1 = gf.Mesh('import', 'gmsh', './static_contact_planetary_3.msh')
m_p2 = gf.Mesh('import', 'gmsh', './static_contact_planetary_4.msh')
m_p3 = gf.Mesh('import', 'gmsh', './static_contact_planetary_5.msh')

# regions definitions for boundary conditions
RG_NEUMANN_1 = 1
RG_DIRICHLET_2 = 2
RG_CONTACT_1 = 3
RG_CONTACT_2 = 4
RG_CONTACT_p1_out = 5
RG_CONTACT_p2_out = 6
RG_CONTACT_p3_out = 7
RG_CONTACT_p1_in = 8
RG_CONTACT_p2_in = 9
RG_CONTACT_p3_in = 10
RG_CONTACT_p1 = 11
RG_CONTACT_p2 = 12
RG_CONTACT_p3 = 13

for i in range(1, z_1 + 1):
   m_1.region_merge(RG_NEUMANN_1, 100043+100*i)
   m_1.region_merge(RG_NEUMANN_1, 100083+100*i)
   m_1.region_merge(RG_CONTACT_1, 100013+100*i)
   m_1.region_merge(RG_CONTACT_1, 100053+100*i)

for i in range(1, abs(z_2) + 1):
   m_2.region_merge(RG_DIRICHLET_2, 200043+100*i)
   m_2.region_merge(RG_DIRICHLET_2, 200083+100*i)
   m_2.region_merge(RG_CONTACT_2, 200013+100*i)
   m_2.region_merge(RG_CONTACT_2, 200053+100*i)

for i in range(1, z_p + 1):
   m_p1.region_merge(RG_CONTACT_p1_in, 300043+100*i)
   m_p1.region_merge(RG_CONTACT_p1_in, 300083+100*i)
   m_p1.region_merge(RG_CONTACT_p1_out, 300013+100*i)
   m_p1.region_merge(RG_CONTACT_p1_out, 300053+100*i)

   m_p2.region_merge(RG_CONTACT_p2_in, 400043+100*i)
   m_p2.region_merge(RG_CONTACT_p2_in, 400083+100*i)
   m_p2.region_merge(RG_CONTACT_p2_out, 400013+100*i)
   m_p2.region_merge(RG_CONTACT_p2_out, 400053+100*i)

   m_p3.region_merge(RG_CONTACT_p3_in, 500043+100*i)
   m_p3.region_merge(RG_CONTACT_p3_in, 500083+100*i)
   m_p3.region_merge(RG_CONTACT_p3_out, 500013+100*i)
   m_p3.region_merge(RG_CONTACT_p3_out, 500053+100*i)

m_p1.region_merge(RG_CONTACT_p1, RG_CONTACT_p1_in)
m_p1.region_merge(RG_CONTACT_p1, RG_CONTACT_p1_out)
m_p2.region_merge(RG_CONTACT_p2, RG_CONTACT_p2_in)
m_p2.region_merge(RG_CONTACT_p2, RG_CONTACT_p2_out)
m_p3.region_merge(RG_CONTACT_p3, RG_CONTACT_p3_in)
m_p3.region_merge(RG_CONTACT_p3, RG_CONTACT_p3_out)

N = m_1.dim()

# displacement meshfems
mfu_1 = gf.MeshFem(m_1, N)
mfu_2 = gf.MeshFem(m_2, N)
mfu_p1 = gf.MeshFem(m_p1, N)
mfu_p2 = gf.MeshFem(m_p2, N)
mfu_p3 = gf.MeshFem(m_p3, N)
mfu_1.set_classical_fem(disp_fem_order)
mfu_2.set_classical_fem(disp_fem_order)
mfu_p1.set_classical_fem(disp_fem_order)
mfu_p2.set_classical_fem(disp_fem_order)
mfu_p3.set_classical_fem(disp_fem_order)

# rhs meshfems
mfout_1 = gf.MeshFem(m_1, 1)
mfout_2 = gf.MeshFem(m_2, 1)
mfout_p1 = gf.MeshFem(m_p1, 1)
mfout_p2 = gf.MeshFem(m_p2, 1)
mfout_p3 = gf.MeshFem(m_p3, 1)
mfout_1.set_classical_discontinuous_fem(disp_fem_order)
mfout_2.set_classical_discontinuous_fem(disp_fem_order)
mfout_p1.set_classical_discontinuous_fem(disp_fem_order)
mfout_p2.set_classical_discontinuous_fem(disp_fem_order)
mfout_p3.set_classical_discontinuous_fem(disp_fem_order)

mfmult_1 = gf.MeshFem(m_1, N)
#mfmult_2 = gf.MeshFem(m_2, N)
mfmult_p1 = gf.MeshFem(m_p1, N)
mfmult_p2 = gf.MeshFem(m_p2, N)
mfmult_p3 = gf.MeshFem(m_p3, N)
mfmult_1.set_classical_discontinuous_fem(mult_fem_order)
#mfmult_2.set_classical_discontinuous_fem(mult_fem_order)
mfmult_p1.set_classical_discontinuous_fem(mult_fem_order)
mfmult_p2.set_classical_discontinuous_fem(mult_fem_order)
mfmult_p3.set_classical_discontinuous_fem(mult_fem_order)
#mfmult_p1.set_classical_fem(mult_fem_order)
#mfmult_p2.set_classical_fem(mult_fem_order)
#mfmult_p3.set_classical_fem(mult_fem_order)

# Integration methods
mim_1 = gf.MeshIm(m_1, integration_degree)
mim_2 = gf.MeshIm(m_2, integration_degree)
mim_p1 = gf.MeshIm(m_p1, integration_degree)
mim_p2 = gf.MeshIm(m_p2, integration_degree)
mim_p3 = gf.MeshIm(m_p3, integration_degree)
mim_contact_1 = gf.MeshIm(m_1, integration_contact_degree)
mim_contact_2 = gf.MeshIm(m_2, integration_contact_degree)
mim_contact_p1 = gf.MeshIm(m_p1, integration_contact_degree)
mim_contact_p2 = gf.MeshIm(m_p2, integration_contact_degree)
mim_contact_p3 = gf.MeshIm(m_p3, integration_contact_degree)


# Model definition
md = gf.Model('real')
md.add_fem_variable('u_1', mfu_1)
md.add_fem_variable('u_2', mfu_2)
md.add_fem_variable('u_p1', mfu_p1)
md.add_fem_variable('u_p2', mfu_p2)
md.add_fem_variable('u_p3', mfu_p3)
w_1_str = w_2_str = w_p1_str = w_p2_str = w_p3_str = ''
if f_coeff > 1e-10:
   w_1_str = 'w_1'
   w_2_str = 'w_2'
   w_p1_str = 'w_p1'
   w_p2_str = 'w_p2'
   w_p3_str = 'w_p3'
   md.add_fem_data(w_1_str, mfu_1)
   md.add_fem_data(w_2_str, mfu_2)
   md.add_fem_data(w_p1_str, mfu_p1)
   md.add_fem_data(w_p2_str, mfu_p2)
   md.add_fem_data(w_p3_str, mfu_p3)

elast_law = 'SaintVenant Kirchhoff'
md.add_initialized_data('elast_params', [Lambda, Mu])
md.add_nonlinear_elasticity_brick(mim_1, 'u_1', elast_law, 'elast_params')
md.add_nonlinear_elasticity_brick(mim_2, 'u_2', elast_law, 'elast_params')
md.add_nonlinear_elasticity_brick(mim_p1, 'u_p1', elast_law, 'elast_params')
md.add_nonlinear_elasticity_brick(mim_p2, 'u_p2', elast_law, 'elast_params')
md.add_nonlinear_elasticity_brick(mim_p3, 'u_p3', elast_law, 'elast_params')

# Dirichlet BC's
F = md.interpolation('[0,0]', mfu_2, RG_DIRICHLET_2)
md.add_initialized_fem_data('dirichlet_2', mfu_2, F)
md.add_Dirichlet_condition_with_multipliers(mim_2, 'u_2', mfu_2, RG_DIRICHLET_2, 'dirichlet_2')

# Load
area_1_in = gf.asm_generic(mim_1, 0, "1", RG_NEUMANN_1)
m1 = torsion / area_1_in
_expr_load_ = "{m}/Norm_sqr(X+u_1)*[-X(2)-u_1(2);X(1)+u_1(1)].Test_u_1".format(m=m1)
md.add_nonlinear_term(mim_1, _expr_load_, RG_NEUMANN_1)


# Add inertia, used temporarily for getting an initial solution
md.add_initialized_data('penalty_param', 1e0)
#ibin_1 = md.add_mass_brick(mim_1, 'u_1', 'penalty_param')
##ibin_2 = md.add_mass_brick(mim_2, 'u_2', 'penalty_param')
#ibin_p1 = md.add_mass_brick(mim_p1, 'u_p1', 'penalty_param')
#ibin_p2 = md.add_mass_brick(mim_p2, 'u_p2', 'penalty_param')
#ibin_p3 = md.add_mass_brick(mim_p3, 'u_p3', 'penalty_param')

ibin_1 = md.add_linear_term(mim_1, 'penalty_param*u_1.Test_u_1')
ibin_p1 = md.add_linear_term(mim_p1, 'penalty_param*u_p1.Test_u_p1')
ibin_p2 = md.add_linear_term(mim_p2, 'penalty_param*u_p2.Test_u_p2')
ibin_p3 = md.add_linear_term(mim_p3, 'penalty_param*u_p3.Test_u_p3')

md.add_filtered_fem_variable('mult_1', mfmult_1, RG_NEUMANN_1)
#md.add_filtered_fem_variable('mult_2', mfmult_2, RG_CONTACT_2)
md.add_filtered_fem_variable('mult_p1', mfmult_p1, RG_CONTACT_p1)
md.add_filtered_fem_variable('mult_p2', mfmult_p2, RG_CONTACT_p2)
md.add_filtered_fem_variable('mult_p3', mfmult_p3, RG_CONTACT_p3)

release_dist = 5.
aug_factor = 0.1;
alpha = 1.
md.add_initialized_data( 'r', aug_factor * Mu * (3*Lambda + 2*Mu) / (Lambda + Mu) )
md.add_initialized_data( 'f_coeff', f_coeff)
md.add_initialized_data('alpha', alpha)
ibc = md.add_integral_large_sliding_contact_brick_raytracing('r', release_dist, 'f_coeff', 'alpha', 0)

md.add_slave_contact_boundary_to_large_sliding_contact_brick\
(ibc, mim_contact_1, RG_NEUMANN_1, 'u_1', 'mult_1', w_1_str)
md.add_slave_contact_boundary_to_large_sliding_contact_brick\
(ibc, mim_contact_p1, RG_CONTACT_p1, 'u_p1', 'mult_p1', w_p1_str)
md.add_slave_contact_boundary_to_large_sliding_contact_brick\
(ibc, mim_contact_p2, RG_CONTACT_p2, 'u_p2', 'mult_p2', w_p2_str)
md.add_slave_contact_boundary_to_large_sliding_contact_brick\
(ibc, mim_contact_p3, RG_CONTACT_p3, 'u_p3', 'mult_p3', w_p3_str)
md.add_master_contact_boundary_to_large_sliding_contact_brick\
(ibc, mim_contact_1, RG_CONTACT_1, 'u_1', w_1_str)
md.add_master_contact_boundary_to_large_sliding_contact_brick\
(ibc, mim_contact_2, RG_CONTACT_2, 'u_2', w_2_str)

bearing_1 = 'Norm(X-[%e;%e])-(%e)' % (0., 0., 25.)
bearing_p1 = 'Norm(X-[%e;%e])-(%e)' % (0., a, R_i)
bearing_p2 = 'Norm(X-[%e;%e])-(%e)' % (a*cos(7*pi/6), a*sin(7*pi/6), R_i)
bearing_p3 = 'Norm(X-[%e;%e])-(%e)' % (a*cos(11*pi/6), a*sin(11*pi/6), R_i)
#md.add_rigid_obstacle_to_large_sliding_contact_brick(ibc, bearing_1, N)
md.add_rigid_obstacle_to_large_sliding_contact_brick(ibc, bearing_p1, N)
md.add_rigid_obstacle_to_large_sliding_contact_brick(ibc, bearing_p2, N)
md.add_rigid_obstacle_to_large_sliding_contact_brick(ibc, bearing_p3, N)

print('nbdof_1', mfu_1.nbdof())
print('nbdof_2', mfu_2.nbdof())
print('nbdof_p1', mfu_p1.nbdof())

for i in range(6):
   print('SOLVING WITH TEMPORARY INERTIA FACTOR=%g' % pow(10.,1-i))
   md.set_variable('penalty_param', [pow(10.,1-i)])
   md.solve('noisy', 'max_iter', 40, 'max_res', 1e-2, #)[0]
            'lsearch', 'simplest', 'alpha max ratio', 10, 'alpha min', 0.6, 'alpha mult', 0.6,
            'alpha threshold res', 1000.)
md.disable_bricks([ibin_1, ibin_p1, ibin_p2, ibin_p3])

for nit in range(steps+1):
   print('SOLVING LOAD STEP %i' % nit)
   if nit > 0:
      rot_angle = nit*anglestep
      cc = cos(rot_angle)
      ss = sin(rot_angle)
      F = md.interpolation('[X(1)*(%e)+X(2)*(%e);X(1)*(%e)+X(2)*(%e)]' % (cc-1.,ss,-ss,cc-1.),
                           mfu_2, RG_DIRICHLET_2)
      md.set_variable('dirichlet_2', F)

   if w_1_str:
      md.set_variable(w_1_str, md.variable("u_1"))
      md.set_variable(w_2_str, md.variable("u_2"))
      md.set_variable(w_p1_str, md.variable("u_p1"))
      md.set_variable(w_p2_str, md.variable("u_p2"))
      md.set_variable(w_p3_str, md.variable("u_p3"))

   md.solve('noisy', 'max_iter', 40, 'max_res', 1e-8, #)[0]
            'lsearch', 'simplest', 'alpha max ratio', 10, 'alpha min', 0.3, 'alpha mult', 0.6,
            'alpha threshold res', 1000.)

   U_1 = md.variable('u_1')
   U_2 = md.variable('u_2')
   U_p1 = md.variable('u_p1')
   U_p2 = md.variable('u_p2')
   U_p3 = md.variable('u_p3')
   VM_1 = md.compute_Von_Mises_or_Tresca('u_1', elast_law, 'elast_params', mfout_1)
   VM_2 = md.compute_Von_Mises_or_Tresca('u_2', elast_law, 'elast_params', mfout_2)
   VM_p1 = md.compute_Von_Mises_or_Tresca('u_p1', elast_law, 'elast_params', mfout_p1)
   VM_p2 = md.compute_Von_Mises_or_Tresca('u_p2', elast_law, 'elast_params', mfout_p2)
   VM_p3 = md.compute_Von_Mises_or_Tresca('u_p3', elast_law, 'elast_params', mfout_p3)

   mfout_1.export_to_vtu('static_contact_planetary_1_%i.vtu' % nit,
                         mfout_1,  VM_1, 'Von Mises Stress', mfu_1, U_1, 'Displacement')
   mfout_2.export_to_vtu('static_contact_planetary_2_%i.vtu' % nit,
                         mfout_2,  VM_2, 'Von Mises Stress', mfu_2, U_2, 'Displacement')
   mfout_p1.export_to_vtu('static_contact_planetary_p1_%i.vtu' % nit,
                          mfout_p1,  VM_p1, 'Von Mises Stress', mfu_p1, U_p1, 'Displacement')
   mfout_p2.export_to_vtu('static_contact_planetary_p2_%i.vtu' % nit,
                          mfout_p2,  VM_p2, 'Von Mises Stress', mfu_p2, U_p2, 'Displacement')
   mfout_p3.export_to_vtu('static_contact_planetary_p3_%i.vtu' % nit,
                          mfout_p3,  VM_p3, 'Von Mises Stress', mfu_p3, U_p3, 'Displacement')
