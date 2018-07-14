#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Python GetFEM++ interface
#
# Copyright (C) 2018-2019 Yves Renard.
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
"""  Large deformation thermo-elastostatic problem with a gas-filled crack.

  This program is used to check that python-getfem is working. This is
  also a good example of use of GetFEM++.

  $Id$
"""
import numpy as np
import math
import getfem as gf

# Parameters:  #########################################################
ny = 10                    # Number of element in each direction
E0 = 1E4                   # Young modulus at 20oC N/cm^2
nu = 0.45                  # Poisson ratio
alpha = 10.0               # Thermal exchange coefficient
beta = 10                  # Expansion coefficient
gamma = 0.1                # Coefficient of thermal weakening
theta0 = 30.0              # Reference temperature
Lambda = 10.0              # Diffusion coefficient
NR = 287*0.1               # 287 * mass in kg
precribed_temp_top = 25.0  # temperature prescribed at top boundary
precribed_temp_bottom = 35.0  # temperature prescribed at bottom boundary
quad = True                # quad mesh or triangle one
# Mesh definition:
L = 40.;
l = 10.;
h = l/ny;
nx = np.floor(L/h);

if (quad):
  m=gf.Mesh('cartesian', -L/2.+np.arange(nx+1)*h, -l/2+np.arange(ny+1)*h)
else:
  m=gf.Mesh('regular_simplices', -L/2+np.arange(nx+1)*h, -l/2+np.arange(ny+1)*h)

DIR_BOUND_LEFT   = 101
fleft   = m.outer_faces_with_direction([-1., 0.], 0.5)
m.set_region(DIR_BOUND_LEFT,   fleft);
DIR_BOUND_RIGHT  = 102
fright  = m.outer_faces_with_direction([1., 0.], 0.5)
m.set_region(DIR_BOUND_RIGHT,  fright);
DIR_BOUND_TOP    = 103
ftop    = m.outer_faces_with_direction([0., 1.], 0.5)
m.set_region(DIR_BOUND_TOP,    ftop);
DIR_BOUND_BOTTOM = 104
fbottom = m.outer_faces_with_direction([0., -1.], 0.5)
m.set_region(DIR_BOUND_BOTTOM, fbottom);

m.export_to_vtk('mesh.vtk')

# Levelset definition:
R1 = 2.5; R2 = 16;
ytip = R1
xtip = np.sqrt(R2*R2-R1*R1)
ls1 = gf.LevelSet(m,2,'y-%g*tanh(x/7.)' % R1,'x*x+y*y-%g' % (R2*R2))
ls2 = gf.LevelSet(m,2,'y+%g*tanh(x/7.)' % R1,'x*x+y*y-%g' % (R2*R2))
mls = gf.MeshLevelSet(m)
mls.add(ls1)
mls.add(ls2)
mls.adapt()
mfls = ls1.mf();

# Basic mesh_fem without enrichment:
mf_pre_u = gf.MeshFem(m)
if (quad):
  mf_pre_u.set_fem(gf.Fem('FEM_QK(2,2)'))
else:
  mf_pre_u.set_fem(gf.Fem('FEM_PK(2,2)'))

# Definition of the enriched finite element method (MeshFemLevelSet):
mfls_u = gf.MeshFem('levelset',mls, mf_pre_u)

# Global functions for asymptotic enrichment:
mf_part_unity = gf.MeshFem(m)
mf_part_unity.set_classical_fem(1)
DOFpts = mf_part_unity.basic_dof_nodes()
Idofs_1 = np.nonzero(np.square(DOFpts[0,:]+xtip) +
                     np.square(DOFpts[1,:]-ytip) <= 0.5**2)[0]
Idofs_2 = np.nonzero(np.square(DOFpts[0,:]-xtip) +
                     np.square(DOFpts[1,:]+ytip) <= 0.5**2)[0]
Idofs_3 = np.nonzero(np.square(DOFpts[0,:]+xtip) +
                     np.square(DOFpts[1,:]+ytip) <= 0.5**2)[0]
Idofs_4 = np.nonzero(np.square(DOFpts[0,:]-xtip) +
                     np.square(DOFpts[1,:]-ytip) <= 0.5**2)[0]

ck0 = gf.GlobalFunction('crack',0)
ck1 = gf.GlobalFunction('crack',1)
ck2 = gf.GlobalFunction('crack',2)
ck3 = gf.GlobalFunction('crack',3)
mf_sing_u1 = gf.MeshFem('global function',m,ls1,[ck0,ck1,ck2,ck3],1)
mf_sing_u2 = gf.MeshFem('global function',m,ls2,[ck0,ck1,ck2,ck3],1)
# mf_sing_u1 = gf.MeshFem('global function',m,ls1,[ck0],1)
# mf_sing_u2 = gf.MeshFem('global function',m,ls2,[ck0],1)
mf_xfem_sing1 = gf.MeshFem('product', mf_part_unity, mf_sing_u1)
mf_xfem_sing1.set_enriched_dofs(np.union1d(Idofs_1, Idofs_2))
mf_xfem_sing2 = gf.MeshFem('product', mf_part_unity, mf_sing_u2)
mf_xfem_sing2.set_enriched_dofs(np.union1d(Idofs_3, Idofs_4))
mf_u = gf.MeshFem('sum', mf_xfem_sing1, mf_xfem_sing2, mfls_u)
# mf_u = gf.MeshFem('sum', mfls_u)
mf_u.set_qdim(2)

mf_theta = gf.MeshFem('sum', mf_xfem_sing1, mf_xfem_sing2, mfls_u)
# mf_theta = gf.MeshFem('sum', mfls_u)

# MeshIm definition (MeshImLevelSet):
if (quad):
  mim = gf.MeshIm('levelset', mls, 'all',
        gf.Integ('IM_STRUCTURED_COMPOSITE(IM_TRIANGLE(6),5)'),
        gf.Integ('IM_STRUCTURED_COMPOSITE(IM_GAUSS_PARALLELEPIPED(2,6),9)'),
        gf.Integ('IM_GAUSS_PARALLELEPIPED(2,4)'))
else:
  mim = gf.MeshIm('levelset', mls, 'all',
        gf.Integ('IM_STRUCTURED_COMPOSITE(IM_TRIANGLE(6),3)'),
        gf.Integ('IM_STRUCTURED_COMPOSITE(IM_GAUSS_PARALLELEPIPED(2,6),9)'),
        gf.Integ('IM_STRUCTURED_COMPOSITE(IM_TRIANGLE(6),5)'))

mim_bound1 = gf.MeshIm('levelset', mls, 'boundary(a)',
                       gf.Integ('IM_STRUCTURED_COMPOSITE(IM_TRIANGLE(6),3)'))
mim_bound2 = gf.MeshIm('levelset', mls, 'boundary(b)',
                       gf.Integ('IM_STRUCTURED_COMPOSITE(IM_TRIANGLE(6),3)'))
                      
surf_crack = gf.asm("generic", mim_bound1, 0, '1', -1)+gf.asm("generic", mim_bound2, 0, '1', -1)
print("surf_crack = %g" % surf_crack)


# Model definition:
md = gf.Model('real')
md.add_fem_variable('u', mf_u)
md.add_fem_variable('theta', mf_theta)
md.add_variable('P', 1)
md.add_variable('T', 1)
# Data
md.add_initialized_data('surf_crack', [surf_crack])
md.add_initialized_data('E0', [E0])
md.add_initialized_data('nu', [nu])
md.add_initialized_data('alpha', [alpha])
md.add_initialized_data('beta', [beta])
md.add_initialized_data('gamma', [gamma])
md.add_initialized_data('NR', [NR])
md.add_initialized_data('lambda', [Lambda])
md.add_initialized_data('theta0', [theta0])
# Jump of a variable accross the crack
md.add_macro('jump(a)', 'Xfem_plus(a)-Xfem_minus(a)')
# Lamé coefficient
md.add_macro('mu', "E0/((2*(1+nu))*(1+gamma*(theta-theta0)))")
# Lamé coefficient
md.add_macro('la', "2*mu*nu/(1-2*nu)")
# Deformation gradient
md.add_macro('F', "Id(meshdim)+Grad_u")
md.add_macro('Fm', "Id(meshdim)+Xfem_minus(Grad_u)")
md.add_macro('Fp', "Id(meshdim)+Xfem_plus(Grad_u)")
md.add_macro('invFt', "Inv(F')")
md.add_macro('invFmt', "Inv(Fm')")
md.add_macro('invFpt', "Inv(Fp')")

# Deformation tensor
md.add_macro('E', "(Grad_u+Grad_u'+Grad_u'*Grad_u)*0.5")
# Second Piola-Kirchhoff Stress tensor, elastic part
md.add_macro('S', "(la*Trace(E)*Id(meshdim)+2*mu*E)")
# Elasticity term
md.add_nonlinear_term(mim, "(F*S):Grad_Test_u")
md.add_nonlinear_term(mim, "-1E30*neg_part(Det(F)-0.1)*Id(meshdim):Grad_Test_u")
# Thermal expansion term
md.add_nonlinear_term(mim, "-((beta*(theta-theta0)*Det(F))*F):Grad_Test_u")
# Diffusion term
md.add_nonlinear_term(mim,
            "lambda*(invFt*Grad_theta).(invFt*Grad_Test_theta)*Det(F)")
# Ideal Gas law
md.add_nonlinear_term(mim_bound1, "(P*(abs(jump(u).Normalized(Inv(Fm')*Normal))+1)-NR*T/surf_crack)*Test_P")
md.add_nonlinear_term(mim_bound2, "(P*(abs(jump(u).Normalized(Inv(Fm')*Normal))+1)-NR*T/surf_crack)*Test_P")

# Heat flux equilibrium
md.add_nonlinear_term(mim_bound1,
                      "(Norm(invFmt)*(T-Xfem_minus(theta))+Norm(invFpt)*(T-Xfem_plus(theta)))*Test_T")
md.add_nonlinear_term(mim_bound2,
                      "(Norm(invFmt)*(T-Xfem_minus(theta))+Norm(invFpt)*(T-Xfem_plus(theta)))*Test_T")

# Thermal exchange on the cracks
md.add_nonlinear_term(mim_bound1, "-alpha*(Norm(invFmt)*(T-Xfem_minus(theta))*Xfem_minus(Test_theta) + Norm(invFpt)*(T-Xfem_plus(theta))*Xfem_plus(Test_theta))")
md.add_nonlinear_term(mim_bound2, "-alpha*(Norm(invFmt)*(T-Xfem_minus(theta))*Xfem_minus(Test_theta) + Norm(invFpt)*(T-Xfem_plus(theta))*Xfem_plus(Test_theta))")

# Follower pressure
md.add_nonlinear_term(mim_bound1,
                      "(( P*invFmt*abs(Det(Fm)))*Normal).Xfem_minus(Test_u)")
md.add_nonlinear_term(mim_bound1,
                      "((-P*invFpt*abs(Det(Fp)))*Normal).Xfem_plus(Test_u)")
md.add_nonlinear_term(mim_bound2,
                      "(( P*invFmt*abs(Det(Fm)))*Normal).Xfem_minus(Test_u)")
md.add_nonlinear_term(mim_bound2,
                      "((-P*invFpt*abs(Det(Fp)))*Normal).Xfem_plus(Test_u)")


# Fixed displacement
md.add_Dirichlet_condition_with_multipliers(mim, 'u', 1, DIR_BOUND_LEFT)
md.add_Dirichlet_condition_with_multipliers(mim, 'u', 1, DIR_BOUND_RIGHT)
# Fixed temperature
md.add_initialized_data("precribed_temp_bottom", precribed_temp_bottom)
md.add_initialized_data("precribed_temp_top", precribed_temp_top)
md.add_Dirichlet_condition_with_multipliers(mim, 'theta', 1, DIR_BOUND_BOTTOM,
                                            "precribed_temp_bottom")
md.add_Dirichlet_condition_with_multipliers(mim, 'theta', 1, DIR_BOUND_TOP,
                                            "precribed_temp_top")




# Solve with fixed pressure (to open the crack lips)
md.disable_variable('P');
md.set_variable('P', [15])
md.solve('max_res', 5E-4, 'max_iter', 100, 'noisy')


# Solve fully coupled problem
md.enable_variable('P');
md.solve('max_res', 5E-4, 'max_iter', 100, 'noisy')


# Increase temperature gap
md.set_variable("precribed_temp_top", [24])
md.set_variable("precribed_temp_bottom", [35])
md.solve('max_res', 5E-4, 'max_iter', 100, 'noisy')
md.set_variable("precribed_temp_top", [22])
md.set_variable("precribed_temp_bottom", [35])
md.solve('max_res', 5E-8, 'max_iter', 100, 'noisy')


U = md.variable('u')
theta = md.variable('theta')
T = md.variable('T')
print("Gas temperature %g" % T)
P = md.variable('P')
print("Gas pressure %g" % P)

# Interpolation of the solution on a cut mesh for the drawing purpose
cut_mesh = mls.cut_mesh()
mfv = gf.MeshFem(cut_mesh, 2)
mfv.set_classical_discontinuous_fem(2, 0.001)
mfw = gf.MeshFem(cut_mesh, 1)
mfw.set_classical_discontinuous_fem(2, 0.001)

V  = gf.compute_interpolate_on(mf_u, U, mfv)
Th = gf.compute_interpolate_on(mf_theta, theta, mfw)

# Computation of the Von Mises stress
mfvm = gf.MeshFem(cut_mesh)
mfvm.set_classical_discontinuous_fem(4, 0.001)
md.add_interpolate_transformation_from_expression("id_trans", cut_mesh, m, 'X')
md.add_macro('graduint', 'Interpolate(Grad_u, id_trans)')
md.add_macro('thetaint', 'Interpolate(theta, id_trans)')
md.add_macro('Fint', "Id(meshdim)+graduint")
md.add_macro('Eint', "graduint+graduint'+graduint'*graduint")
md.add_macro('muint', "E0/((2*(1+nu))*(1+gamma*(thetaint-theta0)))")
md.add_macro('laint', "2*muint*nu/(1-2*nu)")
md.add_macro('Sint', "(laint*Trace(Eint)*Id(meshdim)+2*muint*Eint)-(beta*(thetaint-theta0)*Det(Fint))*Id(meshdim)")

VM = md.interpolation("sqrt(3/2)*Norm(Deviator(Cauchy_stress_from_PK2(Sint, graduint)))", mfvm)


mfv.export_to_pos('cracked_body.pos', V, 'V', mfw, Th, 'Temperature', mfvm, VM, 'Von Mises stress')
mfv.export_to_vtk('cracked_body.vtk', V, 'V', mfw, Th, 'Temperature', mfvm, VM, 'Von Mises stress')

print('You can view the solution with (for example):')
print('paraview cracked_body.vtk')
