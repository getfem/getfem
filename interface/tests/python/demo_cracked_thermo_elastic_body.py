#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Python GetFEM interface
#
# Copyright (C) 2018 Yves Renard, Konstantinos Poulios.
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
"""  Large deformation thermo-elastostatic problem with a gas-filled crack.

  This program is used to check that python-getfem is working. This is
  also a good example of use of GetFEM.

  $Id$
"""
import numpy as np

import getfem as gf

np.set_printoptions(threshold=100000)
gf.util_trace_level(1)

# Parameters:  #########################################################
E0 = 5E2            # [N/mm^2] Young modulus at reference temperature
nu = 0.45           # [-] Poisson ratio
beta = 8e-5         # [1/K] Expansion coefficient
gamma = 0.1         # [1/K] Coefficient of thermal weakening
theta0 = 30.        # [C] Reference temperature
alpha = 0.1         # [mm^2/s] Thermal diffusivity (0.1 for rubber)
alpha_env = 0.02    # [N/mm/s/K] Heat transfer coefficient
alpha_crack = 0.01  # [N/mm/s/K] Heat transfer coefficient
mR = 287e3*5e-9     # [N*mm/K] --> 287000 [N*mm/kg/K] * mass [kg]
room_temp = 30.     # [C] temperature on one side of the block
jump_temp = 100.    # [C] temperature jump between the two sides

# Mesh definition:
L = 40.             # [mm]
l = 10.             # [mm]
ny = 10             # Number of elements in each direction
quad = True         # quad mesh or triangle one

h = l/ny;
nx = np.floor(L/h);

if (quad):
  m=gf.Mesh("cartesian", -L/2.+np.arange(nx+1)*h, -l/2+np.arange(ny+1)*h)
else:
  m=gf.Mesh("regular_simplices", -L/2+np.arange(nx+1)*h, -l/2+np.arange(ny+1)*h)

LEFT_RG   = 101
RIGHT_RG  = 102
BOTTOM_RG = 103
TOP_RG    = 104
fleft   = m.outer_faces_with_direction([-1., 0.], 0.5)
fright  = m.outer_faces_with_direction([1., 0.], 0.5)
fbottom = m.outer_faces_with_direction([0., -1.], 0.5)
ftop    = m.outer_faces_with_direction([0., 1.], 0.5)
m.set_region(LEFT_RG,   fleft);
m.set_region(RIGHT_RG,  fright);
m.set_region(BOTTOM_RG, fbottom);
m.set_region(TOP_RG,    ftop);

m.export_to_vtk("mesh.vtk")

# Levelset definition:
R1 = 2.5; R2 = 16;
ytip = R1
xtip = np.sqrt(R2*R2-R1*R1)
ls1 = gf.LevelSet(m, 2, "y-%g*tanh(x/7.)" % R1, "x*x+y*y-%g" % (R2*R2))
ls2 = gf.LevelSet(m, 2, "y+%g*tanh(x/7.)" % R1, "x*x+y*y-%g" % (R2*R2))
mls = gf.MeshLevelSet(m)
mls.add(ls1)
mls.add(ls2)
mls.adapt()

# Basic mesh_fem without enrichment:
mf_pre = gf.MeshFem(m)
if (quad):
  mf_pre.set_fem(gf.Fem("FEM_QK(2,2)"))
else:
  mf_pre.set_fem(gf.Fem("FEM_PK(2,2)"))

# Definition of the enriched finite element method (MeshFemLevelSet):
mfls = gf.MeshFem("levelset", mls, mf_pre)

# Global functions for asymptotic enrichment:
mf_part_unity = gf.MeshFem(m)
mf_part_unity.set_classical_fem(1)
DOFpts = mf_part_unity.basic_dof_nodes()
ctip_dofs = [np.nonzero(np.linalg.norm(DOFpts-x,axis=0) < 0.5)[0]
             for x in [[[xtip],[-ytip]],[[-xtip],[ytip] ],
                       [[xtip],[ytip]], [[-xtip],[-ytip]]]]
ck = [gf.GlobalFunction("crack",i) for i in range(4)]
mf_sing = [gf.MeshFem("product", mf_part_unity,
                      gf.MeshFem("global function", m, ls, ck, 1))
           for ls in [ls1,ls2]]
mf_sing[0].set_enriched_dofs(np.union1d(ctip_dofs[0],ctip_dofs[1]))
mf_sing[1].set_enriched_dofs(np.union1d(ctip_dofs[2],ctip_dofs[3]))
mf_u = gf.MeshFem("sum", mf_sing[0], mf_sing[1], mfls)
mf_u.set_qdim(2)

mf_theta = gf.MeshFem("sum", mf_sing[0], mf_sing[1], mfls)

# MeshIm definition (MeshImLevelSet):
if (quad):
  mim = gf.MeshIm("levelset", mls, "all",
         gf.Integ("IM_STRUCTURED_COMPOSITE(IM_TRIANGLE(6),5)"),
         gf.Integ("IM_STRUCTURED_COMPOSITE(IM_GAUSS_PARALLELEPIPED(2,6),9)"),
         gf.Integ("IM_GAUSS_PARALLELEPIPED(2,4)"))
else:
  mim = gf.MeshIm("levelset", mls, "all",
         gf.Integ("IM_STRUCTURED_COMPOSITE(IM_TRIANGLE(6),3)"),
         gf.Integ("IM_STRUCTURED_COMPOSITE(IM_GAUSS_PARALLELEPIPED(2,6),9)"),
         gf.Integ("IM_STRUCTURED_COMPOSITE(IM_TRIANGLE(6),5)"))

mim_bound = [gf.MeshIm("levelset", mls, boundary,
                       gf.Integ("IM_STRUCTURED_COMPOSITE(IM_TRIANGLE(6),3)"))
             for boundary in ["boundary(a)", "boundary(b)"]]

surf_crack = gf.asm("generic", mim_bound[0], 0, "1", -1)\
            +gf.asm("generic", mim_bound[1], 0, "1", -1)
print("surf_crack = %g" % surf_crack)


# Model definition:
md = gf.Model("real")
md.add_fem_variable("u", mf_u)
md.add_fem_variable("theta", mf_theta)
md.add_variable("P", 1)
md.add_variable("T", 1); md.set_variable("T", 273+room_temp)
# Data
md.add_initialized_data("surf_crack", surf_crack)
md.add_initialized_data("mu0", E0/(2*(1+nu)))
md.add_initialized_data("la0", E0*nu/((1-2*nu)*(1+nu)))
md.add_initialized_data("ka0", E0/(3*(1-2*nu)))
md.add_initialized_data("alpha", alpha)
md.add_initialized_data("alpha_env", alpha_env)
md.add_initialized_data("alpha_crack", alpha_crack)
md.add_initialized_data("beta", beta)
md.add_initialized_data("gamma", gamma)
md.add_initialized_data("mR", mR)
md.add_initialized_data("theta0", theta0)
# Lamé coefficients + bulk modulus
for p1,p2 in [["la","la0"],["mu","mu0"],["ka","ka0"]]:
  md.add_macro(p1, "%s/(1+gamma*(theta-theta0))" % p2)
# Deformation gradient
md.add_macro("F", "Id(2)+Grad_u")
md.add_macro("Fm", "Id(2)+Xfem_minus(Grad_u)")
md.add_macro("Fp", "Id(2)+Xfem_plus(Grad_u)")

md.add_macro("Nanson(FF)", "Det(FF)*Norm(Inv(FF')*Normal)")
md.add_macro("th", "exp(beta*(theta-theta0))")

# Deformation tensor
md.add_macro("E", "(Grad_u+Grad_u'+Grad_u'*Grad_u)*0.5")
# Second Piola-Kirchhoff Stress tensor, elastic + thermal expansion part
md.add_macro("S", "((la*Trace(E)*Id(2)+2*mu*E)/th-1.5*(th-1/th)*ka*Id(2))")
# Elasticity term
md.add_nonlinear_term(mim, "(F*S):Grad_Test_u")
# Diffusion term
md.add_nonlinear_term(mim,
  "alpha*Det(F)*(Inv(F')*Grad_theta).(Inv(F')*Grad_Test_theta)")
for mimb in mim_bound:
  md.add_nonlinear_term(mimb, # Ideal Gas law
  "(P*((Xfem_plus(u)-Xfem_minus(u)).Normalized(Inv(Fm')*Normal)+1e-6)"
  "-mR*T/surf_crack)*Test_P")
  md.add_nonlinear_term(mimb, # Heat flux equilibrium
                        "(Nanson(Fm)*(T-273-Xfem_minus(theta))"
                        "+Nanson(Fp)*(T-273-Xfem_plus(theta)))*Test_T")
  md.add_nonlinear_term(mimb, # Heat exchange on the cracks
  "-alpha_crack*(Nanson(Fm)*(T-273-Xfem_minus(theta))*Xfem_minus(Test_theta)"
               "+Nanson(Fp)*(T-273-Xfem_plus(theta))*Xfem_plus(Test_theta))")
for var,rg in [["theta_B", BOTTOM_RG],["theta_T",TOP_RG]]: # Heat exchange at
  md.add_initialized_data(var, room_temp)                  # bottom/top
  md.add_nonlinear_term(mim,"-alpha_env*Nanson(F)*("+var+"-theta)*Test_theta",
                        rg)
for mimb in mim_bound: # Follower pressure
  md.add_nonlinear_term(mimb,"(P*Det(Fm)*Inv(Fm')*Normal).Xfem_minus(Test_u)")
  md.add_nonlinear_term(mimb,"(-P*Det(Fp)*Inv(Fp')*Normal).Xfem_plus(Test_u)")

# Fixed displacement
md.add_Dirichlet_condition_with_multipliers(mim, "u", 1, LEFT_RG)
md.add_Dirichlet_condition_with_multipliers(mim, "u", 1, RIGHT_RG)


# Solve with fixed pressure (to open the crack lips)
md.disable_variable("P");
md.set_variable("P", [1e-2])
md.solve("max_res", 5E-4, "max_iter", 100, "noisy")

# Solve fully coupled problem
md.enable_variable("P");

# Interpolation of the solution on a cut mesh for the drawing purpose
cut_mesh = mls.cut_mesh()
mfv = gf.MeshFem(cut_mesh, 2)
mfv.set_classical_discontinuous_fem(2, 0.001)
mfw = gf.MeshFem(cut_mesh, 1)
mfw.set_classical_discontinuous_fem(2, 0.001)

# For the computation of the Von Mises stress
mfvm = gf.MeshFem(cut_mesh)
mfvm.set_classical_discontinuous_fem(4, 0.001)
md.add_interpolate_transformation_from_expression("id_trans", cut_mesh, m, "X")
md.add_macro("graduI", "Interpolate(Grad_u, id_trans)")
md.add_macro("thetaI", "Interpolate(theta, id_trans)")
md.add_macro("FI", "Id(2)+graduI")
md.add_macro("EI", "graduI+graduI'+graduI'*graduI")
md.add_macro("muI", "mu0/(1+gamma*(thetaI-theta0))")
md.add_macro("laI", "la0/(1+gamma*(thetaI-theta0))")
md.add_macro("kaI", "ka0/(1+gamma*(thetaI-theta0))")
md.add_macro("thI", "exp(beta*(thetaI-theta0))")
md.add_macro("SI", "((laI*Trace(EI)*Id(2)+2*muI*EI)/thI-1.5*(thI-1/thI)*kaI*Id(2))")

# Increase temperature gap
it = 0
for fact in [0,0.2,0.4,0.6,0.8,1.]:
  md.set_variable("theta_T", room_temp + fact*jump_temp)
  md.solve("max_res", 5E-8, "max_iter", 100, "noisy")

  U = md.variable("u")
  theta = md.variable("theta")
  T = md.variable("T")
  print("Gas temperature %g C" % (T-273))
  P = md.variable("P")
  print("Gas pressure %g MPa" % P)

  V  = gf.compute_interpolate_on(mf_u, U, mfv)
  Th = gf.compute_interpolate_on(mf_theta, theta, mfw)
  VM = md.interpolation("sqrt(3/2)*Norm(Deviator(Cauchy_stress_from_PK2(SI, graduI)))", mfvm)

  mfv.export_to_pos("cracked_body_%i.pos" % it, V, "V", mfw, Th, "Temperature", mfvm, VM, "Von Mises stress")
  mfv.export_to_vtk("cracked_body_%i.vtk" % it, V, "V", mfw, Th, "Temperature", mfvm, VM, "Von Mises stress")
  it +=1

print("You can view the solution with (for example):")
print("paraview cracked_body_2.vtk")
