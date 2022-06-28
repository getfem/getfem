#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Python GetFEM interface
#
# Copyright (C) 2022 Yves Renard.
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

""" Incompressible Navier-Stokes fluid in interaction with a ball in a cavity.
  Middle point scheme for the fluid, Verlet's scheme for the ball (not the
  best but can of course be changed)

  This program is used to check that python-getfem is working. This is also
  a good example of use of GetFEM.
"""

import os
import numpy as np
import getfem as gf
gf.util('trace level', 0)

version = 1   # 1-without cut-fem and Nitsche's method,
              # 2-with cut-fem

# Phisical parameters (à ajuster pour être physique ...)
nu = 0.002        # cinematic viscosity
rho = 1.          # fluid density      
mu = rho * nu     # dynamic viscosity
g = 9.81          # gravity
in_velocity = 4.  # inward velocity at the center bottom
ball_radius = 0.1
ball_mass = 0.032
ball_init_pos = np.array([0., 0.2])

# Discretization parameters
dt = 0.01       # Time step
T = 40.         # Final time
gamma0 = 1.     # Nitsche's method parameter

# Geometry and Mesh of the cavity
W1 = 0.05         #       ____________
W2 = 0.475        #   H2 |_          _|
W = W1+2.*W2      #      |            |
H1 = 0.8          #      |            |
H2 = 0.2          #   H1 |            |
H = H2+H1         #      |            |
                  #      |____|__|____|
                  #        W2  W1  W2
NH = 40
NW = 40
DH = H / NH; DW = W / NW;
print("Number of elements: ", W1/DW, W2/DW, H1/DH, H2/DH)
if (abs(W1 / DW - round(W1/DW)) +  abs(W2 / DW - round(W2/DW))
    + abs(H1 / DH - round(H1/DH)) + abs(H2 / DH != round(H2/DH)) > 1E-8) :
   print("Discretization error"); exit(1)
   
out_velocity = W1*in_velocity / (2.*H2); 
Reynold = rho*in_velocity*H/mu
print("Reynold = ", Reynold)

geotrans = "GT_QK(2,2)"
m = gf.Mesh("import", "structured",
            "GT='%s';ORG=[%f,%f];SIZES=[%f,%f];NSUBDIV=[%i,%i]"
            % (geotrans, -W/2., 0., W, H, NW, NH))

# Mesh regions
IN_RG = 1
OUT_RG1 = 2
OUT_RG2 = 3
WALL_RG = 4
INTERNAL_EDGES = 5;
in_rg   = m.outer_faces_in_box([-W1/2.-1e-8,-1e-8],[W1/2+1e-8,1e-8])
out_rg1 = m.outer_faces_in_box([-W/2.-1e-8,H1-1e-8],[-W/2+1e-8,H+1e-8])
out_rg2 = m.outer_faces_in_box([W/2.-1e-8,H1-1e-8],[W/2+1e-8,H+1e-8])
wall_rg = m.outer_faces()
m.set_region(IN_RG, in_rg)
m.set_region(OUT_RG1, out_rg1)
m.set_region(OUT_RG2, out_rg2)
m.set_region(WALL_RG, wall_rg)
m.region_subtract(WALL_RG, IN_RG)
m.region_subtract(WALL_RG, OUT_RG1)
m.region_subtract(WALL_RG, OUT_RG2)

innerf = m.inner_faces()
m.set_region(INTERNAL_EDGES, innerf)

# Finite element spaces and integration methods
mfv = gf.MeshFem(m, 2)
mfv.set_classical_fem(2)
nbdofv = mfv.nb_basic_dof()

mfp = gf.MeshFem(m, 1)
mfp.set_classical_fem(1)
nbdofp = mfp.nb_basic_dof()

mim = gf.MeshIm(m, 4)

# Levelset definition and adapted integration methods

ls = gf.LevelSet(m,2)
mf_ls = ls.mf()
P = mf_ls.basic_dof_nodes(); x = P[0,:]; y = P[1,:]
ULS = ((x - ball_init_pos[0])**2 + (y - ball_init_pos[1])**2) - ball_radius**2
ls.set_values(ULS)
mls = gf.MeshLevelSet(m)
mls.add(ls)
mls.adapt()
mim_bound = gf.MeshIm('levelset',mls,'boundary', gf.Integ('IM_TRIANGLE(6)')) #, gf.Integ('IM_QUAD(5)'))
mim_in = gf.MeshIm('levelset',mls,'inside', gf.Integ('IM_TRIANGLE(6)'))
mim_in.set_integ(4)
mim_out = gf.MeshIm('levelset',mls,'outside', gf.Integ('IM_TRIANGLE(6)'))
mim_out.set_integ(4)

# Model definition

md = gf.Model("real")
md.add_fem_variable("v", mfv)
md.add_fem_data("v0", mfv)
md.add_fem_data("ls", mf_ls)
md.add_fem_variable("p", mfp)
md.add_fem_data("p_in", mfp)
md.add_initialized_data("f", [0., -rho*g])
md.add_initialized_data("ball_v", [0., 0.])
md.add_initialized_data("v_in", [0., in_velocity])
md.add_initialized_data("v_out1", [-out_velocity, 0.])
md.add_initialized_data("v_out2", [out_velocity, 0.])
md.add_initialized_data("rho", rho)
md.add_initialized_data("gamma", gamma0/DW)
md.add_initialized_data("dt", [dt])
md.add_initialized_data("mu", [mu])

md.add_Dirichlet_condition_with_multipliers(mim, "p", mfp, IN_RG)
md.add_Dirichlet_condition_with_multipliers(mim, "v", mfv, IN_RG, "v_in")
md.add_Dirichlet_condition_with_multipliers(mim, "v", mfv, OUT_RG1, "v_out1")
md.add_Dirichlet_condition_with_multipliers(mim, "v", mfv, OUT_RG2, "v_out2")
md.add_Dirichlet_condition_with_multipliers(mim, "v", mfv, WALL_RG)

# Diffusive terms
if (version == 2):
   mim_t = mim_out;
else:
   mim_t = mim
md.add_nonlinear_term(mim_t, "(1/dt)*rho*(v-v0).Test_v + mu*Grad_v:Grad_Test_v")
# Nonlinear convective term
md.add_nonlinear_term(mim_t, "0.25*rho*((Grad_v+Grad_v0)*(v+v0)).Test_v")
# Pressure terms
md.add_nonlinear_term(mim_t, "-p*Div_Test_v - 0.5*Test_p*(Div_v+Div_v0)")
# Gravity term
md.add_nonlinear_term(mim_t, "-f.Test_v")
# Small ghost penalty term
if (version == 2):
  md.add_linear_term(mim, "1E-6*(Grad_v-Interpolate(Grad_v, neighbor_element)):(Grad_Test_v-Interpolate(Grad_Test_v, neighbor_element))", INTERNAL_EDGES)
  # md.add_linear_term(mim, "1E-7*(p-Interpolate(p, neighbor_element))*(Test_p-Interpolate(Test_p, neighbor_element))", INTERNAL_EDGES)
  # Penalty term on internal dofs
  Bv = gf.Spmat('empty',nbdofv)
  ibv = md.add_explicit_matrix("v", "v", Bv)
  Bp = gf.Spmat('empty',nbdofp)
  ibp = md.add_explicit_matrix("p", "p", Bp)
  # md.add_linear_term(mim_in, "1E-4*(v - ball_v).Test_v")
# Nitsche's term on the FS interface
if ((version == 1) or (version == 2)):
  md.add_nonlinear_term(mim_bound, "-(mu*Grad_v-p*Id(meshdim))*Normalized(Grad_ls).Test_v + gamma * (v-ball_v).Test_v - (mu*Grad_Test_v)*Normalized(Grad_ls).(v-ball_v)")

t = 0
step = 0
ball_pos = ball_init_pos
ball_pos_prec = ball_init_pos
ball_v = np.array([0., 0.])
os.system('mkdir -p FSI_results');
while t < T+1e-8:
   print("Solving step at t=%f" % t)
   md.set_variable("v0", md.variable("v"))
   md.set_variable("ball_v", ball_v)
   
   # levelset update
   P = mf_ls.basic_dof_nodes(); x = P[0,:]; y = P[1,:]
   ULS = ((x - ball_pos[0])**2 + (y - ball_pos[1])**2) - ball_radius**2
   ls.set_values(ULS)
   md.set_variable('ls', ULS)
   mls.adapt()
   mim_bound.adapt()
   mim_in.adapt()
   mim_out.adapt()

   # Penalization of ball internal dofs with no contribution on the boundary
   if (version == 2):
     idofv = np.setdiff1d(np.arange(nbdofv), mfv.dof_from_im(mim_out))
     idofp = np.setdiff1d(np.arange(nbdofp), mfp.dof_from_im(mim_out))
     Bv = gf.Spmat('empty', nbdofv)
     for i in idofv: Bv.add(i,i,1.)
     md.set_private_matrix(ibv, Bv)
     Bp = gf.Spmat('empty', nbdofp)
     for i in idofp: Bp.add(i,i,1.)
     md.set_private_matrix(ibp, Bp)
   
   # Solve
   # md.solve("noisy", "lsolver", "mumps", "max_res", 1e-8)
   md.solve("max_res", 1e-8, "max_iter", 25)

   # Balance of forces on the ball and Verlet's scheme
   R = gf.asm('generic', mim_bound, 0,
              '(2*mu*Sym(Grad_v)-p*Id(meshdim))*Normalized(Grad_ls)', -1, md)
   # R = gf.asm('generic', mim_bound, 0, 'Normalized(Grad_ls)', -1, md)
   ball_pos_next = 2*ball_pos - ball_pos_prec + dt*dt*(R/ball_mass - [0, g])
   ball_v = (ball_pos_next - ball_pos_prec) / (2*dt)
   ball_pos_prec = ball_pos
   ball_pos = ball_pos_next
   
   # Enforce the ball to remain inside the cavity
   if (ball_pos[0] < -W/2+ball_radius) :
      ball_pos_prec[0] = ball_pos[0] = -W/2+ball_radius; ball_v[0] *= 0.
   if (ball_pos[0] >  W/2-ball_radius) :
      ball_pos_prec[0] = ball_pos[0] =  W/2-ball_radius; ball_v[0] *= 0.
   if (ball_pos[1] <  ball_radius)     :
      ball_pos_prec[1] = ball_pos[1] =  ball_radius; ball_v[1] *= 0.
   if (ball_pos[1] >  H-ball_radius)   :
      ball_pos_prec[1] = ball_pos[1] =  H-ball_radius; ball_v[1] *= 0.
   print ("ball position = ", ball_pos, " ball velocity = ", ball_v)

   # Post-processing
   cut_m = mls.cut_mesh();
   mfd = gf.MeshFem(cut_m, 1)
   mfd.set_classical_discontinuous_fem(0)
   P = mfd.basic_dof_nodes(); x = P[0,:]; y = P[1,:]
   UB = (((x - ball_pos_prec[0])**2 + (y - ball_pos_prec[1])**2)
         - ball_radius**2) >= 0.
   mfd.export_to_vtk("FSI_results/FSI_%i.vtk" % step,
                      mfv, md.variable("v"), "Velocity",
                      mfp, md.variable("p"), "Pressure",
                      mfd, UB, "Ball position")
   t += dt
   step += 1

print("See for instance with ")
print("paraview FSI_results/FSI_..vtk")
print("You can add Glyph for velocity.")
print("The field 'Ball position' allows to locate the ball.")
