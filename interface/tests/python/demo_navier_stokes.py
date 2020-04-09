#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Python GetFEM interface
#
# Copyright (C) 2015 Konstantinos Poulios.
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

""" Incompressible Navier-Stokes equation solved in an L-shaped domain.

  This program is used to check that python-getfem is working. This is also
  a good example of use of GetFEM.
"""

# This example is ported from:
# http://fenicsproject.org/documentation/dolfin/1.0.1/python/demo/pde/navier-stokes/python/documentation.html

import numpy as np

import getfem as gf

W1 = 0.5         #       _____
H1 = 0.5         #      |  W1 |
W2 = 0.5         #      |     |H1
H2 = 0.5         #      |     |_____
NW1 = 10         # H1+H2|        W2 |
NH1 = 10         #      |           |H2
NW2 = 10         #      |___________|
NH2 = 10         #          W1+W2

dt = 0.01
T = 3.
nu = 0.01

p_in_str = "np.sin(3*{0})"

#Mesh and MeshRegion
IN_RG = 1
OUT_RG = 2
INOUT_RG = 3
WALL_RG = 4

geotrans = "GT_QK(2,2)"
m = gf.Mesh("import", "structured",
            "GT='%s';ORG=[%f,%f];SIZES=[%f,%f];NSUBDIV=[%i,%i]"
            % (geotrans, 0, H2, W1, H1, NW1, NH1))
m.merge(gf.Mesh("import", "structured",
                "GT='%s';ORG=[%f,%f];SIZES=[%f,%f];NSUBDIV=[%i,%i]"
                % (geotrans, 0, 0, W1, H2, NW1, NH2)))
m.merge(gf.Mesh("import", "structured",
                "GT='%s';ORG=[%f,%f];SIZES=[%f,%f];NSUBDIV=[%i,%i]"
                % (geotrans, W1, 0, W2, H2, NW2, NH2)))
m.optimize_structure()

l_rg = m.outer_faces_with_direction([-1,0], np.pi/180)
b_rg = m.outer_faces_with_direction([0,-1], np.pi/180)
r_rg = m.outer_faces_with_direction([1,0], np.pi/180)
t_rg = m.outer_faces_with_direction([0,1], np.pi/180)
in_rg = m.outer_faces_in_box([-1e-6,H1+H2-1e-6],[W1+1e-6,H1+H2+1e-6])
out_rg = m.outer_faces_in_box([W1+W2-1e-6,-1e-6],[W1+W2+1e-6,H2+1e-6])
m.set_region(IN_RG, in_rg)
m.set_region(OUT_RG, out_rg)
m.extend_region(INOUT_RG, in_rg)
m.extend_region(INOUT_RG, out_rg)

m.extend_region(WALL_RG, l_rg)
m.extend_region(WALL_RG, b_rg)
m.extend_region(WALL_RG, r_rg)
m.extend_region(WALL_RG, t_rg)
m.region_subtract(WALL_RG, INOUT_RG)

#MeshFem
mfv_ = gf.MeshFem(m, 2)
mfv_.set_classical_fem(2)
kept_dofs = np.setdiff1d(np.arange(mfv_.nb_basic_dof()),
                         mfv_.basic_dof_on_region(WALL_RG))
mfv = gf.MeshFem("partial", mfv_, kept_dofs)

mfp_ = gf.MeshFem(m, 1)
mfp_.set_classical_fem(1)
kept_dofs = np.setdiff1d(np.arange(mfp_.nb_basic_dof()),
                         mfp_.basic_dof_on_region(OUT_RG))
mfp = gf.MeshFem("partial", mfp_, kept_dofs)

mim = gf.MeshIm(m, 5) # 9 gauss points per quad

md = gf.Model("real")
md.add_fem_variable("v", mfv)
md.add_fem_data("v0", mfv)
md.add_fem_variable("p", mfp)
md.add_fem_data("p_in", mfp_)
md.add_initialized_data("f", [0,0])
md.add_initialized_data("dt", [dt])
md.add_initialized_data("nu", [nu])

md.add_Dirichlet_condition_with_multipliers(mim, "p", mfp, IN_RG, "p_in")
md.add_nonlinear_term\
(mim, "1/dt*(v-v0).Test_v + (Grad_v0*v0).Test_v + nu*Grad_v:Grad_Test_v - f.Test_v")
md.add_nonlinear_term\
(mim, "Grad_p.Grad_Test_p + 1/dt*Trace(Grad_v)*Test_p")

mmat_v = gf.asm_mass_matrix(mim, mfv)
#mmat_v = gf.asm_generic(mim, 2, "Test2_v.Test_v", -1, "v", 1, mfv, np.zeros(mfv.nbdof()))
IV = md.interval_of_variable("v")

t = 0
step = 0
while t < T+1e-8:
   print("Solving step at t=%f" % t)
   md.set_variable\
   ("p_in", mfp_.eval(p_in_str.format(t), globals(), locals()).flatten("F"))
   md.set_variable("v0", md.variable("v"))
   md.solve("noisy", "lsolver", "mumps", "max_res", 1e-8)
   vv = (gf.asm_generic(mim, 1, "(v-dt*Grad_p).Test_v", -1, md))[IV[0]:IV[0]+IV[1]]
   md.set_variable("v", gf.linsolve_mumps(mmat_v, vv))

   mfv.export_to_vtk("results_%i.vtk" % step,
                     mfv, md.variable("v"), "Velocity",
                     mfp, md.variable("p"), "Pressure")
   t += dt
   step += 1
