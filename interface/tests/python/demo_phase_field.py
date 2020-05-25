#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Python GetFEM interface
#
# Copyright (C) 2018-2020 Konstantinos Poulios.
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

import os

import numpy as np

import getfem as gf

np.set_printoptions(threshold=100000)
gf.util_trace_level(1)

# Input data
NX = 60         # number of elements in horizontal direction
NY = 60         # number of elements in verical direction

LX = 1.         # [mm] Length
LY = 1.         # [mm] Height

E = 210e3       # [N/mm^2]
nu = 0.3

Gc0 = 2.7e0     # [N/mm]

ll = 0.015      # [mm] length scale

strain_rate = 2e-4
init_time_step = 1.

steps = 100

#------------------------------------
geotrans = "GT_QK(2,2)"  # geometric transformation
disp_fem_order = 2       # displacements finite element order
phi_fem_order = 2        # phase field finite element order

integration_degree = 5   # 9 gauss points per quad

#------------------------------------

NX = int(2*((NX+1)/2))
NY = int(2*((NY+1)/2))

resultspath = "./demo_phase_field_results"
if not os.path.exists(resultspath):
   os.makedirs(resultspath)

# auxiliary constants
B_BOUNDARY = 4
T_BOUNDARY = 6
TB_BOUNDARY = 7

NX_seed1 = np.linspace(-1., 0., int(NX/3+1))
NX_seed2 = np.linspace(0., 1., int((NX-NX/3)+1))[1:]
NY_seed = np.linspace(-1., 1., NY+1)
X_seed1 = LX/2*(0.2*NX_seed1+0.8*np.sign(NX_seed1)*np.power(np.abs(NX_seed1),1.5))
X_seed2 = LX/2*(0.6*NX_seed2+0.4*np.sign(NX_seed2)*np.power(np.abs(NX_seed2),1.5))
X_seed = np.concatenate((X_seed1,X_seed2))
Y_seed = LY/2*(0.2*NY_seed+0.8*np.sign(NY_seed)*np.power(np.abs(NY_seed),1.5))
m1 = gf.Mesh("cartesian", X_seed, Y_seed[0:int(NY/2+1)])
m2 = gf.Mesh("cartesian", X_seed, Y_seed[int(NY/2):])
gap = 1e-5
pts = m1.pts()
for i in range(pts.shape[1]):
  if pts[0,i] < -1e-10:
    pts[1,i] -=  gap/2*(1+pts[1,i]/(LY/2))
m1.set_pts(pts)
pts = m2.pts()
for i in range(pts.shape[1]):
  if pts[0,i] < -1e-10:
    pts[1,i] +=  gap/2*(1-pts[1,i]/(LY/2))
m2.set_pts(pts)
mesh = m1
mesh.merge(m2)
N = mesh.dim()

bottom_faces = mesh.outer_faces_in_box([-LX/2-1e-5,-LY/2-1e-5],[LX/2+1e-5,-LY/2+1e-5])
top_faces = mesh.outer_faces_in_box([-LX/2-1e-5,LY/2-1e-5],[LX/2+1e-5,LY/2+1e-5])

mesh.set_region(B_BOUNDARY, bottom_faces)
mesh.set_region(T_BOUNDARY, top_faces)
mesh.region_merge(TB_BOUNDARY, T_BOUNDARY)
mesh.region_merge(TB_BOUNDARY, B_BOUNDARY)

mesh.export_to_vtk("%s/mesh.vtk" % resultspath)

# FEM
mfu = gf.MeshFem(mesh, N)
mfu.set_classical_fem(disp_fem_order)

mfdir = mfu

mfphi = gf.MeshFem(mesh, 1)
mfphi.set_classical_fem(phi_fem_order)

mfout = gf.MeshFem(mesh)
mfout.set_classical_discontinuous_fem(2)

# Integration method
mim = gf.MeshIm(mesh, integration_degree)
mimd1 = gf.MeshImData(mim, -1)

# Model
md = gf.Model("real")

md.add_fem_variable("u", mfu)          # displacements
md.add_fem_variable("phi", mfphi)      # phase field

md.add_fem_data("u_stored", mfu)
md.add_fem_data("phi_stored", mfphi)

md.add_im_data("psi0_max", mimd1)

md.add_initialized_data("kappa", E/(3.*(1.-2.*nu)))
md.add_initialized_data("mu", E/(2*(1+nu)))

md.add_initialized_data("Gc0", Gc0)
md.add_initialized_data("l", ll)

md.add_macro("damage", "max(1e-7,sqr(1-phi))")
md.add_macro("psi0", "(0.5*kappa*sqr(Trace(Grad_u))+mu*Norm_sqr((Sym(Grad_u)-Trace(Grad_u)/3*Id(2))))")
md.add_macro("Gc", "Gc0")

_sigma_ = "damage*(kappa*Trace(Grad_u)*Id(2)+2*mu*(Sym(Grad_u)-Trace(Grad_u)/3*Id(2)))"
md.add_nonlinear_term(mim, _sigma_+":Grad_Test_u")
md.add_nonlinear_term(mim, "(-2*(1-phi)*max(psi0_max,psi0)*Test_phi+Gc*(phi/l*Test_phi+l*Grad_phi.Grad_Test_phi))")

md.add_initialized_data("MM", 0)
md.add_nonlinear_term(mim, "MM*((u-u_stored).Test_u+1e2*(phi-phi_stored)*Test_phi)")

# Load
md.add_fem_data("dirichlet_data", mfu);
ibdir = md.add_Dirichlet_condition_with_multipliers(mim, "u", mfdir, TB_BOUNDARY, "dirichlet_data")
dirmultname = md.mult_varname_Dirichlet(ibdir)

mass_mat = gf.asm_mass_matrix(mim, mfout)

print("Displacement dofs: %i" % mfu.nbdof())
print("Total model dofs: %i" % md.nbdof())

time_step = init_time_step
eps = 0.
pseudodynamic = False
MM = 0.
with open("%s/demo_phase_field_forces.dat" % resultspath, "w") as f:
  for step in range(steps):
    eps_old = eps
    while True:
      if step > 0:
        eps = eps_old + strain_rate*time_step

      X = md.from_variables()
      print("Step %i with eps=%e" % (step, eps))
      md.set_variable("dirichlet_data", md.interpolation("{eps}*[0;X(2)]".format(eps=eps), mfu, -1))
      if pseudodynamic:
        print("With damping %e" % MM)
      md.set_variable("MM",[MM])

      iters = 20
      try:
        nit = \
        md.solve("noisy", "lsolver", "mumps", "max_iter", iters, "max_res", 1e-7, #)[0]
                 "lsearch", "simplest", "alpha max ratio", 1.5, "alpha min", 0.4, "alpha mult", 0.6,
                  "alpha threshold res", 2)[0]
      except (KeyboardInterrupt, SystemExit):
        raise
      except:
        nit = iters

      if nit >= iters:
        if step == 0:
          break
        md.to_variables(X)
        eps = eps_old
        if time_step > init_time_step/50.:
          time_step *= 0.5
        else:
          if pseudodynamic:
            MM *= 100.
          else:
            pseudodynamic = True
            MM = 10000.
      else:
        md.set_variable("u_stored", md.variable("u"))
        md.set_variable("phi_stored", md.variable("phi"))
        if pseudodynamic:
          if MM < 1e-4:
            MM = 0.
            pseudodynamic = False
          elif nit <= 6:
            MM /= 2.
        else:
          TMP = md.interpolation("max(psi0_max,psi0)", mimd1, -1)
          md.set_variable("psi0_max", TMP)
          break

    out = (mfu, md.variable("u"), "Displacements",
           mfphi, md.variable("phi"), "phi")
    for i,j in [[1,1],[2,2],[1,2]]:
      sigma = gf.asm_generic(mim, 1, "{sigma}({i},{j})*Test_t".format(sigma=_sigma_, i=i, j=j),
                             -1, md, "t", True, mfout, np.zeros(mfout.nbdof()))\
              [md.nbdof():]
      sigma = np.transpose(gf.linsolve_mumps(mass_mat, sigma))
      out += (mfout, sigma, "Cauchy Stress {i}{j}".format(i=i, j=j))

    mfout.export_to_vtk("%s/demo_phase_field_%i.vtk" % (resultspath, step), *out)

    DIRMULT = -md.variable(dirmultname)
    DIRMULT = np.reshape(DIRMULT,[1,DIRMULT.size])
    dfT = gf.asm_boundary_source(T_BOUNDARY, mim, mfu, md.mesh_fem_of_variable(dirmultname), DIRMULT)
    f.write(("step=%i eps=%e fR=(%e,%e)\n") %
            (step, eps, dfT[0::N].sum(), dfT[1::N].sum()))
    f.flush()
