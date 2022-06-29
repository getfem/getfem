#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Python GetFEM interface
#
# Copyright (C) 2020-2020 Konstantinos Poulios.
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
"""  This example simulates necking during uniaxial tension of a plate
     under plane strain, with a hyperelastoplastic constitutive law with
     linear isotropic hardening. Only the cross section of the plate is
     modeled, using quadrilateral 2D elements.

     This is a reference implementation of finite strain plasticity with
     linear hardening in 2D (plane strain) using the generic weak form
     language of GetFEM.
"""
import getfem as gf
import numpy as np
import os, sys, subprocess
import time
import shutil
gf.util_trace_level(1)
gf.util_warning_level(1)

# Input data
L = 2*26.667        # block length
H = 2*6.413         # block height
dH = 0.018*H        # height reduction at the center of the block

N_L = 20            # number of elements in block length direction
N_H = 10            # number of elements in block height direction

E = 210e3           # Young's modulus
nu = 0.3            # Poisson's ratio
pl_sigma_0 = 5e2    # Initial yield stress
pl_H = 21e1         # Plastic modulus (0.1% of E)

disp = 8.           # maximum displacement
steps_t = 200       # number of load steps for increasing the load

#------------------------------------
geotrans = "GT_QK(2,2)"  # geometric transformation

disp_fem_order = 2  # displacements finite element order
mult_fem_order = 2  # dirichlet multipliers finite element order

#integration_degree = 3 # 4 gauss points per quad
integration_degree = 5 # 9 gauss points per quad

export_results = False;

#------------------------------------
resultspath = "./results"
if (export_results and not os.path.exists(resultspath)):
   os.makedirs(resultspath)

tee = subprocess.Popen(["tee", "%s/tension_plane_strain.log" % resultspath],
                       stdin=subprocess.PIPE)
sys.stdout.flush()
os.dup2(tee.stdin.fileno(), sys.stdout.fileno())
sys.stderr.flush()
os.dup2(tee.stdin.fileno(), sys.stderr.fileno())

# auxiliary constants
XM_RG = 1
XP_RG = 2
YM_RG = 3

# Mesh
xmin = 0.
ymin = 0.
dx = L/2
dy = H/2
mesh = gf.Mesh("import", "structured",
               "GT='%s';ORG=[%f,%f];SIZES=[%f,%f];NSUBDIV=[%i,%i]"
               % (geotrans, xmin, ymin, dx, dy, N_L, N_H))
N = mesh.dim()

outer_faces = mesh.outer_faces()
outer_normals = mesh.normal_of_faces(outer_faces)
mesh.set_region(XM_RG,
                outer_faces[:,np.nonzero(outer_normals[0] < -0.95)[0]])
mesh.set_region(XP_RG,
                outer_faces[:,np.nonzero(outer_normals[0] >  0.95)[0]])
mesh.set_region(YM_RG,
                outer_faces[:,np.nonzero(outer_normals[1] < -0.95)[0]])

if dH > 0:
   pts = mesh.pts()
   dy = 0.2   # fine element size ratio w.r.t. uniform mesh size
   y0 = 0.04  # ratio of the finer meshed area w.r.t. the total length
   x0 = y0/dy
   b2 = (1-dy)/(x0-1)**2
   b1 = dy - 2*b2*x0
   b0 = 1 - b1 -b2
   for i in range(pts.shape[1]):
      x = pts[0,i]
      y = pts[1,i]
      t = 2*abs(x)/L
      if t < x0:
        x *= dy;
      else:
        x = (b0 + b1*t + b2*t**2) * np.sign(x)*L/2
      pts[0,i] = x
      pts[1,i] -= (y*dH)/(2*H) * (1 + np.cos(2.*np.pi*x/L))
   mesh.set_pts(pts)
   
if (export_results):
   mesh.export_to_vtu("%s/mesh.vtu" % resultspath)

# FEM
mfN = gf.MeshFem(mesh, N)
mfN.set_classical_fem(disp_fem_order)
#if disp_fem_order == 2:
#   mfN.set_fem(gf.Fem("FEM_Q2_INCOMPLETE(2)"))
#else:
#   mfN.set_classical_fem(disp_fem_order)
keptdofs = np.arange(mfN.nbdof())
keptdofs = np.setdiff1d(keptdofs, mfN.basic_dof_on_region(XM_RG)[0::N])
keptdofs = np.setdiff1d(keptdofs, mfN.basic_dof_on_region(YM_RG)[1::N])
mfu = gf.MeshFem("partial", mfN, keptdofs)

mfmult = gf.MeshFem(mesh, 1)
mfmult.set_classical_fem(mult_fem_order)

mfout1 = gf.MeshFem(mesh)
mfout1.set_classical_discontinuous_fem(disp_fem_order-1)
mfout2 = gf.MeshFem(mesh)
mfout2.set_classical_discontinuous_fem(disp_fem_order)

mim = gf.MeshIm(mesh, integration_degree)

mimd1 = gf.MeshImData(mim)
mimd4 = gf.MeshImData(mim, -1, 4)

# Model
md = gf.Model("real")

md.add_fem_variable("u", mfu)

# Vertical displacement
md.add_initialized_data("disp", [0.])

md.add_initialized_data("K", E/(3.*(1.-2.*nu))) # Bulk modulus
md.add_initialized_data("mu", E/(2*(1+nu)))     # Shear modulus
md.add_macro("F", "Id(2)+Grad_u")
md.add_macro("F3d", "Id(3)+[1,0;0,1;0,0]*Grad_u*[1,0,0;0,1,0]")
md.add_macro("J", "Det(F)")
md.add_macro("tauH", "K*log(J)")
md.add_im_data("gamma0", mimd1)                        # accumulated plastic strain at previous time step
md.add_im_data("invCp0vec", mimd4)                     # Components 11, 22, 33 and 12 of the plastic part of
md.set_variable("invCp0vec",                           # the inverse right Cauchy Green tensor at the previous
                np.tile([1,1,1,0], mimd4.nbpts()))     # step. Symmetric components are omitted.
md.add_macro("invCp0", "[[[1,0,0],[0,0,0],[0,0,0]],"+\
                       " [[0,0,0],[0,1,0],[0,0,0]],"+\
                       " [[0,0,0],[0,0,0],[0,0,1]],"+\
                       " [[0,1,0],[1,0,0],[0,0,0]]].invCp0vec") #Vec4ToMat3x3
md.add_macro("devlogbetr", "Deviator(Logm(F3d*invCp0*F3d'))")
md.add_macro("Y0","{A}+{B}*gamma0".format(A=np.sqrt(2./3.)*pl_sigma_0, B=2./3.*pl_H))
md.add_macro("ksi",
             "(1-1/max(1,mu*pow(J,-5/3)*Norm(devlogbetr)/Y0))/(2+{B}/(mu*pow(J,-5/3)))"
             .format(B=2./3.*pl_H))
md.add_macro("gamma", "gamma0+ksi*Norm(devlogbetr)")
md.add_macro("devlogbe", "(1-2*ksi)*devlogbetr")
md.add_macro("tauD2d", "mu*pow(J,-2/3)*[1,0,0;0,1,0]*devlogbe*[1,0;0,1;0,0]")

md.add_nonlinear_generic_assembly_brick\
   (mim, "((tauH*Id(2)+tauD2d)*Inv(F')):Grad_Test_u")

md.add_macro("sigmaD", "(mu*pow(J,-5/3)*devlogbe)")
md.add_macro("sigma", "tauH/J*Id(3)+mu*pow(J,-5/3)*devlogbe")
md.add_macro("invCp", "(Inv(F3d)*Expm(-ksi*devlogbetr)*(F3d))*invCp0"\
                      "*(Inv(F3d)*Expm(-ksi*devlogbetr)*(F3d))'")

# Dirichlet condition
md.add_filtered_fem_variable("dirmult", mfmult, XP_RG)
md.add_nonlinear_generic_assembly_brick(mim, "(disp-u(1))*dirmult", XP_RG)

print("Model dofs: %i" % md.nbdof())
print("Displacement fem dofs: %i" % mfu.nbdof())
print("Dirichlet mult dofs: %i" % md.mesh_fem_of_variable("dirmult").nbdof())

if (export_results):
   shutil.copyfile(os.path.abspath(sys.argv[0]),resultspath+"/"+sys.argv[0])

starttime_overall = time.process_time()

if (export_results):
   f1 = open("%s/tension_plane_strain.dat" % resultspath, "w")

try:
   for step in range(steps_t+1):
      md.set_variable("disp", disp*step/float(steps_t))
      print('STEP %i: Solving with disp = %g' % (step, md.variable("disp")))

      starttime = time.process_time()
      md.solve('noisy', 'max_iter', 25, 'max_res', 1e-10,
               'lsearch', 'simplest', 'alpha max ratio', 100, 'alpha min', 1, 'alpha mult', 0.6,
               'alpha threshold res', 1e9)
      print('STEP %i COMPLETED IN %f SEC' % (step, time.process_time()-starttime))

      F = gf.asm_generic(mim, 0, "dirmult", XP_RG, md)
      print("Displacement %g, total force %g" % (md.variable("disp"), F))
      A = gf.asm_generic(mim, 0, "Norm(J*Inv(F')*[1;0])", XP_RG, md)
      V = gf.asm_generic(mim, 0, "1", -1, md)
      sigma11 = gf.asm_generic(mim, 0, "sigma(1,1)", -1, md)/V
      gamma = gf.asm_generic(mim, 0, "gamma", -1, md)/V
      if (export_results):
         f1.write("%.10g %.10g %.10g %.10g %10g %10g\n"
                  % (md.variable("disp"), F, A, F/A, sigma11, gamma))
         f1.flush()

         output = (mfout1, md.local_projection(mim, "sqrt(1.5)*Norm(sigmaD)", mfout1), "Von Mises Stress",
                mfout1, md.local_projection(mim, "J", mfout1), "J",
                mfout1, md.local_projection(mim, "sigma(1,1)", mfout1), "Cauchy stress 11",
                mfout1, md.local_projection(mim, "sigma(2,2)", mfout1), "Cauchy stress 22",
                mfout1, md.local_projection(mim, "sigma(1,2)", mfout1), "Cauchy stress 12",
                mfout1, md.local_projection(mim, "sigma(3,3)", mfout1), "Cauchy stress 33",
                mfu, md.variable("u"), "Displacements",
                mfout2, md.interpolation("dirmult", mfout2, XP_RG), "Nominal reaction traction",
                mfout2, md.local_projection(mim, "gamma", mfout2), "plastic strain")
         
         mfout2.export_to_vtu("%s/tension_plane_strain_%i.vtu" % (resultspath, step), *output)

      md.set_variable("gamma0", md.interpolation("gamma", mimd1, -1))
      md.set_variable("invCp0vec",
                      md.interpolation("[[[1,0,0,0]  ,[0,0,0,0.5],[0,0,0,0]],"+\
                                       " [[0,0,0,0.5],[0,1,0,0]  ,[0,0,0,0]],"+\
                                       " [[0,0,0,0]  ,[0,0,0,0]  ,[0,0,1,0]]]:invCp", mimd4, -1))
finally:
   if (export_results) :
      f1.close()

print('OVERALL SOLUTION TIME IS %f SEC' % (time.process_time()-starttime_overall))

