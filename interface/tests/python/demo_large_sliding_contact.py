#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Python GetFEM interface
#
# Copyright (C) 2012-2020 Yves Renard.
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

import numpy as np

import getfem as gf

gf.util_trace_level(1)

test_case = 1 # 0 = 2D punch on a rigid obstacle
              # 1 = 2D punch on a deformable obstacle (one slave, one master)
              # 2 = 2D with two different meshes
              # 3 = 2D with multi-body and only one mesh
              # 4 = 3D case (sphere / parallelepiped) (two meshes)

clambda1 = 1.   # Elasticity parameters
cmu1 = 1.
clambda2 = 1.   # Elasticity parameters
cmu2 = 1.
r = 0.1         # Augmentation parameter
alpha = 1.      # Alpha coefficient for "sliding velocity"
f_coeff = 0.    # Friction coefficient

test_tangent_matrix = False
nonlinear_elasticity = False
max_iter = 50

if test_case in [0,1]:
   vf = 0.
   vf_mult = 1.
   penalty_parameter = 0.
   dirichlet_translation = -0.5
   max_res = 1e-8
   release_dist = 1.5
   self_contact = False
   load_steps = 40
elif test_case == 3:
   vf = 0.01       # Vertical force
   vf_mult = 1.01
   penalty_parameter = 0.1
   release_dist = 0.05
   max_res = 1e-9
   self_contact = True
   load_steps = 250
elif test_case in [2,4]:
   vf = 0.01
   vf_mult = 1.5
   penalty_parameter = 0.01
   max_res = 1e-8
   if test_case == 2:
      release_dist = 0.1
   else:
      release_dist = 5.
   self_contact = True
   load_steps = 10000

if test_case == 0:
   #mesh1 = gf.Mesh("load", "../../../tests/meshes/punch2D_1.mesh")
   mesh1 = gf.Mesh("load", "../../../tests/meshes/punch2D_2.mesh")
elif test_case == 1:
   #mesh1 = gf.Mesh("load", "../../../tests/meshes/punch2D_1.mesh")
   mesh1 = gf.Mesh("load", "../../../tests/meshes/punch2D_2.mesh")
   mesh2 = gf.Mesh("import", "structured", "GT='GT_PK(2,1)';ORG=[-14,-5];SIZES=[28,5];NSUBDIV=[28,5]")
elif test_case == 2:
   mesh1 = gf.Mesh("load", "../../../tests/meshes/disc_with_a_hole.mesh")
   #mesh1 = gf.Mesh("import", "structured", "GT='GT_PK(2,1)';ORG=[-0.5,0.1];SIZES=[1,0.1];NSUBDIV=[20,2]")
   mesh2 = gf.Mesh("import", "structured", "GT='GT_PK(2,1)';ORG=[-0.5,0];SIZES=[1,0.1];NSUBDIV=[20,2]")
elif test_case == 3:
   mesh1 = gf.Mesh("load", "../../../tests/meshes/multi_body.mesh")
elif test_case == 4:
   mesh1 = gf.Mesh("load", "../../../tests/meshes/sphere_with_quadratic_tetra_400_elts.mesh")
   mesh2 = gf.Mesh("import", "structured", "GT='GT_PK(3,1)';ORG=[-15,-15,-4];SIZES=[30,30,4];NSUBDIV=[10,10,2]")
two_meshes = "mesh2" in locals()

N = mesh1.dim()

mfu1 = gf.MeshFem(mesh1, N)
mfu1.set_classical_fem(2)

pre_mflambda1 = gf.MeshFem(mesh1, N)
pre_mflambda1.set_classical_fem(1)

mfvm1 = gf.MeshFem(mesh1)
mfvm1.set_classical_discontinuous_fem(1)

CONTACT_BOUNDARY1 = 1
DIRICHLET_BOUNDARY1 = 3
border = mesh1.outer_faces()
if test_case >= 2:
   mesh1.set_region(CONTACT_BOUNDARY1, border)
else:
   normals = mesh1.normal_of_faces(border)
   contact_boundary = border[:,np.nonzero(normals[N-1] < -0.01)[0]]
   mesh1.set_region(CONTACT_BOUNDARY1, contact_boundary)
   P = mesh1.pts()  # get list of mesh points coordinates
   ctop = (P[N-1,:] > 39.999)  # find those on top of the object
   pidtop = np.compress(ctop, list(range(0, mesh1.nbpts())))
   ftop = mesh1.faces_from_pid(pidtop)
   mesh1.set_region(DIRICHLET_BOUNDARY1, ftop)

# dol1 = pre_mflambda1.basic_dof_on_region(CONTACT_BOUNDARY1)
# mflambda1 = gf.MeshFem("partial", pre_mflambda1, dol1)

mim1 = gf.MeshIm(mesh1, 4)
mim1_contact = gf.MeshIm(mesh1, 4)

if test_case not in [0,3]:
   mfu2 = gf.MeshFem(mesh2, N)
   mfu2.set_classical_fem(2)

   pre_mflambda2 = gf.MeshFem(mesh2, N)
   pre_mflambda2.set_classical_fem(1)

   mfvm2 = gf.MeshFem(mesh2)
   mfvm2.set_classical_discontinuous_fem(1)

   CONTACT_BOUNDARY2 = 2
   border = mesh2.outer_faces()
   if test_case != 1:
      mesh2.set_region(CONTACT_BOUNDARY2, border)
   else:
      normals = mesh2.normal_of_faces(border)
      contact_boundary = border[:,np.nonzero(normals[N-1] > 0.01)[0]]
      mesh2.set_region(CONTACT_BOUNDARY2, contact_boundary)
      dirichlet_boundary = border[:,np.nonzero(normals[N-1] < -0.01)[0]]
      DIRICHLET_BOUNDARY2 = 5
      mesh2.set_region(DIRICHLET_BOUNDARY2, dirichlet_boundary)

   mim2 = gf.MeshIm(mesh2, 4)
   mim2_contact = gf.MeshIm(mesh2, 4)

md = gf.Model("real")

F = np.zeros(N)
F[N-1] = -vf

w1_str = ""
w2_str = ""

md.add_fem_variable("u1", mfu1)
if f_coeff > 1e-10:
   md.add_fem_data("w1", mfu1)
   w1_str = "w1"
md.add_filtered_fem_variable("lambda1", pre_mflambda1, CONTACT_BOUNDARY1)

if nonlinear_elasticity:
   lawname = "Ciarlet Geymonat"
   params1 = [clambda1, cmu1, cmu1/2-clambda1/8]
   md.add_initialized_data("params1", params1)
   md.add_nonlinear_elasticity_brick(mim1, "u1", lawname, "params1")
else:
   md.add_initialized_data("clambda1", clambda1)
   md.add_initialized_data("cmu1", cmu1)
   md.add_isotropic_linearized_elasticity_brick(mim1, "u1", "clambda1", "cmu1")

if test_case == 2:
#   md.add_initialized_data("cpoints1", [0 0.5 0 1.5 0 0.5 0 1.5])
#   md.add_initialized_data("cunitv1", [1 0 1 0 0 1 0 1])
#   md.add_initialized_data("cdata", [0 0 -0.01 -0.01])
#   md.add_pointwise_constraints_with_multipliers("u1", "cpoints1", "cunitv1", "cdata")
   md.add_initialized_data("cpoints1", [0,0.5,0,1.5])
   md.add_initialized_data("cunitv1", [1,0,1,0])
   md.add_initialized_data("cdata", [0,0])
   md.add_pointwise_constraints_with_multipliers("u1", "cpoints1", "cunitv1", "cdata")

md.add_initialized_data("penalty_param1", [penalty_parameter])
md.add_mass_brick(mim1, "u1", "penalty_param1")
md.add_initialized_data("data1", F)
md.add_source_term_brick(mim1, "u1", "data1")

if test_case not in [0,3]:
   md.add_fem_variable("u2", mfu2)
   if f_coeff > 1e-10:
      md.add_fem_data("w2", mfu2)
      w2_str = "w2"
   if self_contact:
      md.add_filtered_fem_variable("lambda2", pre_mflambda2, CONTACT_BOUNDARY2)

   if nonlinear_elasticity:
      lawname = "Ciarlet Geymonat"
      params2 = [clambda2, cmu2, cmu2/2-clambda2/8]
      md.add_initialized_data("params2", params2)
      md.add_nonlinear_elasticity_brick(mim2, "u2", lawname, "params2")
   else:
      md.add_initialized_data("clambda2", clambda2)
      md.add_initialized_data("cmu2", cmu2)

      md.add_isotropic_linearized_elasticity_brick(mim2, "u2", "clambda2", "cmu2")

   if test_case == 2:
      md.add_initialized_data("cpoints2", [0,0])
      md.add_initialized_data("cunitv2", [1,0])
      md.add_pointwise_constraints_with_multipliers("u2", "cpoints2", "cunitv2")

   md.add_initialized_data("penalty_param2", [penalty_parameter])
   md.add_mass_brick(mim2, "u2", "penalty_param2")
   md.add_initialized_data("data2", F)
   md.add_source_term_brick(mim2, "u2", "data2")

   if test_case == 1:
      Ddata = np.zeros(N)
      md.add_initialized_data("Ddata2", Ddata)
      md.add_Dirichlet_condition_with_multipliers(mim2, "u2", 1, DIRICHLET_BOUNDARY2, "Ddata2")

if test_case <= 1:
   Ddata = np.zeros(N)
   Ddata[N-1] = dirichlet_translation
   md.add_initialized_data("Ddata1", Ddata)
   md.add_Dirichlet_condition_with_multipliers(mim1, "u1", 1, DIRICHLET_BOUNDARY1, "Ddata1")

md.add_initialized_data("r", r)
md.add_initialized_data("alpha", alpha)
md.add_initialized_data("f", f_coeff)

direct_generic_assembly = False
if direct_generic_assembly:  # Direct use of high-level generic assembly
  # TODO: account for w1, w2 when f_coeff > 0
  md.add_raytracing_transformation("contact_trans", release_dist)
  if two_meshes: # The definition of a variable group is not mandatory. Just for test.
    md.define_variable_group("u", "u1", "u2")
  else:
    md.define_variable_group("u", "u1")

  if self_contact:
    md.add_master_contact_boundary_to_raytracing_transformation("contact_trans", mesh1, "u", CONTACT_BOUNDARY1)
  else:
    md.add_slave_contact_boundary_to_raytracing_transformation("contact_trans", mesh1, "u", CONTACT_BOUNDARY1)

  if test_case == 0:
    md.add_rigid_obstacle_to_raytracing_transformation("contact_trans", "80-sqrt(sqr(x)+sqr(y-80))", N)
  elif test_case == 1:
    md.add_master_contact_boundary_to_raytracing_transformation("contact_trans", mesh2, "u", CONTACT_BOUNDARY2)
  elif test_case == 2:
    md.add_master_contact_boundary_to_raytracing_transformation("contact_trans", mesh2, "u", CONTACT_BOUNDARY2)
    md.add_rigid_obstacle_to_raytracing_transformation("contact_trans", "y+1", N)
  elif test_case == 3:
    md.add_rigid_obstacle_to_raytracing_transformation("contact_trans", "2-sqrt(sqr(x)+sqr(y-1))", N)
  elif test_case == 4:
    md.add_master_contact_boundary_to_raytracing_transformation("contact_trans", mesh2, "u", CONTACT_BOUNDARY2)
    md.add_rigid_obstacle_to_raytracing_transformation("contact_trans", "z+5", N)

  md.add_nonlinear_term(mim1_contact, "-lambda1.Test_u1", CONTACT_BOUNDARY1) 
  md.add_nonlinear_term(mim1_contact, "Interpolate_filter(contact_trans, lambda1.Interpolate(Test_u,contact_trans), 1)", CONTACT_BOUNDARY1) 
  md.add_nonlinear_term(mim1_contact, "-(1/r)*lambda1.Test_lambda1", CONTACT_BOUNDARY1)
  md.add_nonlinear_term(mim1_contact, "Interpolate_filter(contact_trans, (1/r)*Coulomb_friction_coupled_projection(lambda1, Transformed_unit_vector(Grad_u1, Normal), u1, (Interpolate(X,contact_trans)-X-u1).Transformed_unit_vector(Grad_u1, Normal), f, r).Test_lambda1, 2)", CONTACT_BOUNDARY1)
  md.add_nonlinear_term(mim1_contact, "Interpolate_filter(contact_trans, (1/r)*Coulomb_friction_coupled_projection(lambda1, Transformed_unit_vector(Grad_u1, Normal), u1-Interpolate(u,contact_trans), (Interpolate(X,contact_trans)+Interpolate(u,contact_trans)-X-u1).Transformed_unit_vector(Grad_u1, Normal), f, r).Test_lambda1, 1)", CONTACT_BOUNDARY1)
  
  if two_meshes and self_contact:
    md.add_nonlinear_term(mim2_contact, "-lambda2.Test_u2", CONTACT_BOUNDARY2) 
    md.add_nonlinear_term(mim2_contact, "Interpolate_filter(contact_trans, lambda2.Interpolate(Test_u,contact_trans), 1)", CONTACT_BOUNDARY2) 
    md.add_nonlinear_term(mim2_contact, "-(1/r)*lambda2.Test_lambda2", CONTACT_BOUNDARY2)
    md.add_nonlinear_term(mim2_contact, "Interpolate_filter(contact_trans, (1/r)*Coulomb_friction_coupled_projection(lambda2, Transformed_unit_vector(Grad_u2, Normal), u2, (Interpolate(X,contact_trans)-X-u2).Transformed_unit_vector(Grad_u2, Normal), f, r).Test_lambda2, 2)", CONTACT_BOUNDARY2)
    md.add_nonlinear_term(mim2_contact, "Interpolate_filter(contact_trans, (1/r)*Coulomb_friction_coupled_projection(lambda2, Transformed_unit_vector(Grad_u2, Normal), u2-Interpolate(u,contact_trans), (Interpolate(X,contact_trans)+Interpolate(u,contact_trans)-X-u2).Transformed_unit_vector(Grad_u2, Normal), f, r).Test_lambda2, 1)", CONTACT_BOUNDARY2)  

  u_group = "u"
  contact_trans = "contact_trans"
else: # Use of the new contact brick which uses the high-level generic assembly

  ind = md.add_integral_large_sliding_contact_brick_raytracing("r", release_dist, "f", "alpha", 0)

  if self_contact:
    md.add_master_slave_contact_boundary_to_large_sliding_contact_brick(ind,  mim1_contact, CONTACT_BOUNDARY1, "u1", "lambda1", w1_str)
  else:
    md.add_slave_contact_boundary_to_large_sliding_contact_brick(ind,  mim1_contact, CONTACT_BOUNDARY1, "u1", "lambda1", w1_str)

  if test_case == 0:
    md.add_rigid_obstacle_to_large_sliding_contact_brick(ind, "80-sqrt(sqr(x)+sqr(y-80))", N)
  elif test_case == 1:
    md.add_master_contact_boundary_to_large_sliding_contact_brick(ind, mim2_contact, CONTACT_BOUNDARY2, "u2", w2_str)
  elif test_case == 2:
    md.add_master_slave_contact_boundary_to_large_sliding_contact_brick(ind, mim2_contact, CONTACT_BOUNDARY2, "u2", "lambda2", w2_str)
    md.add_rigid_obstacle_to_large_sliding_contact_brick(ind, "y+1", N)
  elif test_case == 3:
    md.add_rigid_obstacle_to_large_sliding_contact_brick(ind, "2-sqrt(sqr(x)+sqr(y-1))", N)
  elif test_case == 4:
    md.add_master_slave_contact_boundary_to_large_sliding_contact_brick(ind, mim2_contact, CONTACT_BOUNDARY2, "u2", "lambda2", w2_str)
    md.add_rigid_obstacle_to_large_sliding_contact_brick(ind, "z+5", N)

  u_group = md.displacement_group_name_of_large_sliding_contact_brick(ind)
  contact_trans = md.transformation_name_of_large_sliding_contact_brick(ind)


for nit in range(load_steps):

   if test_tangent_matrix:
      errmax = md.test_tangent_matrix(1E-8, 20, 0.0001)
      #errmax = md.test_tangent_matrix_term("lambda1", "u1", 1E-8, 20, 0.0001)
      print("errmax = %g" % errmax)
      if errmax > 1e-3:
         print("bad tangent matrix")

   if w1_str:
      md.set_variable(w1_str, md.variable("u1"))
   if w2_str:
      md.set_variable(w2_str, md.variable("u2"))

   print("SOLVING LOAD STEP %i" % nit)
   md.solve("noisy", "max_iter", max_iter, "max_res", max_res) # , "lsearch", "simplest")

   U1 = md.variable("u1")
   if nonlinear_elasticity:
      VM1 = md.compute_Von_Mises_or_Tresca("u1", lawname, "params1", mfvm1)
   else:
      VM1 = md.compute_isotropic_linearized_Von_Mises_or_Tresca("u1", "clambda1", "cmu1", mfvm1)
   mfvm1.export_to_vtk("lsc_1_%i.vtk" % nit, mfvm1,  VM1,
                       "Von Mises Stresses 1", mfu1, U1, "Displacements 1")

   lambda1 = md.variable("lambda1")
   mf_lambda1 = md.mesh_fem_of_variable("lambda1")
   sl = gf.Slice(("boundary",), mf_lambda1, CONTACT_BOUNDARY1)
   sl.export_to_vtk("lsc_1_boundary_%i.vtk" % nit,
                    mfu1, U1, "BDisplacements 1",
                    mf_lambda1, lambda1, "BMultiplier 1")

   if test_case not in [0,3]:
      U2 = md.variable("u2")
      if nonlinear_elasticity:
         VM2 = md.compute_Von_Mises_or_Tresca("u2", lawname, "params2", mfvm2)
      else:
         VM2 = md.compute_isotropic_linearized_Von_Mises_or_Tresca("u2", "clambda2", "cmu2", mfvm2)
      mfvm2.export_to_vtk("lsc_2_%i.vtk" % nit, mfvm2,  VM2,
                          "Von Mises Stresses 2", mfu2, U2, "Displacements 2")

      sl = gf.Slice(("boundary",), mfu2, CONTACT_BOUNDARY2)
      sl.export_to_vtk("lsc_2_boundary_%i.vtk" % nit,
                       mfu2, U2, "BDisplacements 2")

   vf *= vf_mult
   F[N-1] = -vf
   md.set_variable("data1", F)
   if test_case not in [0,3]:
      md.set_variable("data2", F)

   if test_case <= 1:
      Ddata[N-1] -= 1.
      md.set_variable("Ddata1", Ddata)
