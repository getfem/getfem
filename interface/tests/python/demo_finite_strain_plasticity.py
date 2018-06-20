#!/usr/bin/env python
# -*- coding: UTF8 -*-
#
# Copyright (C) 2014-2018 Konstantinos Poulios.
#
# This file is a part of GetFEM++
#
# GetFEM++  is  free software;  you  can  redistribute  it  and/or modify it
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

import getfem as gf
import numpy as np
import sys
import os
import time
import shutil

np.set_printoptions(threshold=np.nan)
gf.util_trace_level(1)

# Input data
UL = 1.             # unit length
UF = 1.             # unit force

L = 2*26.667*UL     # block length
H = 2*6.413*UL      # block height
dH = 0.018*H        # height reduction at the center of the block

N_L = 20        # number of elements in block length direction
N_H = 10        # number of elements in block height direction

E = 206.9e3           # Elastic Young's modulus
nu = 0.29             # Poisson's ratio
sigma_0 = 0.45e3      # Initial yield stress
sigma_inf = 0.715e3
h1 = 16.93
h2 = -0.012924e3      # Hardening/softening tangent modulus

r_aug = 1e4       # augmentation parameter

steps = 200
#dirichlet_str = '[0,0]'
#dirichlet_str = '[{0}*{1}*y,{0}*{1}*x]'.format('{0}',0.4/steps)  #  pure shearing with 0.4 the shearing angle
#dirichlet_str = '[{0}*{1}*x*y,0]'.format('{0}',0.4/L/steps)  # symmetric xy shearing with 0.2 the shearing angle
#dirichlet_str = '[0,{0}*{1}*x*y]'.format('{0}',0.4/L/steps)  # symmetric yx shearing with 0.2 the shearing angle
#dirichlet_str = '[{0}*{1}*y,0]'.format('{0}',0.4/steps)  #  simple xy shearing with 0.4 the shearing angle
#dirichlet_str = '[0,{0}*{1}*x]'.format('{0}',0.4/steps)  #  simple yx shearing with 0.4 the shearing angle
dirichlet_str = '[{0}*{1}*x,0]'.format('{0}',10./L/steps)      # x stretching with 0.2 the stretching ratio
#dirichlet_str = '[0,{0}*{1}*y]'.format('{0}',0.1/steps)      # y stretching with 0.8 the stretching ratio
#dirichlet_str = '[np.cos({0}*{1})*x-np.sin({0}*{1})*y-x,'\
#                'np.sin({0}*{1})*x+np.cos({0}*{1})*y-y]'.format('{0}',np.pi/4/steps)   # rotation by pi/4
#dirichlet_str = '[np.cos({0}*{1})*(1+{0}*{2})*x-np.sin({0}*{1})*y-x,'\
#                'np.sin({0}*{1})*(1+{0}*{2})*x+np.cos({0}*{1})*y-y]'.format('{0}',np.pi/4/steps,1.2/steps)   # rotation by pi/4 and stretching (0.2)
#dirichlet_str = '[np.cos({0}*{1})*(x+{0}*{2}*y)-np.sin({0}*{1})*y-x,'\
#               'np.sin({0}*{1})*(x+{0}*{2}*y)+np.cos({0}*{1})*y-y]'.format('{0}',np.pi/4/steps,0.4/steps)   # rotation by pi/4 and shearing xy
#dirichlet_str = '[np.cos({0}*{1})*x-np.sin({0}*{1})*(y+{0}*{2}*x)-x,'\
#                'np.sin({0}*{1})*x+np.cos({0}*{1})*(y+{0}*{2}*x)-y]'.format('{0}',np.pi/4/steps,0.4/steps)   # rotation by pi/4 and shearing yx
#dirichlet_str = '[np.cos({0}*{1})*(x+{0}*{2}*y)-np.sin({0}*{1})*(y+{0}*{2}*x)-x,'\
#                'np.sin({0}*{1})*(x+{0}*{2}*y)+np.cos({0}*{1})*(y+{0}*{2}*x)-y]'.format('{0}',np.pi/4/steps,0.4/steps)   # rotation by pi/4 and pure shearing

#dirichlet_str = '[{0}*{1}*x,-{0}*{1}/(1+{0}*{1})*y]'.format('{0}',0.1/steps)      # isochoric stretching

#------------------------------------
geotrans = 'GT_QK(2,2)'  # geometric transformation
#geotrans = 'GT_Q2_INCOMPLETE(2)'

disp_fem_order = 2   # displacements finite element order
press_fem_order = 1  # pressure finite element order
mult_fem_order = 2   # plastic multiplier finite element order

#integration_degree = 3 # 4 gauss points per quad
integration_degree = 5 # 9 gauss points per quad
#integration_degree = 7 # 16 gauss points per quad
#integration_degree = 9 # 25 gauss points per quad

integration_degree_press = 5

symmetric = True
rigid_grip = False

#------------------------------------

resultspath = './results_finite_strain_plasticity_deSouzaNeto'
if not os.path.exists(resultspath):
    os.makedirs(resultspath)


# auxiliary parameters
K = E/(3.*(1.-2.*nu))  # Bulk modulus
mu = E/(2.*(1.+nu))    # Shear modulus

# auxiliary constants
L_BOUNDARY = 3
R_BOUNDARY = 4
B_BOUNDARY = 5
T_BOUNDARY = 6
DIR_BOUNDARY = 7
OUT_BOUNDARY = 8

xmin = -L/2
ymin = -H/2
dx = L
dy = H
if symmetric:
   xmin = 0.
   ymin = 0.
   dx = L/2
   dy = H/2
mesh = gf.Mesh('import', 'structured',
               'GT="%s";ORG=[%f,%f];SIZES=[%f,%f];NSUBDIV=[%i,%i]'
               % (geotrans, xmin, ymin, dx, dy, N_L, N_H))

N = mesh.dim()

outer_faces = mesh.outer_faces()
outer_normals = mesh.normal_of_faces(outer_faces)
left_boundary   = outer_faces[:,np.nonzero(outer_normals[0] < -0.95)[0]]
right_boundary  = outer_faces[:,np.nonzero(outer_normals[0] >  0.95)[0]]
bottom_boundary = outer_faces[:,np.nonzero(outer_normals[1] < -0.95)[0]]
top_boundary    = outer_faces[:,np.nonzero(outer_normals[1] >  0.95)[0]]
mesh.set_region(L_BOUNDARY, left_boundary)
mesh.set_region(R_BOUNDARY, right_boundary)
mesh.set_region(B_BOUNDARY, bottom_boundary)
mesh.set_region(T_BOUNDARY, top_boundary)

mesh.set_region(DIR_BOUNDARY, left_boundary)
mesh.set_region(DIR_BOUNDARY, right_boundary)

mesh.set_region(OUT_BOUNDARY, outer_faces)

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

mesh.export_to_vtk('%s/mesh.vtk' % resultspath)

# FEM
mf1 = gf.MeshFem(mesh, 1)
mfN = gf.MeshFem(mesh, N)
#if disp_fem_order == 2:
#   mf1.set_fem(gf.Fem('FEM_Q2_INCOMPLETE(2)'))
#   mfN.set_fem(gf.Fem('FEM_Q2_INCOMPLETE(2)'))
#else:
mf1.set_classical_fem(disp_fem_order)
mfN.set_classical_fem(disp_fem_order)

if symmetric:
   xdofsL = mfN.basic_dof_on_region(L_BOUNDARY)[0::N]
   ydofsB = mfN.basic_dof_on_region(B_BOUNDARY)[1::N]
   kept_dofs = np.setdiff1d(np.arange(mfN.nb_basic_dof()),
                            np.union1d(xdofsL,ydofsB))
   mfu = gf.MeshFem('partial', mfN, kept_dofs)
else:
   mfu = mfN

mfp = gf.MeshFem(mesh, 1)
mfp.set_classical_fem(press_fem_order)

mfksi = gf.MeshFem(mesh, 1)
mfksi.set_classical_fem(mult_fem_order)

mfout = gf.MeshFem(mesh)
mfout.set_classical_discontinuous_fem(disp_fem_order)

mim = gf.MeshIm(mesh, integration_degree)
mimd1 = gf.MeshImData(mim, -1, [1])
mimd4 = gf.MeshImData(mim, -1, [4])

mimpress = gf.MeshIm(mesh, integration_degree_press)


# Model
md = gf.Model('real')

md.add_fem_variable('u', mfu)        # displacements
md.add_fem_variable('p', mfp)        # pressure
md.add_fem_variable("ksi", mfksi)    # plastic multiplier for the current step

md.add_im_data("gamma0", mimd1)      # accumulated plastic strain at previous time step
md.add_im_data('invCp0', mimd4)      # Components 11, 22, 33, 12, 13, 23 of the plastic part
md.set_variable('invCp0',            # of the inverse right Cauchy Green tensor at the
                np.tile([1,1,1,0],   # previous step.
                mimd4.nbpts()))      #

if symmetric:
   if rigid_grip:
      md.add_fem_data('dirichlet_data', mfN);
      ibdir = md.add_Dirichlet_condition_with_multipliers(mim, 'u', mfN, R_BOUNDARY, 'dirichlet_data')
   else:
      dirichlet_str = dirichlet_str[dirichlet_str.find('[')+1:dirichlet_str.find(',')]
      md.add_fem_data('dirichlet_data', mf1);
      ibdir = md.add_normal_Dirichlet_condition_with_multipliers(mim, 'u', mf1, R_BOUNDARY, 'dirichlet_data')
else:
   if not rigid_grip:
      print('Not rigid grip option is ignored in non-symmetric model')
   md.add_fem_data('dirichlet_data', mfN);
   ibdir = md.add_Dirichlet_condition_with_multipliers(mim, 'u', mfN, DIR_BOUNDARY, 'dirichlet_data')


gf.asm_define_function("hardening_law", 1,
                       "({sigma0}+{dsigma}*(1-exp(-{delta}*t))+{H}*t)"\
                       .format(sigma0=np.sqrt(2./3.)*sigma_0,
                               dsigma=np.sqrt(2./3.)*(sigma_inf-sigma_0),
                               delta=np.sqrt(2./3.)*h1,
                               H=2./3.*h2));

_K_ = "{0:.17g}".format(K)  # Bulk modulus
_mu_ = "{0:.17g}".format(mu)             # Shear modulus
ibplast = \
md.add_finite_strain_elastoplasticity_brick\
  (mim, "Simo_Miehe", "displacement_and_plastic_multiplier_and_pressure",
   "u", "ksi", "p", "gamma0", "invCp0",
   _K_, _mu_, "hardening_law")

print('Displacement dofs: %i' % mfu.nbdof())
print('Total model dofs: %i' % md.nbdof())

dirichlet_multname = md.mult_varname_Dirichlet(ibdir)
mf_dirichlet_mult = md.mesh_fem_of_variable(dirichlet_multname)
mass_mat = gf.asm_mass_matrix(mim, mfout)

shutil.copyfile(os.path.abspath(sys.argv[0]),resultspath+'/'+sys.argv[0])
starttime_overall = time.clock()
with open('%s/finite_strain_plasticity_forces.dat' % resultspath, 'w') as f:
   for step in range(0,steps+1):
      print('Solving step: %i' % step)

      backup = md.from_variables()
      dirichlet_data = mf_dirichlet_mult.eval(dirichlet_str.format(step),globals(),locals())
      md.set_variable('dirichlet_data', dirichlet_data)

      starttime = time.clock()
      iters = 41
      nit = \
      md.solve('noisy', 'lsolver', 'mumps', 'max_iter', iters, 'max_res', 1e-10, #)[0]
               'lsearch', 'simplest', 'alpha max ratio', 1.5, 'alpha min', 0.05, 'alpha mult', 0.6)[0]
      if nit > 40:
        print('STEP %i HAS NOT CONVERGED, TRYING SUBSTEPPING' % (step))
        md.to_variables(backup)
        for substep in range(1,5):
          dirichlet_data = mf_dirichlet_mult.eval(dirichlet_str.format(step-1.+substep/4.),globals(),locals())
          md.set_variable('dirichlet_data', dirichlet_data)
          md.solve('noisy', 'lsolver', 'mumps', 'max_iter', iters, 'max_res', 1e-10, #)[0]
                   'lsearch', 'simplest', 'alpha max ratio', 1.5, 'alpha min', 0.05, 'alpha mult', 0.6)[0]
      print('STEP %i COMPLETED IN %f SEC' % (step, time.clock()-starttime))

      VM = md.compute_finite_strain_elastoplasticity_Von_Mises\
         (mim, mfout, "Simo_Miehe", "displacement_and_plastic_multiplier_and_pressure",
          "u", "ksi", "p", "gamma0", "invCp0", _K_, _mu_, "hardening_law")

      RHS = md.rhs()
      IU = md.interval_of_variable('u')
      IP = md.interval_of_variable('p')
      IKSI = md.interval_of_variable('ksi')

      output = (mfu, md.variable('u'), 'Displacements',
                mfp, md.variable('p'), 'Pressure',
                mfksi, md.variable('ksi'), 'dksi',
                mf_dirichlet_mult, -md.variable(dirichlet_multname), "Reaction stresses",
                mfout, VM, 'Von Mises Stress',
                mfu, RHS[IU[0]:IU[0]+IU[1]], 'residual u',
                mfp, RHS[IP[0]:IP[0]+IP[1]], 'residual p',
                mfksi, RHS[IKSI[0]:IKSI[0]+IKSI[1]], 'residual ksi')

      md.finite_strain_elastoplasticity_next_iter\
         (mim, "Simo_Miehe", "displacement_and_plastic_multiplier_and_pressure",
          "u", "ksi", "p", "gamma0", "invCp0", _K_, _mu_, "hardening_law")

      GAMMA = gf.asm_generic(mim, 1, "gamma*Test_t", -1,
                             "gamma", False, mimd1, md.variable("gamma0"),
                             "t", True, mfout, np.zeros(mfout.nbdof()))
      GAMMA = np.transpose(gf.linsolve_mumps(mass_mat, GAMMA))
      output += (mfout, GAMMA, "plastic strain")

      mf1.export_to_vtk('%s/finite_strain_plasticity_%i.vtk' % (resultspath, step), 'ascii', *output)

      DIRMULT = -md.variable(dirichlet_multname)
      DIRMULT = np.reshape(DIRMULT,[1,DIRMULT.size])
      if symmetric and not rigid_grip:
         dfR = gf.asm_boundary_source(R_BOUNDARY, mim, mf1, mf_dirichlet_mult, DIRMULT)
         f.write(('step=%i fR=(%f,0) sR=(%f,0)\n') %
                 (step, dfR.sum(), dfR.sum()/H))
      else:
         dfR = gf.asm_boundary_source(R_BOUNDARY, mim, mfN, mf_dirichlet_mult, DIRMULT)
         f.write(('step=%i fR=(%f,%f) sR=(%f,%f)\n') %
                 (step, dfR[0::N].sum(), dfR[1::N].sum(),
                        dfR[0::N].sum()/H, dfR[1::N].sum()/H))
      f.flush()

print('OVERALL SOLUTION TIME IS %f SEC' % (time.clock()-starttime_overall))
