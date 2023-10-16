#!/usr/bin/env python3
# -*- coding: UTF8 -*-
# Python GetFEM interface
#
# Copyright (C) 2021-2021 Konstantinos Poulios.
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

"""
This example provides an implementation of large strain quasilinear
viscoelasticity according to (de Pascalis, Abrahams, Parnell 2014). It
implements the compressible case with the Horgan-Murphy law from the paper,
as well as a neo-Hookean law from (Simo, Taylor, Pister 1985). The
implementation only covers Prony stress relaxation functions, which allow
a much simpler time integration, requiring storage of data at only one
previous time step.
"""

import getfem as gf

gf.util_trace_level(1)
gf.util_warning_level(1)

# Material parameters under high strain rate
E =   1e2                    # Young's modulus [Pa]
nu =  0.49                   # Poisson's ratio [-]
kappa = E/(3*(1-2*nu))
mu = E/(2*(1+nu))
MAT_LAW = ('neohookean-Simo','Horgan-Murphy')[1]
gamma = 1./6.  # only used for Horgan-Murphy

# Viscoelastic material constants
tauH = 1e-3
tauD = 1e-3
rH = 0.5 # the low strain rate bulk modulus is rH*kappa
rD = 0.5 # the low strain rate shear modulus is rD*mu

# Dimensions [mm]
W = 0.5     # Block width
H = 2       # Block height
T = 0.5     # Block thikness

# Time
t_max = 28e-3

# Numerical parameters
N_W = 1
N_H = 1
N_T = 1
steps = 112
fem_order = 2

# Mesh
m = gf.Mesh("import",
            "structured","GT='GT_QK(3,2)'; ORG=[0,0,0]; SIZES=[%f,%f,%f]; NSUBDIV=[%i,%i,%i]"
            % (W, H, T, N_W, N_H, N_T))

# Boundaries
RG_LEFT = 11
RG_BOTTOM = 12
RG_BACK = 13
RG_TOP = 14
m.set_region(RG_LEFT, m.outer_faces_in_box([-1e-5,-1e-5,-1e-5], [1e-5,H+1e-5,T+1e-5]))
m.set_region(RG_BOTTOM, m.outer_faces_in_box([-1e-5,-1e-5,-1e-5], [W+1e-5,1e-5,T+1e-5]))
m.set_region(RG_BACK, m.outer_faces_in_box([-1e-5,-1e-5,-1e-5], [W+1e-5,H+1e-5,1e-5]))
m.set_region(RG_TOP, m.outer_faces_in_box([-1e-5,H-1e-5,-1e-5], [W+1e-5,H+1e-5,T+1e-5]))


mfu_ = gf.MeshFem(m, 3)
mfu_.set_classical_fem(fem_order)

kept_dofs = list(range(mfu_.nbdof()))
# remove x-dofs on RG_LEFT, y-dofs on RG_BOTTOM, and z-dofs on RG_BACK
for skipped_dof, RG in enumerate((RG_LEFT,RG_BOTTOM,RG_BACK)):
  for dof in mfu_.basic_dof_on_region(RG):
    if dof % 3 == skipped_dof:
      kept_dofs.remove(dof)
mfu = gf.MeshFem('partial', mfu_, kept_dofs)

mfmult = gf.MeshFem(m, 1)
mfmult.set_classical_fem(fem_order)

mfout = gf.MeshFem(m, 1)
mfout.set_classical_discontinuous_fem(2)

mim = gf.MeshIm(m, 5)

# Model
md = gf.Model('real')
md.add_fem_variable('u', mfu)

mimd3x3 = gf.MeshImData(mim, -1, [3,3])

md.add_initialized_data('epsYY', 0)

md.add_initialized_data('kappa', kappa)
md.add_initialized_data('mu', mu)

md.add_initialized_data('tauH', tauH)
md.add_initialized_data('tauD', tauD)
md.add_initialized_data('rH', rH)
md.add_initialized_data('rD', rD)

dt = t_max/steps
md.add_initialized_data('dt', dt)

md.add_im_data("SeHprev", mimd3x3)   # Elastic hydrostatic stress at previous timestep
md.add_im_data("SeDprev", mimd3x3)   # Elastic deviatoric stress at previous timestep
md.add_im_data("SvHprev", mimd3x3)   # Viscous hydrostatic stress at previous timestep
md.add_im_data("SvDprev", mimd3x3)   # Viscous deviatoric stress at previous timestep

md.add_macro("F", "Id(3)+Grad_u")
md.add_macro("Fprev", "Id(3)+Grad_Previous_u")

if MAT_LAW == 'neohookean-Simo':
  md.add_macro("SeH", "kappa*log(Det(F))*Inv(Right_Cauchy_Green(F))")
  md.add_macro("SeD", "mu*pow(Det(F),-2/3)*Inv(F)*Deviator(Left_Cauchy_Green(F))*Inv(F)'")
elif MAT_LAW == 'Horgan-Murphy':
  md.add_initialized_data('gamma', gamma)
  md.add_macro("Sqr(a)", "a*a");
  md.add_macro("WW1", "mu*(0.5+gamma)")
  md.add_macro("WW2", "mu*(0.5-gamma)")
  md.add_macro("WW3", "kappa*(1-1/Det(F)) - 2*WW2*pow(Det(F),-2/3) - WW1*pow(Det(F),-4/3)")
  md.add_macro("I1", "Norm_sqr(F)")
  md.add_macro("I2", "Sqr(Det(F))*Norm_sqr(Inv(F))")
  md.add_macro("I3", "Sqr(Det(F))")
  md.add_macro("SeH", "(1/3*WW1*I1 + 2/3*WW2*I2 + WW3*I3)*Inv(Right_Cauchy_Green(F))")
  md.add_macro("SeD", "WW1*Id(3) + 1/3*(WW2*I2-WW1*I1)*Inv(Right_Cauchy_Green(F))"
                                "- WW2*I3*Sqr(Inv(Right_Cauchy_Green(F)))")

md.add_macro("SvH", "exp(-dt/tauH)*SvHprev"
                    "-(1-rH)*(tauH/dt-(dt+tauH)/dt*exp(-dt/tauH))*SeHprev"
                    "-(1-rH)*((dt-tauH)/dt+tauH/dt*exp(-dt/tauH))*SeH")
#                    "-(1-rH)*(1-exp(-dt/tauH))*SeHprev"
#                    "-(1-rH)*(1-tauH/dt*(1-exp(-dt/tauH)))*(SeH-SeHprev)")
md.add_macro("SvD", "exp(-dt/tauD)*SvDprev"
                    "-(1-rD)*(tauD/dt-(dt+tauD)/dt*exp(-dt/tauD))*SeDprev"
                    "-(1-rD)*((dt-tauD)/dt+tauD/dt*exp(-dt/tauD))*SeD")
#                    "-(1-rD)*(1-exp(-dt/tauD))*SeDprev"
#                    "-(1-rD)*(1-tauD/dt*(1-exp(-dt/tauD)))*(SeD-SeDprev)")

# Virtual work expression
md.add_nonlinear_term(mim, "(F*(SeH+SvH+SeD+SvD)):Grad_Test_u")

# Cauchy and von Mises definitions
md.add_macro("Cauchy", "(F*(SeH+SvH+SeD+SvD)*F')/Det(F)")
md.add_macro("VM", "sqrt(1.5)*Norm(Deviator(F*(SeH+SvH+SeD+SvD)*F'))/Det(F)")

# Dirichlet condition
md.add_filtered_fem_variable("dirmult", mfmult, RG_TOP)
md.add_linear_term(mim, "(epsYY*X(2)-u(2))*dirmult", RG_TOP)

Vol = gf.asm_generic(mim, 0, "1", -1, md)
with open("QLV_viscoelasticity_results.dat", "w") as f:
  for i in range(steps+1):
    t = i*dt
    epsYY = 0.3*max(0,min(t/(2.*t_max/7.),((6.*t_max/7.)-t)/(2.*t_max/7.),1.))
    md.set_variable("epsYY", epsYY)

    md.solve('noisy', 'lsolver', 'mumps', 'max_res', 1E-7,'max_iter',20,
             'lsearch', 'simplest', 'alpha max ratio', 1., 'alpha min', 0.2, 'alpha mult', 0.6)

    sigmayy = gf.asm_generic(mim, 0, "Cauchy(2,2)", -1, md) / Vol
    f.write("t=%.10g sigma_yy=%.10g\n" % (t,sigmayy))

    U = md.variable('u')
    VM = md.local_projection(mim, "VM", mfout)
    SIGMA11 = md.local_projection(mim, "Cauchy(1,1)", mfout)
    SIGMA22 = md.local_projection(mim, "Cauchy(2,2)", mfout)
    SIGMA33 = md.local_projection(mim, "Cauchy(3,3)", mfout)
    mfout.export_to_vtk('QLV_viscoelasticity%i.vtk' % i, mfu, U, 'Displacements',
                                                     mfout, VM,  'Von Mises Stresses',
                                                     mfout, SIGMA11,  'Sigma 11',
                                                     mfout, SIGMA22,  'Sigma 22',
                                                     mfout, SIGMA33,  'Sigma 33')

    md.set_variable("SeHprev", md.interpolation("SeH", mimd3x3, -1))
    md.set_variable("SeDprev", md.interpolation("SeD", mimd3x3, -1))
    md.set_variable("SvHprev", md.interpolation("SvH", mimd3x3, -1))
    md.set_variable("SvDprev", md.interpolation("SvD", mimd3x3, -1))

