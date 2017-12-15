#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Python GetFEM++ interface
#
# Copyright (C) 2017-2018 Yves Renard, Franz Chouly.
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

""" Dynamic with impact of an elastic road, comparison with the exact solution

  This program is used to check that python-getfem is working. This is also
  a good example of use of GetFEM++.
"""

import getfem as gf
import numpy as np


# Parameters
NX=20                # Number of elements
T = 10               # Simulation time
dt = 0.001           # Time step
u_degree = 1         # Degree of the finite element method for u

gamma0_N = 0.05      # Nitsche parameter gamma
theta_N = 0          # Nitsche parameter theta

beta = 0.25          # Newmark parameter beta
gamma = 0.5          # Newmark parameter gamma
theta = 1.0          # Theta-scheme parameter

singular_mass = 0    # singular mass or not

scheme = 0
version = 0;         # 0 = penalized contact ...

# Mesh
m=gf.Mesh('cartesian', np.arange(0,1+1./NX,1./NX))

# Selection of the contact and Dirichlet boundaries
GAMMAC = 1; GAMMAD = 2

border = m.outer_faces()
normals = m.normal_of_faces(border)
contact_boundary = border[:,np.nonzero(normals[0] < -0.01)[0]]
m.set_region(GAMMAC, contact_boundary)
contact_boundary = border[:,np.nonzero(normals[0] > 0.01)[0]]
m.set_region(GAMMAD, contact_boundary)

# Finite element methods
mfu = gf.MeshFem(m)
mfu.set_classical_fem(u_degree)

mfd = gf.MeshFem(m, 1)
mfd.set_classical_fem(u_degree)

# Integration method
mim = gf.MeshIm(m, 4)

# Mass and stiffness matrices
md = gf.Model('real'); md.add_fem_variable('u', mfu);
M = gf.asm_generic(mim, 2, 'u*Test_u', -1, md)
K = gf.asm_generic(mim, 2, 'Grad_u*Grad_Test_u', -1, md)

# Dirichlet condition on the top
for i in range(0, u_degree+1): K[0,i] = 0.;
for i in range(0, u_degree+1): M[0,i] = 0.;
M[0,0] = K[0,0] = 1.;

print uExact(0,0)


print K
print K[1,1]




def uExact(x, t):
    # The solution is periodic of period 3
    # Shift the time 't' with t=0 : beginning of the period
    tp = rem(t,3.)
    # The solution has 3 phases
    # Shift the time 'tp' with t=0 : beginning of a phase
    # and get also the phase number
    tf = rem(tp,1.)
    nf = floor(tp)
    # Get the index of the zone in each phase : I, II, III, IV
    # (zones are given according to characteristics of the wave equation)
    if (tf<=x):
        if (tf<=(1-x)): zone = 1; else: zone = 3;
    else:
        if (tf<=(1-x)): zone = 2; else: zone = 4;
    # Solution according to the Phase (1,2,3) and the zone (I, II, III, IV)
    if nf == 0:
        if   zone == 1: u = 1./2.-x/2.;   dxu = -1./2.;  dtu = 0.;
        elif zone == 2: u = 1./2.-tf/2.;  dxu = 0.;      dtu = -1./2.;
        elif zone == 3: u = 1./2.-x/2.;   dxu = -1./2.;  dtu = 0;
        elif zone == 4: u = 1./2.-tf/2.;  dxu = 0;       dtu = -1./2.;
    elif nf == 1: 
        if   zone == 1: u = -tf/2.;       dxu = 0;       dtu = -1./2.
        elif zone == 2: u = -x/2.;        dxu = -1./2.;  dtu = 0;
        elif zone == 3: u = -1./2.+x/2.;  dxu = 1./2.;   dtu = 0;
        elif zone == 4: u = -1./2.+tf/2.; dxu = 0;       dtu = 1./2.;
    elif nf == 2:
        if   zone == 1: u = tf/2.;        dxu = 0;       dtu = 1./2.;
        elif zone == 2: u = tf/2.;        dxu = 0;       dtu = 1./2.;
        elif zone == 3: u = 1./2.-x/2.;   dxu = -1./2.;  dtu = 0;
        elif zone == 4: u = 1./2.-x/2.;   dxu = -1./2.;  dtu = 0.
    return u

