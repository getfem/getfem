#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Python GetFEM interface
#
# Copyright (C) 2017-2020 Yves Renard, Franz Chouly.
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

""" Dynamic with impact of an elastic rod.
    Comparison to the closed-form solution.

    This program is used to check that python-getfem is working.
    This is also a good example of use of GetFEM.
"""

import os
import sys

import matplotlib.pyplot as plt
import numpy as np

import getfem as gf

# Numerical parameters
NX = 20               # Number of elements
T = 12                # Simulation time
dt = 0.002            # Time step
u_degree = 1          # Degree of the finite element method for u

gamma0_N = 5.         # Nitsche parameter gamma
theta_N =  0.         # Nitsche parameter theta

gamma0_P = 1.         # Penalization parameter

beta = 0.             # Newmark parameter beta
gamma = 0.5           # Newmark parameter gamma

e_PS = 0.             # Restitution coefficient for Paoli-Schatzman scheme

version = 4           # 0 = pure Signorini contact
                      # 1 = pure Signorini contact with Paoli-Schatzman scheme
                      # 2 = penalized contact
                      # 3 = Nitsche's method
                      # 4 = Taylor-Flanagan method

lump_mass_matrix = 0  # 0 = standard mass matrix
                      # 1 = basic lumped mass matrix

mass_matrix_type = 0  # 0 = standard mass matrix
                      # 1 = redistributed mass matrix
                      # 2 = singular dynamic mass matrix


# Output parameters
dtplot = 0.05         # Time step for intermediate plots
do_inter_plot = False # Intermediate plots or not
do_final_plot = True  # Final plots or not
output_directory = './expe_num'
root_filename = 'dyn1d'
do_export_in_files = False;

# Read optional parameters on the command line 
for i in range(1,len(sys.argv)): exec(sys.argv[i])

print("Begin experiment for", end=' ')
if    (version == 0): print("Pure Signorini contact", end=' ')
elif  (version == 1): print("Paoli-Schatzman scheme", end=' ')
elif  (version == 2): print("Penalized contact", end=' ')
elif  (version == 3): print("Nitsche's method", end=' ')
elif  (version == 4): print("Taylor-Flanagan method", end=' ')
print(" in P%d, with NX = %d, dt = %g" % (u_degree,NX, dt))

if (version == 4 and beta != 0): print('Incompatibility'); exit(1)

# Deduced parameters
h = 1./NX
TT = np.arange(0, T+dt, dt)
NT = TT.size
dt_max_approx = h/(2* u_degree);
if (version == 2): dt_max_approx = min(dt_max_approx, 2*h/(gamma0_P))
if (version == 3): dt_max_approx = min(dt_max_approx, 2*h/(gamma0_N))
print('Approximative dt_max for CFL :', dt_max_approx)
if (beta == 0 and dt > dt_max_approx): print('Time step too large'); exit(1)

# Exact solution. The solution is periodic of period 3
# Return the displacement (d=0), x derivative (d=1) or time derivative (d=2)
def uExact(x, t, d = 0):
    # Shift the time 't' with t=0 : beginning of the period
    tp = t % 3.
    # The solution has 3 phases
    # Shift the time 'tp' with t=0 : beginning of a phase
    # and get also the phase number
    tf = tp % 1.
    nf = np.floor(tp)
    # Get the index of the zone in each phase : I, II, III, IV
    # (zones are given according to characteristics of the wave equation)
    if (tf<x): zone = 1 if (tf<=(1-x)) else 3
    else:      zone = 2 if (tf<=(1-x)) else 4
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
    return (u if (d == 0) else (dxu if (d == 1) else dtu))


def linsolve(M, B): # Call Superlu to solve a sparse linear system
    return (((gf.linsolve_superlu(M, B))[0]).T)[0]

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
mfu.set_classical_fem(u_degree) # Assumed to be a Lagrange FEM in the following

mfd = gf.MeshFem(m, 1)
mfd.set_classical_fem(u_degree)

# Integration method
mim = gf.MeshIm(m, 4)

# GetFEM model
md = gf.Model('real'); md.add_fem_variable('u', mfu)
md.add_fem_data('v', mfu)
md.add_initialized_data('t_N', theta_N)
md.add_initialized_data('g_N', gamma0_N/h)

# Initial conditions
U0 = mfu.eval('0.5-0.5*x')   # Initial displacement
Um1 = np.copy(U0)            # U_{-1} for Paoli-Schatzman scheme
Ndof = U0.size
V0 = np.zeros(Ndof)          # Initial velocity
s0 = 0.                      # Initial stress
md.set_variable('u', U0); md.set_variable('v', V0)

# Mass and stiffness matrices
M = gf.asm_generic(mim, 2, 'u*Test_u', -1, md)
if (lump_mass_matrix == 1):
    assert (u_degree == 1), "Sorry, basic lump only for affine elements"
    for j in range(0, Ndof):
        for i in range(1, u_degree+1):
            if (j+i < Ndof): M[j,j] += M[j,j+i]; M[j,j+i] = 0.;
            if (j-i >= 0): M[j,j] += M[j,j-i]; M[j,j-i] = 0.;
K = gf.asm_generic(mim, 2, 'Grad_u*Grad_Test_u', -1, md)

# Dirichlet condition on the top
for i in range(0, u_degree+1): K[Ndof-1,Ndof-1-i] = 0.
for i in range(0, u_degree+1): M[Ndof-1,Ndof-1-i] = 0.
M[Ndof-1,Ndof-1] = 1.; K[Ndof-1,Ndof-1] = 1.;
M2 = gf.Spmat('copy', M);

if (mass_matrix_type == 1): # Redistributed mass matrix
    M[1,1] += M[0,0]; M[0,0] = 0.;
    for i in range(1, u_degree+1):
        M[i,i] += M[i,0]; M[0,i] = M[i,0] = 0.
    M2 = gf.Spmat('copy', M); M2[0,0] = 1.
elif (mass_matrix_type == 2): # Singular dynamic mass matrix
    assert (u_degree == 1), 'Sorry, implemented for linear element only'
    M[1,1] = 7.*h/12.;
    for i in range(0, u_degree+1): M[0,i] = M[i,0] = 0.
    M2 = gf.Spmat('copy', M); M2[0,0] = 1.

# Matrices for Newmark method
MV0 = M.mult(V0); MV1 = np.copy(MV0)
MA0 = -K.mult(U0);
if (mass_matrix_type >= 1): MA0[0] = 0;

K_m = gf.Spmat('copy', K); K_m.scale(dt*dt*beta);
K_N = gf.Spmat('add', M, K_m);  # Newmark matrix without contact
K_N_C = gf.Spmat('copy', K_N);  # Newmark matrix with effective contact

if (version == 0):     # Pure Signorini contact
    for i in range(0, u_degree+1): K_N_C[0,i] = 0.
    K_N_C[0,0] = 1;
elif (version == 1 or version == 4):
    # Pure Signorini contact with Paoli-Schatzman scheme
    for i in range(0, u_degree+1): K_N_C[0,i] = 0.
    K_N_C[0,0] = 1;
elif (version == 2):   # Penalized contact
    K_N_C[0,0] += dt*dt*beta*(gamma0_P/h)
elif (version == 3):   # Nitsche contact
    asstr_Nitsche = '-(t_N/g_N)*Grad_u*Grad_Test_u'
    K_N1 = gf.asm_generic(mim, 2, asstr_Nitsche, GAMMAC, md)
    K_Nitsche1 = gf.Spmat('add', K, K_N1)
    K_N1.scale(dt*dt*beta);
    K_N = gf.Spmat('add', K_N, K_N1)
    asstr_Nitsche = '(1/g_N)*(g_N*u+Grad_u)*(g_N*Test_u+t_N*Grad_Test_u)'
    K_N2 = gf.asm_generic(mim, 2, asstr_Nitsche, GAMMAC, md)
    K_Nitsche2 = gf.Spmat('add', K_Nitsche1, K_N2)
    K_N2.scale(dt*dt*beta);
    K_N_C = gf.Spmat('add', K_N, K_N2)
    asstr_Nitsche = '(t_N/g_N)*Grad_u*Grad_Test_u ' \
             + '+(1/g_N)*neg_part(g_N*u+Grad_u)*(g_N*Test_u+t_N*Grad_Test_u)'
else: assert(false), "Unvalid version"
              
# Misc initializations
Ndof = U0.size
Fc = np.zeros(Ndof); Uex = np.zeros(Ndof); Vex = np.zeros(Ndof); 
Xdraw = np.arange(0, 1+0.001, 0.001) # Grid for plot
Ndraw = Xdraw.size
tplot = 0
store_E   = np.zeros(NT);
store_VL2 = np.zeros(NT); store_VL2_ex = np.zeros(NT)
store_UL2 = np.zeros(NT); store_UL2_ex = np.zeros(NT)
store_UH1 = np.zeros(NT); store_UH1_ex = np.zeros(NT)
store_u0  = np.zeros(NT); store_u0_ex  = np.zeros(NT)
store_v0  = np.zeros(NT); store_v0_ex  = np.zeros(NT)
store_s0  = np.zeros(NT); store_s0_ex  = np.zeros(NT)
fem_nodes = mfu.basic_dof_nodes(); s0 = s1 = ux = 0.;

# Time steps
for nit in range(0, NT):
    t = TT[nit];

    if (t > 0.):
        if (beta == 0): # Newmark Explicit scheme (beta = 0)
            F = Fc - K.mult(U0)
            FF = linsolve(M2, F)
            U1 = U0 + dt*V0 + (dt*dt/2.)*FF
            V1 = V0 + (dt*(1.-gamma))*FF
            a = 0;
            for i in range(1, u_degree+1): a += K[0,i]*U1[i]
            if (version == 0):   # Pure Signorini contact
                if (mass_matrix_type >= 1):
                    U1[0] = max(0., -a) / K[0,0];
                    s1 = -max(0., a); Fc[0] = -s1
                else: assert (false), 'No explicit scheme for ' \
                        + 'pure Signorini contact and standard mass matrix'
            elif (version == 1): # Pure Signorini contact with P.-S. scheme
                if (mass_matrix_type >= 1):
                    U1[0] = max(-e_PS*Um1[0], -a) / K[0,0];
                    s1 = min(-e_PS*Um1[0], -a);
                else:
                    s1 = 0;
                    if ((1./(1.+e_PS))*U1[0] + (e_PS/(1.+e_PS))*Um1[0] < 0):
                        U1_0 = np.copy(U1);
                        B=M.mult(U0)+dt*M.mult(V0)-(dt*dt/2.)*K.mult(U0)
                        B0 = np.copy(B); B0[0] = -e_PS*Um1[0];
                        U1=linsolve(K_N_C, B0)
                        F=M.mult(U1)-B; s1 = -F[0]/(dt*dt)
                        V1 += (U1-U1_0)/dt; 
                Fc[0] = 0;
            elif (version == 4): # Pure Signorini with Taylor-Flanagan scheme
                if (mass_matrix_type >= 1):
                    U1[0] = max(0., -a) / K[0,0];
                    s1 = min(0., -a);
                else:
                    s1 = 0;
                    if (U1[0] < 0):
                        F *= 0.; F[0] = 1;
                        F2=linsolve(M, F);
                        s1 = V1[0]/(F2[0]*dt);
                        F *= 0.; F[0] = -s1;
                        F2=linsolve(M, F);
                        U1 += dt*dt*F2;
                        V1 += dt*F2;
            elif (version == 2): # Penalized contact
                if (mass_matrix_type >= 1):
                    U1[0] = -a / K[0,0];
                    if (U1[0] < 0): U1[0] = -a / (K[0,0] + gamma0_P/h);
                s1 = (gamma0_P/h)*min(0., U1[0]); Fc[0] = -s1
            else:                # Nitsche contact
                if (mass_matrix_type >= 1):
                    a = 0;
                    for i in range(1, u_degree+1): a += K_Nitsche1[0,i]*U1[i]
                    U1[0] = -a / K_Nitsche1[0,0]
                    if (U1[0] < 0):
                        a = 0;
                        for i in range(1,u_degree+1): a += K_Nitsche2[0,i]*U1[i]
                        U1[0] = -a / K_Nitsche2[0,0]
                md.set_variable('u', U1); ux = md.interpolation('Grad_u',[0.],m)
                Fc = gf.asm_generic(mim, 1, asstr_Nitsche, GAMMAC, md)
                s1 = min(0.,(gamma0_N/h)*U1[0]+ux)
            FF = linsolve(M2, Fc-K.mult(U1))
            V1 += dt*gamma*FF;
            if (mass_matrix_type >= 1): V1[0] = (U1[0] - U0[0])/dt
        else:  # Implicit scheme (beta > 0)
            B = M.mult(U0)+dt*MV0+(dt*dt*(1-2*beta)/2.)*MA0
            U1 = linsolve(K_N, B)
            MV1 = MV0
            if (version == 0):   # Pure Signorini contact
                if (U1[0] < 0):
                    B2 = np.copy(B); B2[0]=0; U1=linsolve(K_N_C, B2)
                    Fc1 = (K_N.mult(U1) - B)/(beta*dt*dt); s1 = -Fc1[0]
                else: Fc1 = 0.*Fc; s1 = 0
            elif (version == 1): # Pure Signorini contact with P.-S. scheme
                if ((1./(1.+e_PS))*U1[0] + (e_PS/(1.+e_PS))*Um1[0] < 0):
                    B2 = np.copy(B); B2[0]= -e_PS*Um1[0];
                    U1 = linsolve(K_N_C, B2)
                    Fc1 = (K_N.mult(U1) - B)/(dt*dt); s1 = -Fc1[0]
                    MV1 += dt*Fc1
                else: s1 = 0;
                Fc1 = 0.*Fc;
            elif (version == 2): # Penalized contact
                if (U1[0] < 0):
                    U1 = linsolve(K_N_C, B)
                    Fc1 = (K_N.mult(U1) - B)/(beta*dt*dt); s1 = -Fc1[0];
                else: Fc1 = 0.*Fc; s1 = 0
            else:                # Nitsche contact
                md.set_variable('u', U1); ux = md.interpolation('Grad_u',[0.],m)
                if ((gamma0_N/h)*U1[0]+ux < 0):
                    U1 = linsolve(K_N_C, B); md.set_variable('u', U1);
                    ux = md.interpolation('Grad_u',[0.],m)
                Fc1 = (K_N.mult(U1) - B)/(beta*dt*dt)
                s1 = min(0.,(gamma0_N/h)*U1[0]+ux)
            MA1 = Fc1-K.mult(U1);
            MV1 += dt*((1.-gamma)*MA0 + gamma*MA1);
            V1 = linsolve(M2, MV1); V1[Ndof-1] = 0.
            if (mass_matrix_type >= 1): V1[0] = (U1[0] - U0[0])/dt
            MA0 = MA1
        
        # End of time step
        Um1 = np.copy(U0); U0 = np.copy(U1);
        V0 = np.copy(V1); MV0 = np.copy(MV1); s0 = s1
        md.set_variable('u', U0); md.set_variable('v', V0)

    # Compute the difference with the exact solution
    for i in range(0,Ndof): Uex[i] = uExact(fem_nodes[0][i], t)
    for i in range(0,Ndof): Vex[i] = uExact(fem_nodes[0][i], t, 2)

    # Compute the Energy
    MMV0 = M.mult(V0); KU0 = K.mult(U0)
    E = (np.vdot(MMV0, V0) + np.vdot(KU0, U0))*0.5
    E_org = E
    if (version == 2): # Energy stored in the penalized contact
        E += (gamma0_P/h)*pow(min(0., U0[0]),2)/2
    if (version == 3): # Nitsche Energy
        E += (0.5*theta_N)*(h/gamma0_N)*(s0*s0-ux*ux);
    
    # Store the data to be ploted at the end
    store_u0[nit] = U0[0]; store_u0_ex[nit] = uExact(0., t)
    store_v0[nit] = V0[0]; store_v0_ex[nit] = uExact(0., t, 2)
    store_s0[nit] = s0
    store_s0_ex[nit] = uExact(0., t, 1);
    store_E[nit] = E;
    store_VL2_ex[nit] = gf.compute_L2_norm(mfu, Vex, mim)
    store_UL2_ex[nit] = gf.compute_L2_norm(mfu, Uex, mim)
    store_UH1_ex[nit] = gf.compute_H1_norm(mfu, Uex, mim)
    store_VL2[nit] = gf.compute_L2_norm(mfu, V0-Vex, mim)
    store_UL2[nit] = gf.compute_L2_norm(mfu, U0-Uex, mim)
    store_UH1[nit] = gf.compute_H1_norm(mfu, U0-Uex, mim)

    
        
    # Draw the approximated and exact solutions
    if (t >= tplot-(1e-10)):
        tplot += dtplot;
        print(("Time %3f"% t), "/", T, end=' ')
        print((" Energy %7f" % E), (" Mech energy %7f" % E_org))
        
        if (do_inter_plot):
            UUex = np.copy(Xdraw)
            plt.figure(1); plt.rc('text', usetex=True)
            plt.subplot(311) # Displacement
            for i in range(0,Ndraw): UUex[i] = uExact(Xdraw[i], t)
            UU = md.interpolation('u', Xdraw, m)
            plt.cla(); plt.axis([0.,1.,-0.6,0.6])
            plt.plot(Xdraw, UUex, 'red'); plt.plot(Xdraw, UU, 'blue')
            plt.ylabel('u'); plt.xlabel('x')
            plt.subplot(312) # Velocity
            for i in range(0,Ndraw): UUex[i] = uExact(Xdraw[i], t, 2)
            UU = md.interpolation('v', Xdraw, m)
            plt.cla(); plt.axis([0.,1.,-0.6,0.6])
            plt.plot(Xdraw, UUex, 'red'); plt.plot(Xdraw, UU, 'blue')
            plt.ylabel('v'); plt.xlabel('x')
            plt.subplot(313) # Stress
            for i in range(0,Ndraw): UUex[i] = uExact(Xdraw[i], t, 1)
            UU = md.interpolation('Grad_u', Xdraw, m)
            plt.cla(); plt.axis([0.,1.,-0.6,0.6])
            plt.plot(Xdraw, UUex, 'red'); plt.plot(Xdraw, UU, 'blue')
            plt.ylabel('$\partial_x u$'); plt.xlabel('x')
            plt.pause(0.01); plt.show(0)
    

    
# print the main relative errors
LinfL2u = np.amax(store_UL2) / np.amax(store_UL2_ex)
print('L^\intfy(0,T,L^2)-norm of the error on u: ', LinfL2u)
LinfH1u = np.amax(store_UH1) / np.amax(store_UH1_ex)
print('L^\intfy(0,T,H^1)-norm of the error on u: ', LinfH1u)
LinfL2v = np.amax(store_VL2) / np.amax(store_VL2_ex)
print('L^\intfy(0,T,L^2)-norm of the error on v: ', LinfL2v)
Nor = np.sqrt(np.sum(np.square(store_UL2_ex))*dt)
L2L2u = np.sqrt(np.sum(np.square(store_UL2))*dt) / Nor
print('L^2(0,T,L^2)-norm of the error on u: ', L2L2u)
Nor = np.sqrt(np.sum(np.square(store_UH1_ex))*dt)
L2H1u = np.sqrt(np.sum(np.square(store_UH1))*dt) / Nor
print('L^2(0,T,H^1)-norm of the error on u: ', L2H1u)
Nor = np.sqrt(np.sum(np.square(store_VL2_ex))*dt)
L2L2v = np.sqrt(np.sum(np.square(store_VL2))*dt) / Nor
print('L^2(0,T)-norm of the error on v: ', L2L2v)
Nor = np.sqrt(np.sum(np.square(store_s0_ex))*dt)
L2sn = np.sqrt(np.sum(np.square(store_s0-store_s0_ex))*dt) / Nor
print('L^2(0,T)-norm of the error on contact stress: ', L2sn)


if (do_export_in_files):
    if (not os.path.exists(output_directory)):
        os.makedirs(output_directory)

    # Parameter file
    pf = open(output_directory+'/' + root_filename + '.params', 'w')
    pf.write('NX = %d;' % NX);            pf.write('T = %g;' % T)
    pf.write('dt = %g;' % dt);            pf.write('u_degree = %d;' % u_degree)
    pf.write('gamma0_N = %g;' % gamma0_N);pf.write('theta_N =  %g;' % theta_N)
    pf.write('gamma0_P = %g;' % gamma0_P);pf.write('beta = %g;' % beta)
    pf.write('gamma = %g;' % gamma);      pf.write('e_PS = %g;' % e_PS);
    pf.write('version = %d;' % version)
    pf.write('mass_matrix_type = %d;' % mass_matrix_type)
    pf.close()

    # Export solutions
    np.save(output_directory+'/' + root_filename + '.u0.npy', store_u0);
    np.save(output_directory+'/' + root_filename + '.v0.npy', store_v0);
    np.save(output_directory+'/' + root_filename + '.s0.npy', store_s0);
    np.save(output_directory+'/' + root_filename + '.E.npy' , store_E );
    
    # write error norms
    pf = open(output_directory+'/' + root_filename + '.err_norms', 'w')
    pf.write('LinfL2u = %g;' % LinfL2u); pf.write('LinfH1u = %g;' % LinfH1u)
    pf.write('LinfL2v = %g;' % LinfL2v); pf.write('L2L2u = %g;' % L2L2u)
    pf.write('L2H1u = %g;' % L2H1u);     pf.write('L2L2v = %g;' % L2L2v)
    pf.write('L2sn = %g;' % L2sn)
    pf.close()
    
    



if (do_final_plot):
    plt.figure(2); # Displacement evolution at the contact boundary
    plt.cla(); plt.axis([0.,T,-0.15,0.6])
    plt.plot(TT, store_u0_ex, 'red')
    plt.plot(TT, store_u0, 'blue')
    plt.ylabel('u(0)'); plt.xlabel('t')
    
    plt.figure(3); # velocity evolution at the contact boundary
    plt.cla(); plt.axis([0.,T,-0.6,0.6])
    plt.plot(TT, store_v0_ex, 'red')
    plt.plot(TT, store_v0, 'blue')
    plt.ylabel('v(0)'); plt.xlabel('t')

    plt.figure(4); # stress evolution at the contact boundary
    plt.cla(); plt.axis([0.,T,-0.6,0.6])
    plt.plot(TT, store_s0_ex, 'red')
    plt.plot(TT, store_s0, 'blue')
    plt.ylabel('contact pressure'); plt.xlabel('t')

    plt.figure(5); # Energy evolution
    plt.cla(); plt.axis([0.,T,0,0.6])
    plt.plot([0.,T], [1./8., 1./8.], 'red')
    plt.plot(TT, store_E, 'blue')
    plt.ylabel('Energy'); plt.xlabel('t')

    plt.figure(6);
    plt.subplot(311) # L2-norm of the error evolution
    plt.cla(); plt.axis([0.,T,0,0.6])
    plt.plot(TT, store_UL2, 'blue')
    plt.ylabel('$L^2$-error on $u$'); plt.xlabel('t')
    plt.subplot(312) # L2-norm of the error evolution
    plt.cla(); plt.axis([0.,T,0,0.6])
    plt.plot(TT, store_UH1, 'blue')
    plt.ylabel('$H^1$-error on $u$'); plt.xlabel('t')
    plt.subplot(313) # L2-norm of the error evolution
    plt.cla(); plt.axis([0.,T,0,0.6])
    plt.plot(TT, store_VL2, 'blue')
    plt.ylabel('$L^2$-error on $v$'); plt.xlabel('t')
    
    plt.show()
