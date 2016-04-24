#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Python GetFEM++ interface
#
# Copyright (C) 2016-2016 Yves Renard, Farshid Dabaghi.
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
"""  This program performs a convergence test for several implementations of
     small strain isotropic plasticity in GetFEM++

  $Id: conv_test_small_strain_plasticity 5189 2015-12-15 10:24:07Z renard $
"""

import getfem as gf
import numpy as np
import time
import matplotlib.pyplot as plt
import os
import sys

NT = 128; NX = 256; option = 3; Hi = 0; Hk = 0; load_type = 1; theta = 0.5;
LX=100.; order = 2;
resultspath = './exported_solutions'


def call_test_plasticity():
    return os.system("python test_small_strain_plasticity.py option"+
        ("=%d NX=%d NT=%d Hi=%g Hk=%g theta=%g resultspath=%s load_type=%d order=%d"
         % (option, NX, NT, Hi, Hk, theta, resultspath, load_type, order)))
  


#
# Convergence tests for perfect plasticity, options 1, 2, 3, 4, 5, P1 and P2
#


# Computation of the reference solution if necessary
refname_U  = resultspath+'/ref_perfect_plasticity_U.dat'
refname_mf = resultspath+'/ref_perfect_plasticity_mf.mf'
NT = 256; NX = 256; option = 4; Hi = 0; Hk = 0; load_type = 1;
if (not(os.path.exists(refname_U)) or not(os.path.isfile(refname_U))):
  if (call_test_plasticity() != 0):
      print ('Error in the computation of the reference solution'); exit(1)
  filename = resultspath+('/mf_%d.mf' % (NT))
  os.system("cp %s %s" % (filename, refname_mf))
  filename = resultspath+('/U_%d.dat' % (NT))
  os.system("cp %s %s" % (filename, refname_U))

hrange=[LX/16, LX/22.5, LX/32, LX/45, LX/64, LX/90., LX/128]

for order in [1, 2]:

    theta = 1.
    if (order == 2): theta = 0.5
    errors1 = np.zeros((6,len(hrange),2))-1.
    resname = resultspath+('/conv_test%d.dat' % order)
    if (os.path.exists(resname) and os.path.isfile(resname)):
        errors1 = (np.loadtxt(resname)).reshape(6,len(hrange),2)

    for option in [1, 2, 3, 4, 5]:
        print 'Experiments for option %d and order %d' % (option, order)

        for i in range(0, len(hrange)):
            NX = round(LX/hrange[i]); NT = NX;
            if (errors1[option, i, 1] < 0.):
                if (call_test_plasticity() != 0):
                    print 'Not converged solution'
                    err_L2 = 400
                    err_H1 = 400
                    exit(1)
                else:
                    # Load the final step
                    filename = resultspath+('/mf_%d.mf' % (NT))
                    m = gf.Mesh('load', filename)
                    mf_u = gf.MeshFem('load', filename, m)
                    filename = resultspath+('/U_%d.dat' % (NT))
                    U = np.loadtxt(filename)
            
                    # Load the reference solution
                    m_ref = gf.Mesh('load', refname_mf)
                    mf_u_ref = gf.MeshFem('load', refname_mf, m_ref)
                    U_ref = np.loadtxt(refname_U)
                    mim_ref=gf.MeshIm(m_ref, 6)
                
                    # Estimate of the difference in L2 and H1 norms
                    Ui = gf.compute_interpolate_on(mf_u, U, mf_u_ref)
                    norm_L2 = gf.compute_L2_norm(mf_u_ref, U_ref, mim_ref);
                    err_L2 = gf.compute_L2_dist(mf_u_ref, Ui, mim_ref,
                                                mf_u_ref, U_ref)
                    norm_H1 = gf.compute_H1_semi_norm(mf_u_ref, U_ref, mim_ref);
                    norm_H1 = np.sqrt(pow(norm_L2, 2)+ pow(norm_H1, 2))
                    err_H1 = gf.compute_H1_semi_dist(mf_u_ref, Ui, mim_ref,
                                                     mf_u_ref, U_ref)

                print 'Error in L2 norm: %g' % err_L2
                print 'Error in H1 semi-norm: %g' % err_H1
                errors1[option, i, 0] = err_L2 / norm_L2
                errors1[option, i, 1] = np.sqrt(pow(err_H1,2)
                                                +pow(err_L2,2))/norm_H1

                # Store the result
                np.savetxt(resname, errors1.reshape(errors1.size))

    # Draw the convergence experiment
    fig1, ax1 = plt.subplots()
    plt.show(block=False)
    plt.rc('text', usetex=True)
    ax1.clear()
    l1, = ax1.loglog(hrange, 100.*errors1[1,:,0],linewidth=2,label='option = 1')
    l2, = ax1.loglog(hrange, 100.*errors1[2,:,0],linewidth=2,label='option = 2')
    l3, = ax1.loglog(hrange, 100.*errors1[3,:,0],linewidth=2,label='option = 3')
    l4, = ax1.loglog(hrange, 100.*errors1[4,:,0],linewidth=2,label='option = 4')
    l5, = ax1.loglog(hrange, 100.*errors1[5,:,0],linewidth=2,label='option = 5')
    ax1.legend(loc='best')
    ax1.set_xlabel('h')
    ax1.set_ylabel('$L_2(\Omega)$ error in %')
    ax1.set_title('$L_2(\Omega)$ error for perfect plasticity in $P_%d$'
                  % order)
    # ax1.axis([-0.2, 0.35, -15000, 25000])
    plt.draw()

    fig2, ax2 = plt.subplots()
    plt.show(block=False)
    ax2.clear()
    ax2.loglog(hrange, 100.*errors1[1,:,1], linewidth=2, label='option = 1')
    ax2.loglog(hrange, 100.*errors1[2,:,1], linewidth=2, label='option = 2')
    ax2.loglog(hrange, 100.*errors1[3,:,1], linewidth=2, label='option = 3')
    ax2.loglog(hrange, 100.*errors1[4,:,1], linewidth=2, label='option = 4')
    ax2.loglog(hrange, 100.*errors1[5,:,1], linewidth=2, label='option = 5')
    ax2.legend(loc='best')
    ax2.set_xlabel('h')
    ax2.set_ylabel('$H_1(\Omega)$ error in %')
    ax2.set_title('$H_1(\Omega)$ error for perfect plasticity in $P_%d$'
                  % order)
    # ax2.axis([-0.2, 0.35, -15000, 25000])
    plt.draw()


#
# Convergence tests for plasticity with hardening, options 3, 4, P2
#


# Computation of the reference solution if necessary
refname_U  = resultspath+'/ref_hardening_plasticity_U.dat'
refname_mf = resultspath+'/ref_hardening_plasticity_mf.mf'
NT = 256; NX = 256; option = 4; Hi = 12000; Hk = 12000; load_type = 2;
if (not(os.path.exists(refname_U)) or not(os.path.isfile(refname_U))):
  if (call_test_plasticity() != 0):
      print ('Error in the computation of the reference solution'); exit(1)
  filename = resultspath+('/mf_%d.mf' % (NT))
  os.system("cp %s %s" % (filename, refname_mf))
  filename = resultspath+('/U_%d.dat' % (NT))
  os.system("cp %s %s" % (filename, refname_U))

hrange=[LX/16, LX/22.5, LX/32, LX/45, LX/64, LX/90., LX/128]

for order in [1, 2]:

    theta = 1.
    if (order == 2): theta = 0.5
    errors1 = np.zeros((2,len(hrange),2))-1.
    resname = resultspath+('/conv_test%d.dat' % (order+2))
    if (os.path.exists(resname) and os.path.isfile(resname)):
        errors1 = (np.loadtxt(resname)).reshape(2,len(hrange),2)

    for option in [3, 4]:
        print 'Experiments for option %d and order %d' % (option, order)

        for i in range(0, len(hrange)):
            NX = round(LX/hrange[i]); NT = NX;
            if (errors1[option-3, i, 1] < 0.):
                if (call_test_plasticity() != 0):
                    print 'Not converged solution'
                    err_L2 = 400
                    err_H1 = 400
                    exit(1)
                else:
                    # Load the final step
                    filename = resultspath+('/mf_%d.mf' % (NT))
                    m = gf.Mesh('load', filename)
                    mf_u = gf.MeshFem('load', filename, m)
                    filename = resultspath+('/U_%d.dat' % (NT))
                    U = np.loadtxt(filename)
            
                    # Load the reference solution
                    m_ref = gf.Mesh('load', refname_mf)
                    mf_u_ref = gf.MeshFem('load', refname_mf, m_ref)
                    U_ref = np.loadtxt(refname_U)
                    mim_ref=gf.MeshIm(m_ref, 6)
                
                    # Estimate of the difference in L2 and H1 norms
                    Ui = gf.compute_interpolate_on(mf_u, U, mf_u_ref)
                    norm_L2 = gf.compute_L2_norm(mf_u_ref, U_ref, mim_ref);
                    err_L2 = gf.compute_L2_dist(mf_u_ref, Ui, mim_ref,
                                                mf_u_ref, U_ref)
                    norm_H1 = gf.compute_H1_semi_norm(mf_u_ref, U_ref, mim_ref);
                    norm_H1 = np.sqrt(pow(norm_L2, 2)+ pow(norm_H1, 2))
                    err_H1 = gf.compute_H1_semi_dist(mf_u_ref, Ui, mim_ref,
                                                     mf_u_ref, U_ref)

                print 'Error in L2 norm: %g' % err_L2
                print 'Error in H1 semi-norm: %g' % err_H1
                errors1[option-3, i, 0] = err_L2 / norm_L2
                errors1[option-3, i, 1] = np.sqrt(pow(err_H1,2)
                                                +pow(err_L2,2))/norm_H1

                # Store the result
                np.savetxt(resname, errors1.reshape(errors1.size))

    # Draw the convergence experiment
    fig1, ax1 = plt.subplots()
    plt.show(block=False)
    plt.rc('text', usetex=True)
    ax1.clear()
    l1, = ax1.loglog(hrange, 100.*errors1[0,:,0],linewidth=2,label='option = 3')
    l2, = ax1.loglog(hrange, 100.*errors1[1,:,0],linewidth=2,label='option = 4')
    ax1.legend(loc='best')
    ax1.set_xlabel('h')
    ax1.set_ylabel('$L_2(\Omega)$ error in %')
    ax1.set_title('$L_2(\Omega)$ error for plasticity with hardening in $P_%d$'
                  % order)
    # ax1.axis([-0.2, 0.35, -15000, 25000])
    plt.draw()

    fig2, ax2 = plt.subplots()
    plt.show(block=False)
    ax2.clear()
    ax2.loglog(hrange, 100.*errors1[0,:,1], linewidth=2, label='option = 3')
    ax2.loglog(hrange, 100.*errors1[1,:,1], linewidth=2, label='option = 4')
    ax2.legend(loc='best')
    ax2.set_xlabel('h')
    ax2.set_ylabel('$H_1(\Omega)$ error in %')
    ax2.set_title('$H_1(\Omega)$ error for plasticity with hardening in $P_%d$'
                  % order)
    # ax2.axis([-0.2, 0.35, -15000, 25000])
    plt.draw()



plt.show()
