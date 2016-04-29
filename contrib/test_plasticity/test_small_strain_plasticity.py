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
"""  This program performs a test for several implementations of
     small strain isotropic plasticity in GetFEM++

  $Id: test_small_strain_plasticity.py 5189 2015-12-15 10:24:07Z renard $
"""

import getfem as gf
import numpy as np
import time
import matplotlib.pyplot as plt
import os
import sys

option = 2     # 1 : without hardening, im_data and plastic multiplier
               # 2 : without hardening, with im_data and plastic multiplier
               # 3 : with kinematic and isotropic hardening,
               #     with plastic multiplier
               # 4 : with kinematic and isotropic hardening, with im_data,
               #     without plastic multiplier
               # 5 : with kinematic hardening, with im_data, without plastic
               #     multiplier
               # 6 : Souza-Auricchio model with plastic multiplier (not working)
            
load_type = 1  # 1 : vertical
               # 2 : horizontal
               
bi_material = False
test_tangent_matrix = False
do_plot = True

resultspath = './exported_solutions'




# Physical parameters
lambda_top = 84605.      # Iron
mu_top = 77839.          # Iron
lambda_bottom = 121150.  # Steel
mu_bottom = 80769.       # Steel
von_mises_threshold_top = 8000.
von_mises_threshold_bottom = 7000.

Hk = mu_top/5.; Hi = Hk; Hi = 0. # Kinematic and isotropic hardening parameters


# Numerical parameters
T = 10.
NT = 50
LX = 40.
LY = 20.
NX = 40
theta = 1.; # Parameter for the generalized mid point scheme.
order = 2;

# Arguments from the command line if any
for i in range(1,len(sys.argv)):
  a = sys.argv[i]
  if (a[0:7] == 'option='):
      option = int(a[7:]); print 'option set to %d from argv' % option; continue
  if (a[0:10] == 'load_type='):
      load_type = int(a[10:]); print 'load_type set to %d from argv' % load_type
      continue
  if (a[0:3] == 'NX='):
      NX = int(a[3:]); print 'NX set to %d from argv' % NX; continue
  if (a[0:3] == 'NT='):
      NT = int(a[3:]); print 'NT set to %d from argv' % NT; continue
  if (a[0:6] == 'order='):
      order = int(a[6:]); print 'order set to %d from argv' % order; continue
  if (a[0:3] == 'Hk='):
      Hk = float(a[3:]); print 'Hk set to %g from argv' % Hk; continue
  if (a[0:3] == 'Hi='):
      Hi = float(a[3:]); print 'Hi set to %g from argv' % Hi; continue
  if (a[0:6] == 'theta='):
      theta = float(a[6:]); print 'theta set to %g from argv' % theta; continue
  if (a[0:12] == 'resultspath='):
      resultspath=a[12:]; print 'resultspath set to %s from argv' % resultspath
      continue
  print "Unknow argument '%s' from the command line, exiting" % a; exit(1);

NY = int(np.ceil(NX * LY / (2 * LX))*2)
DT = T/NT


if (do_plot):
    if (not os.path.exists(resultspath)):
        os.makedirs(resultspath)
    print('You can vizualize the optimization steps by launching')
    print('mayavi2 -d %s/von_mises_1.vtk -f WarpVector -m Surface' %resultspath)
    print('mayavi2 -d %s/plast_1.vtk -f WarpVector -m Surface' % resultspath)

if (load_type == 2 and option < 3):
    print('Will not work with this load : will break the body')
    exit(1)

if (load_type == 1):
    f = np.array([0., -400.]);
    t = np.sin(2*np.pi*(np.arange(0, T+DT/2, DT))/20.);
else:
    f = np.array([15000., 0.]);
    t = np.sin(2*np.pi*(np.arange(0, T+DT/2, DT))/10.) + 0.1;


# Create the mesh
# m = gfMesh('triangles grid', [0:(LX/NX):LX], [0:(LY/NY):LY])
m = gf.Mesh('import','structured',
            'GT="GT_PK(2,1)";SIZES=[%d,%d];NOISED=0;NSUBDIV=[%d,%d];'
            % (LX, LY, NX, NY))
N = m.dim()

# Define used MeshIm
mim=gf.MeshIm(m, gf.Integ('IM_TRIANGLE(6)'))

# Define used MeshFem
mf_u=gf.MeshFem(m,2)
mf_u.set_fem(gf.Fem('FEM_PK(2,%d)' % order))

if (option == 1):
    mf_sigma=gf.MeshFem(m,2,2)
    mf_sigma.set_fem(gf.Fem('FEM_PK_DISCONTINUOUS(2,1)'))

# mf_xi = gf.MeshFem(m); mf_xi.set_fem(gf.Fem('FEM_PK(2,2)'));
mf_xi = gf.MeshFem(m); mf_xi.set_fem(gf.Fem('FEM_PK_DISCONTINUOUS(2,1)'))
mf_data=gf.MeshFem(m); mf_data.set_fem(gf.Fem('FEM_PK_DISCONTINUOUS(2,0)'))
mf_vm = gf.MeshFem(m); mf_vm.set_fem(gf.Fem('FEM_PK_DISCONTINUOUS(2,1)'))

# Find the boundary of the domain
P=m.pts()
pidleft=np.compress((abs(P[0,:])<1e-6), range(0, m.nbpts()))
pidright=np.compress((abs(P[0,:] - LX)<1e-6), range(0, m.nbpts()))
fleft  = m.faces_from_pid(pidleft)
fright = m.faces_from_pid(pidright)
m.set_region(1, fleft)
m.set_region(2, fright)


# Decompose the mesh into two regions with different Lame coefficients
if (bi_material):
    separation = LY/2.
else:
    separation = 0

pidtop     = np.compress((abs(P[1,:])>=separation-1E-6), range(0, m.nbpts()))
pidbottom  = np.compress((abs(P[1,:])<=separation+1E-6), range(0, m.nbpts()))
cvidtop    = m.cvid_from_pid(pidtop)
cvidbottom = m.cvid_from_pid(pidbottom)
CVtop      = (mf_data.basic_dof_from_cvid(cvidtop))[0]
CVbottom   = (mf_data.basic_dof_from_cvid(cvidbottom))[0]


# Definition of Lame coeff
cmu = np.zeros(mf_data.nbdof()); clambda = cmu.copy()
von_mises_threshold = cmu.copy()

for i in CVbottom:
    clambda[i] = lambda_bottom
    cmu[i] = mu_bottom
    von_mises_threshold[i] = von_mises_threshold_bottom

for i in CVtop:
    clambda[i] = lambda_top
    cmu[i] = mu_top
    von_mises_threshold[i] = von_mises_threshold_top

# Create the model
md = gf.Model('real')
md.add_fem_variable('u', mf_u)
md.add_fem_data('Previous_u', mf_u);
md.add_initialized_fem_data('lambda', mf_data, clambda)
md.add_initialized_fem_data('mu', mf_data, cmu)
if (option == 1):
    md.add_initialized_fem_data('von_mises_threshold', mf_data, von_mises_threshold*np.sqrt(2./3.))
else:
    md.add_initialized_fem_data('von_mises_threshold', mf_data, von_mises_threshold)
md.add_Dirichlet_condition_with_multipliers(mim, 'u', mf_u, 1)
md.add_initialized_fem_data('VolumicData', mf_data, np.zeros(N*mf_data.nbdof()))
md.add_source_term_brick(mim, 'u', 'VolumicData', 2)
mim_data = gf.MeshImData(mim, -1, [N, N])
mim_data_scal = gf.MeshImData(mim, -1, 1)

if (option == 1):
    md.add_fem_data('sigma', mf_sigma);
    md.add_elastoplasticity_brick(mim, 'VM', 'u', 'Previous_u', 'lambda', 'mu',
                                  'von_mises_threshold', 'sigma');

if (option == 2):
    md.add_fem_variable('xi', mf_xi)
    md.add_initialized_data('theta', [theta])
    md.add_initialized_data('r', [1e-8])
    md.add_im_data('Epn', mim_data)

    # md.add_small_strain_elastoplasticity_brick(mim, 'u', 'Epn', 'lambda', 'mu', 'von_mises_threshold', 'theta');
    # exit(1);
    
    
    if (theta == 1):
        Etheta = '(Sym(Grad_u))';
        Eptheta = '((Epn+2*mu*xi*Deviator('+Etheta+'))/(1+2*mu*xi))'
        Epnp1 = Eptheta
        sigma_np1 = ('(lambda*Trace(Sym(Grad_u)-'+Epnp1
                     +')*Id(meshdim) + 2*mu*(Sym(Grad_u)-'+Epnp1+'))')
        sigma_theta = sigma_np1;
    else:
        Etheta = '(Sym(theta*Grad_u+(1-theta)*Grad_Previous_u))'
        Eptheta = '((Epn+2*mu*theta*xi*Deviator('+Etheta+'))/(1+2*mu*theta*xi))'
        Epnp1 = '(('+Eptheta+'-(1-theta)*Epn)/theta)'
        sigma_np1 = ('(lambda*Trace(Sym(Grad_u)-'+Epnp1
                     +')*Id(meshdim) + 2*mu*(Sym(Grad_u)-'+Epnp1+'))')
        sigma_theta = ('(lambda*Trace('+Etheta+'-'+Eptheta
                       +')*Id(meshdim) + 2*mu*('+Etheta+'-'+Eptheta+'))')
    
    fbound = '(Norm(Deviator('+sigma_theta+'))-sqrt(2/3)*von_mises_threshold)'
    # fbound = '(Norm('+Eptheta+'-Epn)-theta*xi*von_mises_threshold)'
    expr = sigma_np1+':Grad_Test_u+(1/r)*(xi-pos_part(xi+r*'+fbound+'))*Test_xi'
    # expr = sigma_np1+':Grad_Test_u+('+fbound+'+pos_part(-xi/r-'+fbound+
    #        '))*Test_xi'
    md.add_nonlinear_generic_assembly_brick(mim, expr)
        
if (option == 3):
    md.add_fem_variable('xi', mf_xi)
    md.add_initialized_data('theta', [theta])
    md.add_initialized_data('r', [1e-8])
    md.add_im_data('Epn', mim_data)
    md.add_im_data('alphan', mim_data_scal)
    md.add_initialized_data('Hk', [Hk])
    md.add_initialized_data('Hi', [Hi])
    
    if (theta == 1):
        Etheta = '(Sym(Grad_u))'
        Eptheta = '((Epn+2*mu*xi*Deviator('+Etheta+'))/(1+(2*mu+Hk)*xi))'
        Epnp1 = Eptheta
        sigma_np1 = ('(lambda*Trace(Sym(Grad_u)-'+Epnp1
                     +')*Id(meshdim) + 2*mu*(Sym(Grad_u)-'+Epnp1+'))')
        sigma_theta = sigma_np1
        alpha_theta = ('(alphan+xi*(Norm(2*mu*Deviator('+Etheta
                       +')-(2*mu+Hk)*'+Eptheta+')))')
        alpha_np1 = alpha_theta
    else:
        Etheta = '(Sym(theta*Grad_u+(1-theta)*Grad_Previous_u))'
        Eptheta = '((Epn+2*mu*theta*xi*Deviator('+Etheta+'))/(1+(2*mu+Hk)*theta*xi))'
        Epnp1 = '(('+Eptheta+' - (1-theta)*Epn)/theta)'
        sigma_np1 = ('(lambda*Trace(Sym(Grad_u)-'+Epnp1
                     +')*Id(meshdim) + 2*mu*(Sym(Grad_u)-'+Epnp1+'))')
        sigma_theta = ('(lambda*Trace('+Etheta+'-'+Eptheta
                       +')*Id(meshdim) + 2*mu*('+Etheta+'-'+Eptheta+'))')
        alpha_theta = '(alphan+sqrt(2/3)*theta*xi*(Norm(2*mu*Deviator('+Etheta+')-(2*mu+Hk)*'+Eptheta+')))'
        # alpha_np1 = '(('+alpha_theta+' - (1-theta)*alphan)/theta)'
        alpha_np1 = ('(alphan+sqrt(2/3)*xi*(Norm(2*mu*Deviator('+Etheta
                     +')-(2*mu+Hk)*'+Eptheta+')))')
        # alpha_theta = '(alphan+sqrt(2/3)*Norm('+Eptheta+'-Epn))'  # do not work
        # alpha_np1 = '(alphan+sqrt(2/3)*Norm('+Eptheta+'-Epn)/theta)' # do not work
    
    # fbound = ('(Norm(Deviator('+sigma_theta+')-Hk*'+Eptheta
    #           +') - von_mises_threshold - Hi*'+alpha_theta+')')
    fbound = ('(Norm(2*mu*Deviator('+Etheta+')-(2*mu+Hk)*'+Eptheta
              +') - sqrt(2/3)*(von_mises_threshold + Hi*'+alpha_theta+'))')
    expr = (sigma_np1+':Grad_Test_u + (1/r)*(xi - pos_part(xi+r*'+fbound
            +'))*Test_xi')
    # expr = (sigma_np1+':Grad_Test_u + ('+fbound+' + pos_part(-xi/r-'+fbound
    #         +'))*Test_xi')
    md.add_nonlinear_generic_assembly_brick(mim, expr)
    
if (option == 4):
    md.add_initialized_data('theta', [theta])
    md.add_im_data('Epn', mim_data)
    md.add_im_data('alphan', mim_data_scal)
    md.add_initialized_data('Hk', [Hk])
    md.add_initialized_data('Hi', [Hi])
    
    Etheta = '(Sym(theta*Grad_u+(1-theta)*Grad_Previous_u))'
    Btheta = '((2*mu)*Deviator('+Etheta+')-(2*mu+Hk)*Epn)'
    alpha_theta = ('(max(alphan, (sqrt(3/2)*(2*mu+Hk)*alphan+Norm('+Btheta+
                   ') - sqrt(2/3)*von_mises_threshold)/(sqrt(3/2)*(2*mu+Hk)'+
                   '+sqrt(2/3)*Hi)))')
    alpha_np1 = '(('+alpha_theta+' - (1-theta)*alphan)/theta)'
    Eptheta= ('(Epn+(1/(2*mu+Hk))*pos_part(1-sqrt(2/3)*(von_mises_threshold+Hi*'
              +alpha_theta+')/(Norm('+Btheta+')+1e-25))*'+Btheta+')')
    Epnp1 = '(('+Eptheta+' - (1-theta)*Epn)/theta)'
    sigma_np1 = ('(lambda*Trace(Sym(Grad_u)-'+Epnp1
                 +')*Id(meshdim) + 2*mu*(Sym(Grad_u)-'+Epnp1+'))')
    sigma_theta = ('(lambda*Trace('+Etheta+'-'+Eptheta
                   +')*Id(meshdim) + 2*mu*('+Etheta+'-'+Eptheta+'))')
    
    expr = sigma_np1+':Grad_Test_u'
    md.add_nonlinear_generic_assembly_brick(mim, expr)
    
      
if (option == 5):
    md.add_im_data('sigma', mim_data)
    md.add_initialized_data('theta', [theta])
    md.add_initialized_data('Hk', [Hk])

    Is = 'Reshape(Id(meshdim*meshdim),meshdim,meshdim,meshdim,meshdim)'
    IxI = 'Id(meshdim)@Id(meshdim)'
    coeff_long ='((lambda)*(Hk))/((2*(mu)+(Hk))*(meshdim*(lambda)+2*(mu)+(Hk)))'
    B_inv = '((2*(mu)/(2*(mu)+(Hk)))*('+Is+')+('+coeff_long+')*('+IxI+'))'
    B = ('((1+(Hk)/(2*(mu)))*('+Is
         +')-(((lambda)*(Hk))/(2*(mu)*(meshdim*(lambda)+2*(mu))))*('+IxI+'))')
    ApH = '((2*(mu)+(Hk))*('+Is+') + (lambda)*('+IxI+'))'
    Enp1 = 'Sym(Grad_u)'
    En = 'Sym(Grad_Previous_u)'
  
    # Expression de sigma for generalized mid-point scheme
    expr_sigma = ('('+B_inv+'*(Von_Mises_projection(('+B
                  +'*((1-theta)*sigma))+(-(Hk)*(((1-theta)*'+En+')+(theta*'+Enp1
                  +')))+(theta*'+ApH+'*('+Enp1+'-'+En+')) + (theta*'+B
                  +'*sigma), sqrt(2/3)*von_mises_threshold)+(Hk)*(((1-theta)*'
                  +En+')+(theta*'+Enp1+'))))')
  
    md.add_nonlinear_generic_assembly_brick(mim, expr_sigma+':Grad_Test_u')
    
if (option == 6):
    md.add_fem_variable('xi', mf_xi)
    md.add_fem_variable('delta', mf_xi)
    md.add_initialized_data('theta', [theta])
    md.add_initialized_data('r1', [1e-8])
    md.add_initialized_data('r2', [1])
    md.add_im_data('Epn', mim_data)
    md.add_initialized_data('c1', [0])
    md.add_initialized_data('c2', [Hk])
    md.add_initialized_data('c3', [0.3])
    
    # Simplified model without c1 and for theta = 1 and delta as
    # additional multiplier (not working ...)
    
    Etheta = '(Sym(theta*Grad_u+(1-theta)*Grad_Previous_u))'
    Btheta = '(Epn+theta*xi*2*mu*Deviator('+Etheta+'))'
    Eptheta = '('+Btheta+'/(1+(2*mu+c2+delta)*theta*xi))'
    # Eptheta = '('+Btheta+'*min(c3/(Norm(',Btheta,')+1e-5),'
    # +'1/(1+2*mu*theta*xi)))')
    Epnp1 = '(('+Eptheta+' - (1-theta)*Epn)/theta)'
    sigma_np1 = ('(lambda*Trace(Sym(Grad_u)-'+Epnp1
                 +')*Id(meshdim) + 2*mu*(Sym(Grad_u)-'+Epnp1+'))')
    sigma_theta = ('(lambda*Trace('+Etheta+'-'+Eptheta
                   +')*Id(meshdim) + 2*mu*('+Etheta+'-'+Eptheta+'))')
    
    # fbound = ('(Norm(2*mu*Deviator('+Etheta+')-(2*mu+pos_part((Norm('+Btheta
    #           +')-c3)/(c3*theta*xi+1e-25)-2*mu))*'+Eptheta
    #           +')-sqrt(2/3)*von_mises_threshold)')
    fbound = ('(Norm(2*mu*Deviator('+Etheta+')-(2*mu+c2+delta)*'+Eptheta
              +') - sqrt(2/3)*von_mises_threshold)')
    fbound_delta = '(Norm('+Eptheta+')-c3)'
    # fbound = ('(Norm(2*mu*Deviator('+Etheta+')-(2*mu+c2)*'+Eptheta
    #           +') - sqrt(2/3)*von_mises_threshold)')
    expr = (sigma_np1+':Grad_Test_u + (1/r1)*(xi - pos_part(xi+r1*'+fbound
            +'))*Test_xi-(1/r2)*(delta-pos_part(delta+r2*'+fbound_delta
            +'))*Test_delta')
    md.add_nonlinear_generic_assembly_brick(mim, expr)
    
    if (False):
        Etheta = '(Sym(theta*Grad_u+(1-theta)*Grad_Previous_u))'
        Btheta = '(Epn+theta*xi*2*mu*Deviator('+Etheta+'))'
        # Eptheta = ('(Print('+Btheta+')*min(c3/(max(Norm('+Btheta
        #            +'), c3/2)), pos_part(1-theta*xi*c1/(Norm('+Btheta
        #            +')+0.001))/(1+(2*mu+c2)*theta*xi)))')
        # Eptheta = ('('+Btheta+'*min(c3/(max(Norm('+Btheta
        #            +'), c3/2)), 1/(1+(2*mu+c2)*theta*xi)))')
        Eptheta = ('('+Btheta+'*min(c3/(Norm('+Btheta
                   +')+1e-10), 1/(1+(2*mu+c2)*theta*xi)))')
        Epnp1 = '(('+Eptheta+' - (1-theta)*Epn)/theta)'
        sigma_np1 = ('(lambda*Trace(Sym(Grad_u)-'+Epnp1
                     +')*Id(meshdim) + 2*mu*(Sym(Grad_u)-'+Epnp1+'))')
        sigma_theta = ('(lambda*Trace('+Etheta+'-'+Eptheta
                       +')*Id(meshdim) + 2*mu*('+Etheta+'-'+Eptheta+'))')
      
        # fbound = ('(Norm(2*mu*Deviator('+Etheta+')-(2*mu+c2)*'+Eptheta
        #           +'-max(c1, (Norm('+Btheta+')-c3)/(theta*xi+1e-25)'
        #           +'-(2*mu+c2)*c3)*Normalized('+Eptheta
        #           +'))-sqrt(2/3)*von_mises_threshold)')
        # fbound = ('(Norm(2*mu*Deviator('+Etheta+')-(2*mu+c2)*'+Eptheta
        #           +'-pos_part(pos_part(Norm('+Btheta
        #           +')/(theta*xi+1e-10)-c1)/c3-(1/(theta*xi+1e-10)+2*mu+c2))*'
        #           +Eptheta+') - sqrt(2/3)*von_mises_threshold)')
        fbound = ('(Norm(2*mu*Deviator('+Etheta+')-(2*mu+c2)*'+Eptheta
                  +'-pos_part( Norm('+Btheta
                  +')/(c3*(theta*xi+1e-10)) - (1/(theta*xi+1e-10)+2*mu+c2))*'
                  +Eptheta+') - sqrt(2/3)*von_mises_threshold)')
        # fbound = ('(Norm(2*mu*Deviator('+Etheta+')-(2*mu+c2)*'+Eptheta
        #           +') - sqrt(2/3)*von_mises_threshold)')
        expr = (sigma_np1+':Grad_Test_u + (1/r)*(xi - pos_part(xi+r*'+fbound
                +'))*Test_xi')
        # expr = (sigma_np1+':Grad_Test_u + ('+fbound+' + pos_part(-xi/r-'
        #         +fbound+'))*Test_xi')
        md.add_nonlinear_generic_assembly_brick(mim, expr)
      



VM=np.zeros(mf_vm.nbdof())
sigma_fig = np.zeros(len(t))
Epsilon_u_fig = np.zeros(len(t))


for step in range(0, len(t)):

    print('step %d / %d, coeff = %g' % (step, len(t)-1, t[step]))
    source = mf_data.eval('[f[0]*t[step], f[1]*t[step]]', globals(), locals())
    md.set_variable('VolumicData', source)

    if (test_tangent_matrix):
        md.test_tangent_matrix(1E-8, 10, 0.000001)
   
    # Solve the system
    md.solve('noisy', 'lsearch', 'simplest',  'alpha min', 0.4, 'max_iter', 50,
             'max_res', 1e-6)
    # md.solve('noisy', 'max_iter', 80)

    U = md.variable('u')

    # Compute Von Mises and plastic part for graphical post-treatment
    if (do_plot):
        if (option == 1):
            sigma_np1 = 'sigma';
        else:
            if (option == 5):
                sigma_np1 = expr_sigma;
                coeff1='-lambda/(2*mu*(meshdim*lambda+2*mu))'
                coeff2='1/(2*mu)'
                Ainv = ('(%s)*(%s) + (%s)*(%s)' % (coeff1, IxI, coeff2, Is))
                Ep = ('Norm(Sym(Grad_u) - (%s)*sigma)' % Ainv)
            else:
                Ep ='Norm('+Epnp1+')'
        
            von_mises = 'sqrt(3/2)*Norm(Deviator('+sigma_np1+'))'
            VM = md.local_projection(mim, von_mises, mf_vm);
            plast = md.local_projection(mim, Ep, mf_vm);
    
        sigma = md.interpolation(sigma_np1, mim_data);
        Epsilon_u = md.interpolation('Sym(Grad_u)', mim_data);
        ind_gauss_pt = 22500
        if (sigma.shape[0] <= N*N*(ind_gauss_pt + 1)):
            ind_gauss_pt = int(3.*sigma.shape[0] / (4*N*N))
            
        sigma_fig[step]=sigma[N*N*ind_gauss_pt]
        Epsilon_u_fig[step]=Epsilon_u[N*N*ind_gauss_pt]
    
    # Compute new plastic internal variables
    if (option == 1):
        md.elastoplasticity_next_iter(mim, 'u', 'Previous_u', 'VM', 'lambda',
                                      'mu', 'von_mises_threshold', 'sigma')
        plast = md.compute_plastic_part(mim, mf_vm, 'u', 'Previous_u', 'VM',
                                        'lambda', 'mu', 'von_mises_threshold',
                                        'sigma')
        # Compute Von Mises or Tresca stress
        VM = md.compute_elastoplasticity_Von_Mises_or_Tresca('sigma', mf_vm,
                                                             'Von Mises')

    if (option == 2):
        NewEpn = md.interpolation(Epnp1, mim_data)
        md.set_variable('Epn', NewEpn)
        md.set_variable('Previous_u', U)
        
    if (option == 3 or option == 4):
        NewEpn = md.interpolation(Epnp1, mim_data)
        Newalphan = md.interpolation(alpha_np1, mim_data_scal)
        md.set_variable('Epn', NewEpn)
        md.set_variable('alphan', Newalphan)
        md.set_variable('Previous_u', U)
        
    if (option == 5):
        sigma_0 = md.variable('sigma')
        sigma = md.interpolation(expr_sigma, mim_data)
        U_0 = md.variable('Previous_u')
        U_ntheta = theta*U + (1-theta)*U_0
      
        sigma = (sigma - (1-theta)*sigma_0)/theta
        md.set_variable('sigma', sigma)
        md.set_variable('Previous_u', U)
        
    if (option == 6):
        NewEpn = md.interpolation(Epnp1, mim_data)
        md.set_variable('Epn', NewEpn)
        md.set_variable('Previous_u', U)
      
       
    if (do_plot):
        # Export Von Mises and plastic part
        filename = resultspath+('/von_mises_%d.vtk' % (step))
        mf_vm.export_to_vtk(filename, mf_vm, VM.reshape((len(VM))),
                            'Von Mises Stresses', mf_u, U, 'Displacements')
        filename = resultspath+('/plast_%d.vtk' % (step))
        mf_vm.export_to_vtk(filename, mf_vm, plast.reshape((len(plast))),
                            'Norm of plastic strain tensor', mf_u, U,
                            'Displacements')
        filename = resultspath+('/mf_%d.mf' % (step))
        mf_u.save(filename, 'with_mesh')
        filename = resultspath+('/U_%d.dat' % (step))
        np.savetxt(filename, U)
        

        # Draw the stress/strain evolution on a selected Gauss point
        if (step == 0):
            fig, ax = plt.subplots()
            plt.show(block=False)
                
        ax.clear()
        ax.plot(Epsilon_u_fig[0:(step+1)], sigma_fig[0:(step+1)], linewidth=2)
        ax.set_xlabel('Strain')
        ax.set_ylabel('Stress')
        ax.set_title('stress/strain evolution')
        ax.axis([-0.2, 0.35, -15000, 25000])
        plt.draw()

 

# time.sleep(1)
# plt.show()
