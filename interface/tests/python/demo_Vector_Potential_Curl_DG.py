#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Python GetFEM interface
#
# Copyright (C) 2019-2020 Egor Vtorushin.
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
"""
Curl based problem in 3d cube test, including DG.
We solve equation for vector potential approach for the magnetostatic problem as it described in(and many other papers)
"Mixed formulations for finite element analysis of magnetostatic and electrostatic problems"
Fumio Kikuchi, Japan J. Appl. Math. (1989) 6: 209.
For a Interior penalty terms understanding we refer to 
"Discontinuous galerkin approximation of the maxwell eigenproblem",
A. Buffa AND I. Perugia, SIAM Journal on Numerical Analysis, Vol. 44, No. 5 : pp. 2198-222, 2006. 
"Unified analysis of discontinuous Galerkin methods for elliptic problems", 
D.N. Arnold, F. Brezzi, B. Cockburn, L.D. Marini, SIAM J. Numer. Anal. vol. 39:5, pp 1749-1779, 2002.
@author: Egor Vtorushin
@mail: vtorushin@gmail.com
"""

import datetime
import sys
import time
from enum import Enum

import numpy as np

try: gf
except NameError: gf = None
if gf is None:
    import getfem as gf
Exact = Enum('Exact', 'BdmExample DotNormalZero CrossNormalZero')

######################################################################
# Parameters 
CUBE_H=1.001
MESHSIZE=0.25
ALPHAF=1.e3/MESHSIZE #Interior penalty(alpha) factor
INNER_FACES= 42
ALL_SIDES=1
VERIFY_NEIGHBOR_COMPUTATION=True
DS_GA = True if INNER_FACES>0 else False

EXACT_TYPE=Exact.BdmExample
# Parameters END
######################################################################

######################################################################
# Exact functions 
V1 = np.pi/4.0 #for a case
if (EXACT_TYPE==Exact.BdmExample):
    Bexact_exp=f'[{V1}*(sin(pi*X(2)) - sin(pi*X(3))); {V1}*(-sin(pi*X(1)) + sin(pi*X(3))); {V1}*(sin(pi*X(1)) - sin(pi*X(2)))]'
    CurlBexact_exp=f'[{V1}*(-pi*cos(pi*X(2)) - pi*cos(pi*X(3)));{V1}*(-pi*cos(pi*X(1)) - pi*cos(pi*X(3))); {V1}*(-pi*cos(pi*X(1)) - pi*cos(pi*X(2)))]'
elif(EXACT_TYPE==Exact.CrossNormalZero):    
    Bexact_exp=f'{V1}*[sin(pi*X(2)) * sin(pi*X(3)); -sin(pi*X(1)) * sin(pi*X(3)); sin(pi*X(1)) * sin(pi*X(2))]'
    CurlBexact_exp=f'{V1}*[pi*sin(pi*X(1))*cos(pi*X(2)) + pi*sin(pi*X(1))*cos(pi*X(3)); -pi*sin(pi*X(2))*cos(pi*X(1)) + pi*sin(pi*X(2))*cos(pi*X(3)); -pi*sin(pi*X(3))*cos(pi*X(1)) - pi*sin(pi*X(3))*cos(pi*X(2))]'
elif(EXACT_TYPE==Exact.DotNormalZero):
    Bexact_exp=f'{V1}*[sin(pi*X(1))*cos(2*pi*X(2)); -sin(pi*X(2))*cos(2*pi*X(3)); sin(pi*X(3))*cos(2*pi*X(1)) ]'
    CurlBexact_exp=f'{V1}*[-2*pi*sin(pi*X(2))*sin(2*pi*X(3)); 2*pi*sin(2*pi*X(1))*sin(pi*X(3)); 2*pi*sin(pi*X(1))*sin(2*pi*X(2))]'
# Exact functions END
######################################################################

######################################################################
# Creation of a simple 3d simplex mesh 
xa = np.arange(0, CUBE_H, MESHSIZE)
ya = np.arange(0, CUBE_H, MESHSIZE)
za = np.arange(0, CUBE_H, MESHSIZE)
m = gf.Mesh('regular simplices', xa, ya, za)
#m.export_to_vtk('3dregular_simplices.vtk', 'ascii')
# Creation of a simple 3d simplex mesh END
######################################################################

######################################################################
# Create a MeshFem for b and la mult fields
DS='_DISCONTINUOUS' if DS_GA else ""
ALPHA=',0.5' if DS_GA else ""

mfb = gf.MeshFem(m,3)
mfb.set_qdim(3)
mfb.set_fem(gf.Fem(f'FEM_PK{DS}(3,1{ALPHA})'))

mfdiv_mult = gf.MeshFem(m,1)
mfdiv_mult.set_qdim(1)
mfdiv_mult.set_fem(gf.Fem(f'FEM_PK{DS}(3,1{ALPHA})')) # MeshFem FEM_PK FEM_PK_DISCONTINUOUS
# Create a MeshFem for b and la mult fields END
######################################################################

# Integration method used
mim = gf.MeshIm(m, gf.Integ('IM_TETRAHEDRON(2)'))

#########################################################################
# Make it boundary 
m.set_region(ALL_SIDES, m.outer_faces())
if(DS_GA):
    m.set_region(INNER_FACES,m.inner_faces())
# Make it boundary END	
###########################################################################
    
VPRODUCT_VAR_N_='Cross_product(Normal,{var})'
JUMP_VPRODUCT_VAR_='Cross_product(Normal,{var})-Cross_product(Normal,Interpolate({var},neighbor_element))'
MEAN_VAR_='(0.5*{M}+0.5*{N})'
JUMP_VVAR_N_ = '({v}-Interpolate({v},neighbor_element)).Normal'
JUMP_SVAR_N_ = '({V}-Interpolate({V},neighbor_element))*Normal'
CURL_VAR_='[Grad_{var}(3,2)-Grad_{var}(2,3); Grad_{var}(1,3)-Grad_{var}(3,1); Grad_{var}(2,1)-Grad_{var}(1,2)]'
CURL_INT_VAR_='[Interpolate(Grad_{var}, neighbor_element)(3,2)-Interpolate(Grad_{var}, neighbor_element)(2,3); Interpolate(Grad_{var}, neighbor_element)(1,3)-Interpolate(Grad_{var}, neighbor_element)(3,1); Interpolate(Grad_{var}, neighbor_element)(2,1)-Interpolate(Grad_{var}, neighbor_element)(1,2)]'

mean_test_p = '((Test_p+Interpolate(Test_p,neighbor_element))*0.5)'
jump_p = JUMP_SVAR_N_.format(V='p')
jump_test_p = JUMP_SVAR_N_.format(V='Test_p')
jump_b = JUMP_VVAR_N_.format(v='b')
jump_test_b = JUMP_VVAR_N_.format(v='Test_b')
cross_jump_b = JUMP_VPRODUCT_VAR_.format(var='b')
cross_jump_test_b = JUMP_VPRODUCT_VAR_.format(var='Test_b')
curl_int_b=CURL_INT_VAR_.format(var='b')
curl_b=CURL_VAR_.format(var='b')
mean_curl_b=MEAN_VAR_.format(M=curl_b,N=curl_int_b)
curl_int_test_b=CURL_INT_VAR_.format(var='Test_b')
curl_test_b=CURL_VAR_.format(var='Test_b')
mean_curl_test_b=MEAN_VAR_.format(M=curl_test_b,N=curl_int_test_b)
curl_testb_curl_b_exp='({L}).({M})'.format(L=CURL_VAR_.format(var='Test_b'), M=CURL_VAR_.format(var= 'b'))
cross_jump_a = JUMP_VPRODUCT_VAR_.format(var='a')
cross_jump_test_a = JUMP_VPRODUCT_VAR_.format(var='Test_a')
jump_a = JUMP_VVAR_N_.format(v='a')
jump_test_a = JUMP_VVAR_N_.format(v='Test_a')

#magentec potential model     
mdmag=gf.Model('real');

######################################################################################################################
#Interpolate and export exact fields 
intonme = mfb
CurlBexact = mdmag.interpolation(CurlBexact_exp, intonme)
intonme.export_to_vtk('curl_exact3d.vtk', 'ascii', CurlBexact)
print('exact curl B norm l2', gf.compute_L2_norm(intonme, CurlBexact, mim))
Bexact=mdmag.interpolation(Bexact_exp, intonme)
intonme.export_to_vtk('field_exact3d.vtk', 'ascii', Bexact)
print('exact B norm l2', gf.compute_L2_norm(intonme, Bexact, mim))
#Interpolate and export exact fields END
######################################################################################################################

WITH_B_CROSS_N=True

# Main potential unknown
mdmag.add_fem_variable('b', mfb);
# lagrange augmentation unknown for potential
mdmag.add_fem_variable('p', mfdiv_mult)

tm = time.time()
# building (curl b, curl \tau) term 
CURL_B=gf.asm('generic', mim, 2, curl_testb_curl_b_exp, -1, 'b', 1, mfb, 0)
tm = time.time()-tm
print("CURL.CURL elapsed time", tm)

# (curl b, curl \tau) adding to the model
magi=mdmag.add_explicit_matrix('b', 'b', CURL_B)
 
#####################################################################################################
# Interior penalty terms
if(DS_GA):
    J_CROSSTAU_DOT_M_CURLB = gf.asm('generic', mim, 2, '({K}).{L}'.format(K=cross_jump_test_b, L=mean_curl_b), INNER_FACES, 'b', 1, mfb, 0)
    #J_CROSSB_DOT_M_CURLTAU=gf.asm('generic', mim, 2, '({K}).{L}'.format(K=cross_jump_b, L=mean_curl_test_b), INNER_FACES, 'b', 1, mfb, 0)
    CROSS_CROSS_IP_EXP = "({alpha})*({L}).({M})".format(alpha=ALPHAF,L=cross_jump_b,M=cross_jump_test_b)
    CROSS_CROSS_IP=gf.asm('generic', mim, 2, CROSS_CROSS_IP_EXP, INNER_FACES, 'b', 1, mfb, 0)
    bi=mdmag.add_explicit_matrix('b', 'b', CROSS_CROSS_IP)
    bi=mdmag.add_explicit_matrix('b', 'b', -J_CROSSTAU_DOT_M_CURLB)
    #bi=mdmag.add_explicit_matrix('b', 'b', -J_CROSSB_DOT_M_CURLTAU)
# Interior penalty terms END
#####################################################################################################

#####################################################################################################
# Lagrange multipliers for div B=0 condition
DIV_B1=gf.asm('generic', mim, 2, 'Grad_Test_p.Test2_b', -1, 'b', 1, mfb, 0, 'p', 1, mfdiv_mult, 0)
DIV_B1 = gf.Spmat('copy', DIV_B1 ,range(mfb.nbdof(), mfb.nbdof()+mfdiv_mult.nbdof()),range(mfb.nbdof()))
magi=mdmag.add_explicit_matrix('p', 'b', DIV_B1)
DIV_B1.transpose()
magi=mdmag.add_explicit_matrix('b', 'p', DIV_B1)
if(DS_GA): # Interior penalty terms for div B=0 part
    FLUX_IP_EXP = "({alpha})*({J}).({K})".format(alpha=ALPHAF,J=jump_b,K=jump_test_b)
    FLUX_IP=gf.asm('generic', mim, 2, FLUX_IP_EXP, INNER_FACES, 'b', 1, mfb, 0)
    bi=mdmag.add_explicit_matrix('b', 'b', FLUX_IP)
    IP_ADD=gf.asm('generic', mim, 2, jump_b+'.'+mean_test_p, INNER_FACES, 'b', 1, mfb, 0, 'p', 1, mfdiv_mult, 0)
    IP_ADD = gf.Spmat('copy', IP_ADD ,range(mfb.nbdof(), mfb.nbdof()+mfdiv_mult.nbdof()),range(mfb.nbdof()))
    magi=mdmag.add_explicit_matrix('p', 'b', -IP_ADD)
    IP_ADD.transpose()
    magi=mdmag.add_explicit_matrix('b', 'p', -IP_ADD)
    PDIV_IP_EXP="({alpha})*({L}).({M})".format(alpha=ALPHAF,L=jump_p,M=jump_test_p)
    PDIV_IP=gf.asm('generic', mim, 2, PDIV_IP_EXP, INNER_FACES, 'p', 1, mfdiv_mult, 0)
    bi=mdmag.add_explicit_matrix('p', 'p', PDIV_IP)
# Lagrange multipliers for div B=0 condition END
#####################################################################################################

#####################################################################################################
# building rhs with given curlB field and B x N BC  
B_RHS=gf.asm('generic', mim, 1, 'vec.Test_b', -1, 'vec', 0, intonme, CurlBexact, 'b', 1, mfb, 0)
print("B_RHS L2_norm", gf.compute_L2_norm(mfb, B_RHS, mim))
mdmag.add_explicit_rhs('b', B_RHS)
if(WITH_B_CROSS_N):
    B_CROSS_N_RHS_EXP='vec.({S})'.format(S=VPRODUCT_VAR_N_.format(var='Test_b')) # compute B x N
    tm = time.time()
    B_CROSS_N_RHS = gf.asm('generic', mim, 1, B_CROSS_N_RHS_EXP, ALL_SIDES, 'vec', 0, mfb, Bexact, 'b', 1, mfb, 0)
    tm = time.time()-tm
    print("B . (N x tau) term elapsed time", tm)
    B_CROSS_N_RHS_NORM=gf.compute_L2_norm(mfb, B_CROSS_N_RHS, mim)
    print("B . (N x tau) L2_norm", B_CROSS_N_RHS_NORM)
    if(B_CROSS_N_RHS_NORM>MESHSIZE**4):
        print('B_CROSS_N_RHS added')
        mdmag.add_explicit_rhs('b', B_CROSS_N_RHS)
# building rhs with given curlB field and B x N BC END		
#####################################################################################################

#####################################################################################################
# Assembly of the linear system for potential and solve.
try:
    mdmag.assembly(option='build_matrix')
    mdmag.assembly(option='build_rhs')
except:
    print("Something went wrong when building magnitic filed matrix")
finally:
    print("vector potential started at", time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime()), 'total dof number', mdmag.rhs().size)
    tm = time.time()
    #mdmag.solve()
    bb=gf.linsolve_superlu(mdmag.tangent_matrix(),mdmag.rhs())[0][:,0]
    mdmag.set_variable('b', bb[:B_RHS.size])
    mdmag.set_variable('p', bb[B_RHS.size:])
    tm = time.time()-tm
    print("vector potential elapsed time %s " % datetime.timedelta(seconds=tm))
    print("vector potential L2_norm", gf.compute_L2_norm(mfb, mdmag.variable('b'), mim))
# Assembly of the linear system for potential and solve. END
#####################################################################################################
mfb.export_to_vtk('vector_potential3d.vtk', 'ascii', mdmag.variable('b')) # esport magnetic potential

####################################################################
# done vector potential part
####################################################################

if(VERIFY_NEIGHBOR_COMPUTATION):    
    m.set_region(42,m.inner_faces())
    print ("jump A x N conv", np.linalg.norm(gf.asm('generic', mim, 1, '({A}).({B})'.format(A=cross_jump_b, B=cross_jump_test_b), 42, mdmag)))
    print ("jump A . N conv", np.linalg.norm(gf.asm('generic', mim, 1, '({A}).({B})'.format(A=jump_b,B=jump_test_b), 42, mdmag)))
    print ("jump p * N conv", np.linalg.norm(gf.asm('generic', mim, 1, '({A}).({B})'.format(A=jump_p,B=jump_test_p), 42, mdmag)))

################################################################################################################
# start restoring target field from the vector potential 
# Technically we solve equation
# find b \in H_{curl}(\Omega) such as
# (b,\tau) = (curl A, \tau), for any \tau\in H_{curl}(\Omega)
# for given field A\in V, where space V can be H_{div} or H^1 
# To do it one should apply green formula for (curl A, \tau) term and carefully account BC term (\tau x N).A
################################################################################################################

# Model for restoration 
mdpot=gf.Model('real');
# Main unknown
mdpot.add_fem_variable('a', mfb);
# mass matrix create and add
AA=gf.asm('mass matrix', mim, mfb, mfb, -1)
poti=mdpot.add_explicit_matrix('a', 'a', AA)
#####################################################################################################
# Interior penalty terms
if(DS_GA):
    CROSSA_CROSSA_IP_EXP = "({alpha})*({J}).({K})".format(alpha=ALPHAF,J=cross_jump_a,K=cross_jump_test_a)
    CROSSA_CROSSA_IP=gf.asm('generic', mim, 2, CROSSA_CROSSA_IP_EXP, INNER_FACES, 'a', 1, mfb, 0)
    poti=mdpot.add_explicit_matrix('a', 'a', CROSSA_CROSSA_IP)
# Interior penalty terms END
#####################################################################################################

#####################################################################################################
# building rhs with vector potential field found before

# convolution of curl with test function and verse (curl(\tau), a) and (curl(a), \tau)
TAU_CURL_VEC_EXP='[Grad_vec(3,2)-Grad_vec(2,3);Grad_vec(1,3)-Grad_vec(3,1); Grad_vec(2,1)-Grad_vec(1,2)].Test_a'
CURL_TAU_VEC_EXP='[Grad_Test_a(3,2)-Grad_Test_a(2,3);Grad_Test_a(1,3)-Grad_Test_a(3,1); Grad_Test_a(2,1)-Grad_Test_a(1,2)].vec'
tm = time.time()
A_RHS=gf.asm('generic', mim, 1, CURL_TAU_VEC_EXP if WITH_B_CROSS_N else TAU_CURL_VEC_EXP, -1, 'vec', 0, mfb, mdmag.variable('b'), 'a', 1, mfb, 0)
tm = time.time()-tm
print("A . Curl(tau) term elapsed time", tm)
mdpot.add_explicit_rhs('a', A_RHS)
WITH_A_CROSS_N=True
if(WITH_A_CROSS_N):
    A_CROSS_N_RHS_EXP='vec.({S})'.format(S=VPRODUCT_VAR_N_.format(var='Test_a'))
    tm = time.time()
    A_CROSS_N_RHS = gf.asm('generic', mim, 1, A_CROSS_N_RHS_EXP, ALL_SIDES, 'vec', 0, mfb, mdmag.variable('b'), 'a', 1, mfb, 0)
    tm = time.time()-tm
    print("A . (N x tau) term elapsed time", tm)
    A_CROSS_N_RHS_NORM=gf.compute_L2_norm(mfb, A_CROSS_N_RHS, mim)
    print("A . (N x tau) L2_norm", A_CROSS_N_RHS_NORM)
    mdpot.add_explicit_rhs('a', -A_CROSS_N_RHS)
# building rhs with vector potential field found before END
#####################################################################################################

#####################################################################################################
# Assembly of the restoration linear system and solve.		
try:
    mdpot.assembly(option='build_matrix')
    mdpot.assembly(option='build_rhs')
except:
    print("Something went wrong when building matrix")
finally:
    print("curl restoration started at", time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime()), 'total dof number', mdmag.rhs().size)
    tm = time.time()
    #mdmag.solve()
    aa=gf.linsolve_superlu(mdpot.tangent_matrix(),mdpot.rhs())[0][:,0]
    mdpot.set_variable('a', aa[:A_RHS.size])
    tm = time.time()-tm
    print("curl restoration elapsed time %s " % datetime.timedelta(seconds=tm))
# Assembly of the restoration linear system and solve END
#####################################################################################################

#export final solution
mfb.export_to_vtk('restored3d.vtk', 'ascii', mdpot.variable('a'))

print("B L2_norm", gf.compute_L2_norm(mfb, mdpot.variable('a'), mim))
print("B L2_norm ratio", gf.compute_L2_norm(mfb, mdpot.variable('a')-Bexact, mim)/gf.compute_L2_norm(mfb, Bexact, mim))
if(VERIFY_NEIGHBOR_COMPUTATION):
    m.set_region(42,m.inner_faces())
    cross_n_conv=gf.asm('generic', mim, 1, '({P}).({Q})'.format(P=cross_jump_a, Q=cross_jump_test_a), 42, mdpot)
    print('(B+ x N+)-(B- x N-)',np.linalg.norm(cross_n_conv))
    dot_n_conv=gf.asm('generic', mim, 1, '({S}).({T})'.format(S=jump_a,T=jump_test_a), 42, mdpot)
    print('(B+ . N+)-(B- . N-)',np.linalg.norm(dot_n_conv))

####################################################################
# done with restoring target field from the vector potential 
####################################################################
input("Press Enter to continue...")
