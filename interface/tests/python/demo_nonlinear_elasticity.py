#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Python GetFEM++ interface
#
# Copyright (C) 2012-2016 Yves Renard.
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
from numpy import linalg as npla

gf.util_trace_level(1)

dirichlet_version = 2          # 1 = simplification, 2 = penalisation
test_tangent_matrix = False    # Test or not tangent system validity
incompressible = True;         # Incompressibility option

# lawname = 'Ciarlet Geymonat'
# params = [1.,1.,0.25]
lawname = 'SaintVenant Kirchhoff'
params = [1.,1.]
if (incompressible):
    lawname = 'Incompressible Mooney Rivlin'
    params = [1.,1.]

N1 = 2; N2 = 4; h = 20.; DX = 1./N1; DY = (1.*h)/N2;
m = gf.Mesh('cartesian', np.arange(-0.5, 0.5+DX,DX), np.arange(0., h+DY,DY),
            np.arange(-1.5, 1.5+3*DX,3*DX))
mfu  = gf.MeshFem(m, 3)            # mesh-fem supporting a 3D-vector field
mfdu = gf.MeshFem(m,1)
# The mesh_im stores the integration methods for each tetrahedron
mim = gf.MeshIm(m, gf.Integ('IM_GAUSS_PARALLELEPIPED(3,4)'))
# We choose a P2 fem for the main unknown
mfu.set_fem(gf.Fem('FEM_QK(3,2)'))

if (dirichlet_version == 1):
  mfd = mfu;
else:
  mfd = gf.MeshFem(m,1)
  mfd.set_fem(gf.Fem('FEM_QK(3,1)'))

# The P2 fem is not derivable across elements, hence we use a discontinuous
# fem for the derivative of U.
mfdu.set_fem(gf.Fem('FEM_QK_DISCONTINUOUS(3,2)'));

# Display some information about the mesh
print('nbcvs=%d, nbpts=%d, nbdof=%d' % (m.nbcvs(), m.nbpts(), mfu.nbdof()))

# Assign boundary numbers

ftop = m.outer_faces_with_direction([0., 0., 1.], 0.5)
fbot = m.outer_faces_with_direction([0., 0., -1.], 0.5)

m.set_region(1, ftop);
m.set_region(2, ftop);
m.set_region(3, np.append(ftop,fbot,axis=1));

# Model definition

md=gf.Model('real')
md.add_fem_variable('u', mfu)
md.add_initialized_data('params', params)
md.add_finite_strain_elasticity_brick(mim, 'u', lawname, 'params')

# md.add_nonlinear_generic_assembly_brick(mim, 'sqr(Trace(Green_Lagrangian(Id(meshdim)+Grad_u)))/8 + Norm_sqr(Green_Lagrangian(Id(meshdim)+Grad_u))/4')
# md.add_nonlinear_generic_assembly_brick(mim, '((Id(meshdim)+Grad_u)*(params(1)*Trace(Green_Lagrangian(Id(meshdim)+Grad_u))*Id(meshdim)+2*params(2)*Green_Lagrangian(Id(meshdim)+Grad_u))):Grad_Test_u')
# md.add_nonlinear_generic_assembly_brick(mim, 'Saint_Venant_Kirchhoff_potential(Grad_u,params)')
    
if (incompressible):
    mfp = gf.MeshFem(m,1)
    mfp.set_classical_discontinuous_fem(1)
    md.add_fem_variable('p', mfp)
    md.add_finite_strain_incompressibility_brick(mim, 'u', 'p')
    # md.add_nonlinear_generic_assembly brick(mim, 'p*(1-Det(Id(meshdim)+Grad_u))')
    # md.add_nonlinear_generic_assembly_brick(mim, '-p*Det(Id(meshdim)+Grad_u)*(Inv(Id(meshdim)+Grad_u))'':Grad_Test_u + Test_p*(1-Det(Id(meshdim)+Grad_u))')


if (dirichlet_version == 1):
    md.add_fem_data('DirichletData', mfu)
    md.add_Dirichlet_condition_with_simplification('u', 3, 'DirichletData')
else:
    md.add_fem_data('DirichletData', mfd, 3)
    md.add_Dirichlet_condition_with_penalization(mim, 'u', 1e4, 3, 'DirichletData')


VM=np.zeros(mfdu.nbdof());
UU=np.zeros(0);
VVM=np.zeros(0);
nbstep=40;


P = mfd.basic_dof_nodes()
r = np.sqrt(np.square(P[0 ,:]) + np.square(P[2, :]))
theta = np.arctan2(P[2,:], P[0,:]);


def axrot_matrix(A, B, theta):
    n=(np.array(B)-np.array(A)); n = n/npla.norm(n)
    a=n[0]; b=n[1]; c=n[2]
    d=np.sqrt(b*b+c*c)
    T=np.eye(4); T[0:3,3]=-np.array(A)
    Rx=np.eye(4)
    if (npla.norm(n[1:3])>1e-6):
        Rx[1:3,1:3]=np.array([[c,-b],[b,c]]/d)
    Ry=np.eye(4); Ry[[[0],[2]],[0,2]]=[[d,-a],[a,d]]
    Rz=np.eye(4); Rz[0:2,0:2]=np.array([[np.cos(theta),np.sin(theta)],[-np.sin(theta),np.cos(theta)]])
    R = npla.inv(T)*npla.inv(Rx)*npla.inv(Ry)*Rz*Ry*Rx*T
    return R;
  



for step in range(1,nbstep+1):
    w = 3.*step/nbstep

    # Computation of the rotation for Dirichlet's condition
    dtheta =  np.pi
    dtheta2 = np.pi/2
      
    if (dirichlet_version == 1):
        R=np.zeros(mfd.nbdof())
    else:
        R=np.zeros((3, mfd.nbdof()))
    
    i_top = mfd.basic_dof_on_region(1)
    i_bot = mfd.basic_dof_on_region(2)
    
    dd = np.amax(P[1,i_top]*np.sin(w*dtheta))
    
    if (w < 1): 
        RT1 = axrot_matrix([0, h*.75, 0], [0, h*.75, 1], w*dtheta)
        RT2 = axrot_matrix([0, 0, 0], [0, 1, 0], np.sqrt(w)*dtheta2)
        RB1 = axrot_matrix([0, h*.25, 0], [0, h*.25, 1], -w*dtheta)
        RB2 = RT2.transpose()
    elif (w < 2):
        RT1 = axrot_matrix([0, h*.75, 0], [0, h*.75, 1], (2-w)*dtheta);
        RT2 = axrot_matrix([0, 0, 0], [0, 1, 0], w*dtheta2);
        RB1 = axrot_matrix([0, h*.25, 0], [0, h*.25, 1], -(2-w)*dtheta);
        RB2 = RT2.transpose()
    else:
        RT1 = axrot_matrix([0, h*.75, 0], [0, h*.75, 1], 0);
        RT2 = axrot_matrix([0, 0, 0], [0, 1, 0], (3-w)*2*dtheta2);
        RB1 = axrot_matrix([0, h*.25, 0], [0, h*.25, 1], 0);
        RB2 = RT2.transpose()


    if (dirichlet_version == 1):
##       for i=i_top,
##         ro = RT1*RT2*[P(:,i);1];
##         R(i) = ro(1+mod(i-1,3)) - P(1+mod(i-1,3),i);
##       end
##       for i=i_bot,
##         ro = RB1*RB2*[P(:,i);1];
##         R(i) = ro(1+mod(i-1,3)) - P(1+mod(i-1,3),i);
##       end 
    else:
##       for i=i_top,
##         ro = RT1*RT2*[P(:,i);1];
##         R(:, i) = ro(1:3) - P(:,i);
##       end
##       for i=i_bot,
##         ro = RB1*RB2*[P(:,i);1];
##         R(:, i) = ro(1:3) - P(:,i);
##       end
    
    md.set_variable('DirichletData', R)
    md.solve('very noisy', 'max_iter', 100, 'max_res', 1e-5, 'lsearch', 'simplest')

    print("Iteration %d done" % step)

    if (test_tangent_matrix):
       gf_model_get(md.test_tangent_matrix(1E-8, 10, 0.0001)

       
##     U = gf_model_get(md, 'variable', 'u');
##     VM0 = gf_model_get(md, 'compute Von Mises or Tresca', 'u', lawname, 'params', mfdu);
##     # sigma = gf_model_get(md, 'compute second Piola Kirchhoff tensor', 'u', lawname, 'params', mfdu);
    
##     # Direct interpolation of the Von Mises stress
##     # VM = gf_model_get(md, 'interpolation', '(sqrt(3/2)/Det(Id(meshdim)+Grad_u))*Norm((Id(meshdim)+Grad_u)*Saint_Venant_Kirchhoff_sigma(Grad_u,params)*(Id(meshdim)+Grad_u'') - Id(meshdim)*Trace((Id(meshdim)+Grad_u)*Saint_Venant_Kirchhoff_sigma(Grad_u,params)*(Id(meshdim)+Grad_u''))/meshdim)', mfdu);
##     # VM = gf_model_get(md, 'interpolation', '(sqrt(3/2)/Det(Id(meshdim)+Grad_u))*Norm(Deviator((Id(meshdim)+Grad_u)*Saint_Venant_Kirchhoff_sigma(Grad_u,params)*(Id(meshdim)+Grad_u'')))', mfdu);
##     # VM = gf_model_get(md, 'interpolation', 'sqrt(3/2)*Norm(Deviator(Cauchy_stress_from_PK2(Saint_Venant_Kirchhoff_sigma(Grad_u,params),Grad_u)))', mfdu);
##     VM = gf_model_get(md, 'finite strain elasticity Von Mises', 'u', lawname, 'params', mfdu);
##     norm(VM-VM0)
    
##     UU = [UU;U]; 
##     VVM = [VVM;VM];
##     save demo_nonlinear_elasticity_U.mat UU VVM m_char mfu_char mfdu_char;
##   else
##     U=UU(step,:);
##     VM=VVM(step,:);
##   end;
##   disp(sprintf('step %d/%d : |U| = %g',step,nbstep,norm(U)));

##   if (drawing)
##     gf_plot(mfdu,VM,'mesh','off', 'cvlst',gf_mesh_get(mfdu,'outer faces'), 'deformation',U,'deformation_mf',mfu,'deformation_scale', 1, 'refine', 8); colorbar;
##     axis([-3     6     0    20    -2     2]); caxis([0 .3]);
##     view(30+20*w, 23+30*w);  
##     campos([50 -30 80]);
##     camva(8);
##     camup
##     camlight; 
##     axis off;
##     pause(1);
    
##   end
##   end;
















## K = 1.2; mu = 3.0;
## _F_ = "(Id(3)+Grad_u)"
## _J_= "Det{F}".format(F=_F_)
## _be_ = "(Left_Cauchy_Green{F})".format(F=_F_)

## _expr_1 = "{K_over_2}*sqr(log({J}))+{mu_over_2}*(Matrix_j1{be}-3)"\
##          .format(K_over_2=K/2., J=_J_, mu_over_2=mu/2., be=_be_)


## _expr_2 = "{K_over_2}*sqr(log({J}))+{mu_over_2}*(pow(Det{be},-1./3.)*Trace{be}-3)"\
##          .format(K_over_2=K/2., J=_J_, mu_over_2=mu/2., be=_be_)

## # Mettre des prints sur la différence entre les deux ... + version avec potentiel de la loi d'élasticité non linéaire ... plutot recopier exemple elast. nonlinéaire


## print("_expr1_ = %s" % _expr_1);

## print("_expr2_ = %s" % _expr_2);


## exit(1);
