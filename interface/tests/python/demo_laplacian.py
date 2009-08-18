# Python GetFEM++ interface
#
# Copyright (C) 2004-2009 Yves Renard, Julien Pommier.
#                                                       
# This file is a part of GETFEM++                                         
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

## 2D Poisson problem test. 
from numpy import *
from getfem import *

## Parameters
NX=10                               # Mesh parameter.
Dirichlet_with_multipliers = True;  # Dirichlet condition with multipliers
                                    # or penalization
dirichlet_coefficient = 1e10;       # Penalization coefficient

m=Mesh('regular simplices', arange(0,1.01,1./NX), arange(0,1.01,1./NX))
mfu = MeshFem(m, 1)
mfu.set_fem(Fem('FEM_PK(2,2)'))
mfrhs = MeshFem(m, 1)
mfrhs.set_fem(Fem('FEM_PK(2,2)'))
mim = MeshIm(m, Integ('IM_TRIANGLE(4)'))

## Boundary selection
flst = m.outer_faces();
fnor = m.normal_of_faces(flst);
tleft = abs(fnor[1,:]+1) < 1e-14
ttop  = abs(fnor[0,:]-1) < 1e-14
fleft = compress(tleft, flst, axis=1);
ftop  = compress(ttop, flst, axis=1);
fneum = compress(True - ttop - tleft, flst, axis=1);
DIRICHLET_BOUNDARY_NUM1 = 1;
m.set_region(DIRICHLET_BOUNDARY_NUM1, fleft);
DIRICHLET_BOUNDARY_NUM2 = 2;
m.set_region(DIRICHLET_BOUNDARY_NUM2, ftop);
NEUMANN_BOUNDARY_NUM = 3
m.set_region(NEUMANN_BOUNDARY_NUM, fneum);


## Model
md = Model('real')

## Main unknown
md.add_fem_variable("u", mfu)

## Interpolate the exact solution (Assuming mfu is a Lagrange fem)
Uexact = mfu.eval('x[1]*(x[1]-1)*x[0]*(x[0]-1)+x[0]*x[0]*x[0]*x[0]*x[0]');

## Laplacian term on u.
md.add_Laplacian_brick(mim, "u")

## Volumic source term.
md.add_initialized_fem_data("VolumicData", mfrhs,
   mfrhs.eval('-(2*(x[0]*x[0]+x[1]*x[1])-2*x[0]-2*x[1]+20*x[0]*x[0]*x[0])'))
md.add_source_term_brick(mim, "u", "VolumicData")

## Neumann condition.
md.add_initialized_fem_data("NeumannData", mfrhs,
   mfrhs.eval('[x[1]*(x[1]-1)*(2*x[0]-1) + 5*x[0]*x[0]*x[0]*x[0], x[0]*(x[0]-1)*(2*x[1]-1)]'))
md.add_normal_source_term_brick(mim, "u", "NeumannData", NEUMANN_BOUNDARY_NUM)

## Dirichlet condition on the left.
md.add_initialized_fem_data("DirichletData", mfrhs,
   mfrhs.eval('x[1]*(x[1]-1)*x[0]*(x[0]-1)+x[0]*x[0]*x[0]*x[0]*x[0]'))
if (Dirichlet_with_multipliers):
    md.add_Dirichlet_condition_with_multipliers(mim, "u", mfu,
                                                DIRICHLET_BOUNDARY_NUM1,
                                                "DirichletData")
else:
    md.add_Dirichlet_condition_with_penalization(mim, "u",
                                                 dirichlet_coefficient,
                                                 DIRICHLET_BOUNDARY_NUM1,
                                                 "DirichletData")
## Dirichlet condition on the top.
## Two Dirichlet brick in order to test the multiplier selection in the
##     intersection.
if (Dirichlet_with_multipliers):
    md.add_Dirichlet_condition_with_multipliers(mim, "u", mfu,
                                                DIRICHLET_BOUNDARY_NUM2,
                                                "DirichletData")
else:
    md.add_Dirichlet_condition_with_penalization(mim, "u",
                                                 dirichlet_coefficient,
                                                 DIRICHLET_BOUNDARY_NUM2,
                                                 "DirichletData")

# md.listvar();
# md.listbricks();

## Assembly of the linear system and solve.
md.solve()

## Main unknown
U = md.variable("u")
L2error = compute(mfu, U-Uexact, 'L2 norm', mim)
H1error = compute(mfu, U-Uexact, 'H1 norm', mim)

if (H1error > 0.01):
    print 'Error in L2 norm : ', L2error
    print 'Error in H1 norm : ', H1error
    print 'Error too large !'
    exit(1);
