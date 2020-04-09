#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Python GetFEM interface
#
# Copyright (C) 2004-2020 Yves Renard, Julien Pommier.
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
from numpy import *

from getfem import *

print('3D stokes demonstration on a quadratic mesh')

viscosity = 10


m=Mesh('import','GiD','../meshes/tank_quadratic_2500.GiD.msh')
print('mesh loaded!')
mfu=MeshFem(m,3) # velocity
mfulag=MeshFem(m,3)
mfp=MeshFem(m,1) # pressure
mfd=MeshFem(m,1) # data
mfe=MeshFem(m,1)
mim=MeshIm(m, Integ('IM_TETRAHEDRON(5)'))

mfu.set_fem(Fem('FEM_PK(3,2)'))
mfd.set_fem(Fem('FEM_PK(3,2)'))
mfp.set_fem(Fem('FEM_PK(3,1)'))
mfe.set_fem(Fem('FEM_PK_DISCONTINUOUS(3,1,0.01)'))

print('nbcvs=%d, nbpts=%d, qdim=%d, fem = %s, nbdof=%d' % \
      (m.nbcvs(), m.nbpts(), mfu.qdim(), mfu.fem()[0].char(), mfu.nbdof()))


P=m.pts()
r = list(range(0, m.nbpts()));
INpid=compress(abs(P[0,:]+25) < 1e-4, r)
OUTpid=compress(abs(P[0,:]-25) < 1e-4, r)
TOPpid=compress(abs(P[2,:]-20) < 1e-4, r)
INfaces=m.faces_from_pid(INpid)
OUTfaces=m.faces_from_pid(OUTpid)
TOPfaces=m.faces_from_pid(TOPpid)

m.set_region(1, INfaces)
m.set_region(2, OUTfaces)
m.set_region(3, TOPfaces)
m.set_region(4, m.outer_faces())
m.region_subtract(4, 1)
m.region_subtract(4, 2)
m.region_subtract(4, 3)



md=Model('real');
md.add_fem_variable('u', mfu);
md.add_initialized_data('lambda', [0]);
md.add_initialized_data('mu', [viscosity]);
md.add_isotropic_linearized_elasticity_brick(mim, 'u', 'lambda', 'mu');
md.add_fem_variable('p', mfp);
md.add_linear_incompressibility_brick(mim, 'u', 'p');
md.add_variable('mult_spec', 1);
M = Spmat('empty', 1, mfp.nbdof());
M.add(list(range(1)), list(range(mfp.nbdof())), ones((1, mfp.nbdof())));
md.add_constraint_with_multipliers('p', 'mult_spec', M, [0]);
md.add_initialized_data('NeumannData', [0, -10, 0]);
md.add_source_term_brick(mim, 'u', 'NeumannData', 1);

D = mfd.basic_dof_nodes();
x = D[0,:]
y = D[1,:]
z = D[2,:]
md.add_initialized_fem_data('Dir1data', mfd, [9-(y*y+(z-6)*(z-6)),0*x,0*x]);
md.add_Dirichlet_condition_with_multipliers(mim, 'u', mfu, 1, 'Dir1data');
md.add_initialized_fem_data('Dir2data',  mfd, [9-(y*y+(z-6)*(z-6)),0*x,0*x]);
md.add_Dirichlet_condition_with_multipliers(mim, 'u', mfu, 2, 'Dir2data');
md.add_Dirichlet_condition_with_multipliers(mim, 'u', mfu, 3);
md.add_Dirichlet_condition_with_multipliers(mim, 'u', mfu, 4);

print('running solve...')
md.solve('noisy', 'lsolver','superlu')
print('solve done!')


U = md.variable('u');
P = md.variable('p');
VM = md.compute_isotropic_linearized_Von_Mises_or_Tresca('u', 'lambda',
                                                         'mu', mfe);

mfu.save('tank_3D.mfu', 'with_mesh')
mfp.save('tank_3D.mfp', 'with_mesh')
U.tofile('tank_3D.U')
P.tofile('tank_3D.P')

mfe.save('tank_3D.mfe')
VM.tofile('tank_3D.VM')
#memstats()
