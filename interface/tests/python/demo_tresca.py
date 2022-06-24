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
"""  2D Poisson problem test with a Tresca friction condition.

  This program is used to check that python-getfem is working. This is
  also a good example of use of GetFEM.

  $Id$
"""
import numpy as np

# Import basic modules
import getfem as gf

## Parameters
NX = 100                           # Mesh parameter.
thr = 1                            # Tresca threshold
gamma0 = 5                         # Nitsche's parameter
theta = 1                          # Nitsche's version

# Create a simple cartesian mesh
m = gf.Mesh('regular_simplices', np.arange(0,2+1./NX,1./NX),
            np.arange(0,1+1./NX,1./NX))

# Create a MeshFem for u and rhs fields of dimension 1 (i.e. a scalar field)
mfu   = gf.MeshFem(m, 1)
mfrhs = gf.MeshFem(m, 1)
# assign the P2 fem to all convexes of the both MeshFem
mfu.set_fem(gf.Fem('FEM_PK(2,2)'))
mfrhs.set_fem(gf.Fem('FEM_PK(2,2)'))

#  Integration method used
mim = gf.MeshIm(m, gf.Integ('IM_TRIANGLE(4)'))

# Boundary selection
flst  = m.outer_faces()
fnor  = m.normal_of_faces(flst)
ftop  = np.compress(abs(fnor[1,:]-1) < 1e-14, flst, axis=1)
fbottom  = np.compress(abs(fnor[1,:]+1) < 1e-14, flst, axis=1)

# Mark it as boundary
DIRICHLET_BOUNDARY = 1
CONTACT_BOUNDARY = 3
m.set_region(DIRICHLET_BOUNDARY, ftop)
m.set_region(CONTACT_BOUNDARY, fbottom)


# Interpolate the source term
F = mfrhs.eval('x+3*y/2')

# Model
md = gf.Model('real')
md.add_initialized_data("theta", [theta])
md.add_initialized_data("thr", [thr])
md.add_initialized_data("gamma0", [gamma0])
md.add_fem_variable('u', mfu)
md.add_macro("gh", "gamma0/element_size");
md.add_macro("dnu", "Grad_u.Normal");
md.add_macro("dnv", "Grad_Test_u.Normal");
# md.add_macro("P(x,y)", "Ball_projection(x,y)"); # Do not work in old versions
md.add_macro("P(x,y)", "max(min(x,y),-y)"); # Projection definition

# Laplacian term on u
md.add_linear_term(mim, 'Grad_u.Grad_Test_u')

# Nitsche terms for the Tresca condition
md.add_nonlinear_term(mim, '-(theta/gh)*dnu*dnv',
                      CONTACT_BOUNDARY)
md.add_nonlinear_term(mim, '(1/gh)*P(dnu-gh*u,thr)*(theta*dnv-gh*Test_u)',
                      CONTACT_BOUNDARY)



# Volumic source term
md.add_initialized_fem_data('VolumicData', mfrhs, F)
md.add_source_term_brick(mim, 'u', 'VolumicData')

# Dirichlet condition on the top.
md.add_Dirichlet_condition_with_multipliers(mim, 'u', mfu, DIRICHLET_BOUNDARY)


# Assembly of the linear system and solve.
md.solve('max_res', 1e-8, 'max_iter', 100, 'noisy')

# Main unknown
U = md.variable('u')

# Export data
mfu.export_to_pos('tresca.pos', U,'Computed solution')
mfu.export_to_vtk('tresca.vtk', U,'Computed solution')
print('You can view the solution with (for example):')
print('gmsh tresca.pos')
print("mayavi2 -d tresca.vtk -f WarpScalar -m Surface")
