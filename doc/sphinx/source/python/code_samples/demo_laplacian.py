#!/usr/bin/env python
# -*- coding: UTF8 -*-

# import basic modules
import getfem
from numpy import *

# creation of a simple cartesian mesh
m = getfem.Mesh('cartesian', arange(0,1.1,0.1), arange(0,1.1,0.1))

#  create a MeshFem for a field of dimension 1 (i.e. a scalar field)
mf = getfem.MeshFem(m, 1)
# assign the Q2 fem to all convexes of the MeshFem
mf.set_fem(getfem.Fem('FEM_QK(2,2)'))

# an exact integration will be used
mim = getfem.MeshIm(m, getfem.Integ('IM_EXACT_PARALLELEPIPED(2)'))

# detect the border of the mesh
border = m.outer_faces()
# mark it as boundary #42
m.set_region(42, border)

# interpolate the exact solution
R = mf.eval('x[0]*(x[0]-1)*x[1]*(x[1]-1)/2.0')
# interpolate the source term
F = mf.eval('-( x[0]*(x[0]-1) + x[1]*(x[1]-1) )')

# generic elliptic brick
b0 = getfem.MdBrick('generic elliptic', mim, mf)
# add a Dirichlet condition on the domain boundary
b1 = getfem.MdBrick('dirichlet', b0, 42, mf, 'penalized')
b1.set_param('R', mf, R)
# add a source term
b2 = getfem.MdBrick('source term', b1)
b2.set_param('source_term', mf, F)

# model state
mds = getfem.MdState(b2)
b2.solve(mds)

# extracted solution
sol = mds.get('state')

# export data
mf.export_to_pos('sol.pos', R,'Exact solution',
                            sol,'Computed solution',
                            abs(sol-R),'abs differences')
