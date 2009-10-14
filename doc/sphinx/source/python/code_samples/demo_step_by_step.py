#!/usr/bin/env python
# -*- coding: UTF8 -*-

# import basic modules
import getfem as gf
import numpy as np

# creation of a simple cartesian mesh
m = gf.Mesh('cartesian', np.arange(0,1.1,0.1), np.arange(0,1.1,0.1))

# create a MeshFem of for a field of dimension 1 (i.e. a scalar field)
mf = gf.MeshFem(m, 1)
# assign the Q2 fem to all convexes of the MeshFem
mf.set_fem(gf.Fem('FEM_QK(2,2)'))

# view the expression of its basis functions on the reference convex
print gf.Fem('FEM_QK(2,2)').poly_str()

# an exact integration will be used
mim = gf.MeshIm(m, gf.Integ('IM_EXACT_PARALLELEPIPED(2)'))

# detect the border of the mesh
border = m.outer_faces()
# mark it as boundary #42
m.set_region(42, border)

# generic elliptic brick
b0 = gf.MdBrick('generic elliptic', mim, mf)

# list of parameters of the brick
print b0.param_list()

# add a Dirichlet condition on the domain boundary
b1 = gf.MdBrick('dirichlet', b0, 42, mf, 'penalized')

# change Dirichlet condition
R = mf.eval('x[0]*(x[0]-1) - x[1]*(x[1]-1)')
b1.set_param('R', mf, R)

# created model state
mds = gf.MdState('real')
b1.solve(mds)

# extracted solution
sol = mds.get('state')

# export computed solution
mf.export_to_pos('sol.pos',sol,'Computed solution')
