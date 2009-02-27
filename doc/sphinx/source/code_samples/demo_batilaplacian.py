#!/usr/bin/env python
# -*- coding: UTF8 -*-

# import basic modules
import getfem
from numpy import *

# boundary names
up    = 101 # m.region(101)
down  = 102 # m.region(102)
left  = 103 # m.region(103)
right = 104 # m.region(104)

# importing the mesh and boundary
m = getfem.Mesh('import','gmsh','quad.msh')

# create a MeshFem of for a field of dimension 1 (i.e. a scalar field)
mf = getfem.MeshFem(m,1)
# assign the P1 fem to all convexes of the MeshFem
mf.set_fem(getfem.Fem('FEM_PK(2,1)'))

# an exact integration will be used
mim = getfem.MeshIm(m,getfem.Integ('IM_EXACT_SIMPLEX(2)'))

# interpolate the exact solution
R = mf.eval('x[0]*(x[0]-1)*x[1]*(x[1]-1)')
# its Normal derivative on the domain boundary right (left is -DU)
DU = mf.eval('(2*x[0]-1)*x[1]*(x[1]-1)')
# its laplacian
F = mf.eval('-(2*(x[0]**2+x[1]**2)-2*(x[0]+x[1]))')

# generic elliptic brick
b0 = getfem.MdBrick('generic elliptic',mim,mf)
# add a Dirichlet condition on the domain boundary up
b1 = getfem.MdBrick('dirichlet',b0,up,mf,'penalized')
b1.set_param('R',mf,Uexact)
# add a Dirichlet condition on the domain boundary down
b2 = getfem.MdBrick('dirichlet',b1,down,mf,'penalized')
b2.set_param('R',mf,Uexact)
# add a source term
b3 = getfem.MdBrick('source term',b2)
b3.set_param('source_term',mf,F)
# add a Neumann condition on the domain boundary left
b4 = getfem.MdBrick('source_term',b3,left)
b4.set_param('source term',mf,-uN)
# add a Neumann condition on the domain boundary right
b5 = getfem.MdBrick('source_term',b4,right)
b5.set_param('source term',mf,uN)

# model state
mds = getfem.MdState('real')
b5.solve(mds)

# extracted solution
csol = array([mds.get('state')])
vsol = concatenate((0*csol,0*csol,csol),axis=0)

# export
mf.export_to_pos('sol.pos',Uexact,'Exact solution',
                           csol,'Computed solution',
                           vsol,'Displacement',
                 abs(csol-Uexact),'abs difference with exact solution')
