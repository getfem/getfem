#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Python GetFEM interface
#
# Copyright (C) 2018-2020 Yves Renard.
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
"""  Test of the assembly on the direct product of two domains.

  This program is used to check that Python-GetFEM interface, and more
  generally GetFEM are working. It focuses on testing some operations
  of the high generic assembly language.

  $Id$
"""
import numpy as np

import getfem as gf

NX = 4

m1 =  gf.Mesh('cartesian', np.arange(0,1+1./NX,1./NX)) # Structured 1D mesh
mfu1 = gf.MeshFem(m1, 1); mfu1.set_fem(gf.Fem('FEM_PK(1,1)'))
mim1 = gf.MeshIm(m1, gf.Integ('IM_GAUSS1D(4)'))
U1 = mfu1.eval('x')


m2 = gf.Mesh('triangles grid', np.arange(0,1+1./NX,1./NX),
            np.arange(0,1+1./NX,1./NX))     # Structured 2D mesh
mfu2 = gf.MeshFem(m2, 1); mfu2.set_fem(gf.Fem('FEM_PK(2,1)'))
mim2 = gf.MeshIm(m2, gf.Integ('IM_TRIANGLE(4)'))
U2 = mfu2.eval('x+y')



md = gf.Model('real')

md.add_fem_variable('u1', mfu1)
md.set_variable('u1', U1)
md.add_fem_variable('u2', mfu2)
md.set_variable('u2', U2)

md.add_standard_secondary_domain('secdom', mim2);


result = gf.asm('generic', mim1, 0, "u1*Secondary_Domain(u2)", -1, md,
                "Secondary_domain", "secdom")
if (abs(result-0.5) > 1e-8) : print("Bad value"); exit(1)

result = gf.asm('generic', mim1, 0, "u1*(Secondary_Domain(Grad_u2)(1))", -1, md,
                "Secondary_domain", "secdom")
if (abs(result-0.5) > 1e-8) : print("Bad value"); exit(1)

result = gf.asm('generic', mim1, 0,
                "u1*(Secondary_Domain(X)(1)+Secondary_Domain(X)(2))", -1, md,
                "Secondary_domain", "secdom")
if (abs(result-0.5) > 1e-8) : print("Bad value"); exit(1)

result = gf.asm('generic', mim1, 0, "u1*sqr(Secondary_Domain(u2))", -1, md,
                "Secondary_domain", "secdom")
if (abs(result-7./12.) > 1e-8) : print("Bad value"); exit(1)





M = gf.asm('generic', mim1, 2, "Test_u1*Secondary_Domain(Test2_u2)", -1, md,
           "Secondary_domain", "secdom")
V1 = gf.asm('generic', mim1, 1, "Test_u1", -1, md)
V2 = gf.asm('generic', mim2, 1, "Test_u2", -1, md)

for i in range(0, V1.size):
  for j in range(0, V2.size):
    if (abs(M[i,j] - V1[i]*V2[j]) > 1E-8):
      print("Bad value for matrix assembly"); exit(1)
