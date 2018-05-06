#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Python GetFEM++ interface
#
# Copyright (C) 2018-2018 Yves Renard.
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
"""  test high generic assembly language.

  This program is used to check that Python-GetFEM interface, and more
  generally GetFEM are working. It focuses on testing some operations
  of the high generic assembly language.

  $Id$
"""
import numpy as np
import getfem as gf
import os


NX = 4
m = gf.Mesh('triangles grid', np.arange(0,1+1./NX,1./NX),
            np.arange(0,1+1./NX,1./NX))     # Structured mesh
fem = gf.Fem('FEM_PK(2,1)')
mfu = gf.MeshFem(m, 1); mfu.set_fem(fem)    # Lagrange P1 scalar fem
mim = gf.MeshIm(m, gf.Integ('IM_TRIANGLE(4)'))
mfv = gf.MeshFem(m, 3); mfv.set_fem(fem)    # Lagrange P1 vector fem


U = mfu.eval('x+y')
V = mfv.eval('[x*y, x*y, x*y]')


md = gf.Model('real')

md.add_fem_variable('u', mfu)
md.set_variable('u', U)
md.add_fem_variable('v', mfv)
md.set_variable('v', V)

# Simple test on the integral of u
result = gf.asm('generic', mim, 0, "u", -1, md)
if (abs(result-1) > 1e-8) : print "Bad value"; exit(1)

# Single contraction and comparison with Trace
result1 = gf.asm('generic', mim, 0,
                 "Def P(a):=a*(a'); Contract(P(Grad_v), 1, 2)", -1, md)
result2 = gf.asm('generic', mim, 0,
                 "Def P(a):=a*(a'); Trace(P(Grad_v))", -1, md)
if (abs(result1-result2) > 1e-8) : print "Bad value"; exit(1)

# Constant order 3 tensor contraction test
result1 = gf.asm('generic', mim, 0,
                 "Contract([[[1,1],[2,2]],[[1,1],[2,2]]], 1, 2)", -1, md)
result2 = np.array([3., 3.]);
if (np.linalg.norm(result1-result2) > 1e-8) : print "Bad value"; exit(1)

# Single contraction, comparison with "*"
result1 = gf.asm('generic', mim, 0, "Contract(Grad_v, 2, Grad_u, 1)", -1, md)
result2 = gf.asm('generic', mim, 0, "Grad_v * Grad_u", -1, md)
if (np.linalg.norm(result1-result2) > 1e-8) : print "Bad value"; exit(1)

# Double contraction order one expression, comparison with ":"
result1 = gf.asm('generic', mim, 1,
                 "Contract(Grad_v, 1, 2, Grad_Test_v, 1, 2)", -1, md)
result2 = gf.asm('generic', mim, 1, "Grad_v : Grad_Test_v", -1, md)
if (np.linalg.norm(result1-result2) > 1e-8) : print "Bad value"; exit(1)

# Double contraction order two expression, comparison with ":"
result1 = gf.asm('generic', mim, 2,
                 "Contract(Grad_Test2_v, 1, 2, Grad_Test_v, 1, 2)", -1, md)
result2 = gf.asm('generic', mim, 2, "Grad_Test2_v : Grad_Test_v", -1, md)
if (np.linalg.norm(result1.full()-result2.full()) > 1e-8) :
  print "Bad value"; exit(1)
result1 = gf.asm('generic', mim, 2,
                 "Contract(Grad_Test_v, 2, 1, Grad_Test2_v, 2, 1)", -1, md)
if (np.linalg.norm(result1.full()-result2.full()) > 1e-8) :
  print "Bad value"; exit(1)




