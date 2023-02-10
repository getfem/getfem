#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Python GetFEM interface
#
# Copyright (C) 2023-2023 Konstantinos Poulios
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
"""  Test of the gf.MeshFem("bspline_uniform", ...) object.

  This program is used to check that Python-GetFEM interface, and more
  generally GetFEM are working. It focuses on testing the creation of mesh_fem
  with bspline basis functions in 1D and 2D.

  $Id$
"""
import numpy as np
import getfem as gf

NX = 8                            # Number of bspline intervals in x direction
NY = 5                            # Number of bspline intervals in y direction

use_quad = True                  # Quadrilaterals or triangles
export_vtu = False

if (use_quad):
  m1 = gf.Mesh('import','structured','GT="GT_QK(1,1)";SIZES=[12];NOISED=0;NSUBDIV=[%d];' % (2*NX))
  m2 = gf.Mesh('import','structured','GT="GT_QK(2,1)";SIZES=[12,7];NOISED=0;NSUBDIV=[%d,%d];' % (2*NX,2*NY))
else:
  m1 = gf.Mesh('import','structured','GT="GT_PK(1,1)";SIZES=[12];NOISED=0;NSUBDIV=[%d];' % (2*NX))
  m2 = gf.Mesh('import','structured','GT="GT_PK(2,1)";SIZES=[12,7];NOISED=0;NSUBDIV=[%d,%d];' % (2*NX,2*NY))

mim1 = gf.MeshIm(m1, 5)
mim2 = gf.MeshIm(m2, 5)

for order in [3, 4, 5]:
  mf1_free_free = gf.MeshFem("bspline_uniform", m1, NX, order)
  mf1_free_sym = gf.MeshFem("bspline_uniform", m1, NX, order, "free", "symmetry")
  mf1_sym_free = gf.MeshFem("bspline_uniform", m1, NX, order, "symmetry", "free")
  mf1_sym_sym = gf.MeshFem("bspline_uniform", m1, NX, order, "symmetry")
  mf1_periodic = gf.MeshFem("bspline_uniform", m1, NX, order, "periodic")

  mf2_free_free_free_free = gf.MeshFem("bspline_uniform", m2, NX, NY, order)
  mf2_free_sym_sym_free = gf.MeshFem("bspline_uniform", m2, NX, NY, order, "free", "symmetry", "symmetry", "free")
  mf2_sym_sym_sym_sym = gf.MeshFem("bspline_uniform", m2, NX, NY, order, "symmetry", "symmetry")
  mf2_periodic_free_periodic_free = gf.MeshFem("bspline_uniform", m2, NX, NY, order, "periodic", "free")
  mf2_sym_periodic_free_periodic = gf.MeshFem("bspline_uniform", m2, NX, NY, order, "symmetry", "periodic", "free", "periodic")

  print("order %d" % order)
  print("mf1_free_free.nbdof()=", mf1_free_free.nbdof())
  print("mf1_free_sym.nbdof()=", mf1_free_sym.nbdof())
  print("mf1_sym_free.nbdof()=", mf1_sym_free.nbdof())
  print("mf1_sym_sym.nbdof()=", mf1_sym_sym.nbdof())
  print("mf1_periodic.nbdof()=", mf1_periodic.nbdof())
  if export_vtu:
    for dof in range(mf1_free_free.nbdof()):
      mf1_free_free.export_to_vtu(f"basis_funcs_order{order}_free_free_{dof}.vtu",
                                  (np.arange(mf1_free_free.nbdof())==dof).astype(float), "basis function")
    for dof in range(mf1_free_sym.nbdof()):
      mf1_free_sym.export_to_vtu(f"basis_funcs_order{order}_free_sym_{dof}.vtu",
                                 (np.arange(mf1_free_sym.nbdof())==dof).astype(float), "basis function")
    for dof in range(mf1_sym_free.nbdof()):
      mf1_sym_free.export_to_vtu(f"basis_funcs_order{order}_sym_free_{dof}.vtu",
                                 (np.arange(mf1_sym_free.nbdof())==dof).astype(float), "basis function")
    for dof in range(mf1_sym_sym.nbdof()):
      mf1_sym_sym.export_to_vtu(f"basis_funcs_order{order}_sym_sym_{dof}.vtu",
                                (np.arange(mf1_sym_sym.nbdof())==dof).astype(float), "basis function")
    for dof in range(mf1_periodic.nbdof()):
      mf1_periodic.export_to_vtu(f"basis_funcs_order{order}_periodic_{dof}.vtu",
                                 (np.arange(mf1_periodic.nbdof())==dof).astype(float), "basis function")

  print("mf2_free_free_free_free.nbdof()=", mf2_free_free_free_free.nbdof())
  print("mf2_free_sym_sym_free.nbdof()=", mf2_free_sym_sym_free.nbdof())
  print("mf2_sym_sym_sym_sym.nbdof()=", mf2_sym_sym_sym_sym.nbdof())
  print("mf2_periodic_free_periodic_free.nbdof()=", mf2_periodic_free_periodic_free.nbdof())
  print("mf2_sym_periodic_free_periodic.nbdof()=", mf2_sym_periodic_free_periodic.nbdof())
  if export_vtu:
    for dof in range(mf2_free_free_free_free.nbdof()):
      mf2_free_free_free_free.export_to_vtu(f"basis_funcs_order{order}_free_free_free_free_{dof}.vtu",
                                  (np.arange(mf2_free_free_free_free.nbdof())==dof).astype(float), "basis function")
    for dof in range(mf2_free_sym_sym_free.nbdof()):
      mf2_free_sym_sym_free.export_to_vtu(f"basis_funcs_order{order}_free_sym_sym_free_{dof}.vtu",
                   (np.arange(mf2_free_sym_sym_free.nbdof())==dof).astype(float), "basis function")
    for dof in range(mf2_sym_sym_sym_sym.nbdof()):
      mf2_sym_sym_sym_sym.export_to_vtu(f"basis_funcs_order{order}_sym_sym_sym_sym_{dof}.vtu",
                   (np.arange(mf2_sym_sym_sym_sym.nbdof())==dof).astype(float), "basis function")
    for dof in range(mf2_periodic_free_periodic_free.nbdof()):
      mf2_periodic_free_periodic_free.export_to_vtu(f"basis_funcs_order{order}_periodic_free_periodic_free_{dof}.vtu",
                   (np.arange(mf2_periodic_free_periodic_free.nbdof())==dof).astype(float), "basis function")
    for dof in range(mf2_sym_periodic_free_periodic.nbdof()):
      mf2_sym_periodic_free_periodic.export_to_vtu(f"basis_funcs_order{order}_sym_periodic_free_periodic_{dof}.vtu",
                   (np.arange(mf2_sym_periodic_free_periodic.nbdof())==dof).astype(float), "basis function")

  assert(mf1_free_free.nbdof() == (10,11,12)[order-3])
  assert(mf1_free_sym.nbdof() == (9,10,10)[order-3])
  assert(mf1_sym_free.nbdof() == (9,10,10)[order-3])
  assert(mf1_sym_sym.nbdof() == (8,9,8)[order-3])
  assert(mf1_periodic.nbdof() == (8,8,8)[order-3])

  assert(mf2_free_free_free_free.nbdof() == (70,88,108)[order-3]);
  assert(mf2_free_sym_sym_free.nbdof() == (54,70,70)[order-3]);
  assert(mf2_sym_sym_sym_sym.nbdof() == (40,54,40)[order-3]);
  assert(mf2_periodic_free_periodic_free.nbdof() == (56,64,72)[order-3]);
  assert(mf2_sym_periodic_free_periodic.nbdof() == (45,50,50)[order-3]);

  # partition of unity check
  L2error = gf.asm_generic(mim1, 0, "Norm_sqr(1-t)", -1,
                           "t", False, mf1_free_free, np.ones(mf1_free_free.nbdof()))
  assert(L2error < 1e-16);
  print("mf1_free_free partition of unity error:", L2error)
  L2error = gf.asm_generic(mim1, 0, "Norm_sqr(1-t)", -1,
                           "t", False, mf1_free_sym, np.ones(mf1_free_sym.nbdof()))
  assert(L2error < 1e-16);
  print("mf1_free_sym partition of unity error:", L2error)
  L2error = gf.asm_generic(mim1, 0, "Norm_sqr(1-t)", -1,
                           "t", False, mf1_sym_free, np.ones(mf1_sym_free.nbdof()))
  assert(L2error < 1e-16);
  print("mf1_sym_free partition of unity error:", L2error)
  L2error = gf.asm_generic(mim1, 0, "Norm_sqr(1-t)", -1,
                           "t", False, mf1_sym_sym, np.ones(mf1_sym_sym.nbdof()))
  assert(L2error < 1e-16);
  print("mf1_sym_sym partition of unity error:", L2error)
  L2error = gf.asm_generic(mim1, 0, "Norm_sqr(1-t)", -1,
                           "t", False, mf1_periodic, np.ones(mf1_periodic.nbdof()))
  assert(L2error < 1e-16);
  print("mf1_periodic partition of unity error:", L2error)

  L2error = gf.asm_generic(mim2, 0, "Norm_sqr(1-t)", -1,
                           "t", False, mf2_free_free_free_free, np.ones(mf2_free_free_free_free.nbdof()))
  assert(L2error < 1e-16);
  print("mf2_free_free_free_free partition of unity error:", L2error)
  L2error = gf.asm_generic(mim2, 0, "Norm_sqr(1-t)", -1,
                           "t", False, mf2_free_sym_sym_free, np.ones(mf2_free_sym_sym_free.nbdof()))
  assert(L2error < 1e-16);
  print("mf2_free_sym_sym_free partition of unity error:", L2error)
  L2error = gf.asm_generic(mim2, 0, "Norm_sqr(1-t)", -1,
                           "t", False, mf2_sym_sym_sym_sym, np.ones(mf2_sym_sym_sym_sym.nbdof()))
  assert(L2error < 1e-16);
  print("mf2_sym_sym_sym_sym partition of unity error:", L2error)
  L2error = gf.asm_generic(mim2, 0, "Norm_sqr(1-t)", -1,
                           "t", False, mf2_periodic_free_periodic_free, np.ones(mf2_periodic_free_periodic_free.nbdof()))
  assert(L2error < 1e-16);
  print("mf2_periodic_free_periodic_free partition of unity error:", L2error)
  L2error = gf.asm_generic(mim2, 0, "Norm_sqr(1-t)", -1,
                           "t", False, mf2_sym_periodic_free_periodic, np.ones(mf2_sym_periodic_free_periodic.nbdof()))
  assert(L2error < 1e-16);
  print("mf2_sym_periodic_free_periodic partition of unity error:", L2error)

exit(0)
