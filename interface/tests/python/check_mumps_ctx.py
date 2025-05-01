#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Python GetFEM interface
#
# Copyright (C) 2025-2025 Konstantinos Poulios.
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
"""  Test of the mumps context object.

  This program is used to check that Python-GetFEM interface, and more
  generally GetFEM are working. It focuses on testing the interface of
  the MUMPS linear system solver in GetFEM.

  $Id$
"""

import numpy as np
import getfem as gf

nprocs = 1
rank = 0
if gf.util_mpi_parallelism_level() >= 2:
  from mpi4py import MPI
  nprocs = MPI.COMM_WORLD.Get_size()
  rank = MPI.COMM_WORLD.Get_rank()

ctx1 = gf.MumpsContext("symmetric")
ctx2 = gf.MumpsContext("unsymmetric")
ctx1.set_ICNTL(4, 0) # silence MUMPS output
ctx2.set_ICNTL(4, 0) # silence MUMPS output

K1 = gf.Spmat('empty', 3, 3)
K2 = gf.Spmat('empty', 3, 3)
rhs = np.ones(3)

K1.add(0, 0, 1)
K1.add(1, 1, 2)
K1.add(2, 2, 2)
K1.add(1, 2, 0.2)
K1.add(2, 1, 0.2)

K2.add(0, 0, 1)
K2.add(1, 1, 2)
K2.add(2, 2, 3)
K2.add(1, 2, 0.2)
K2.add(2, 1, -0.2)

ctx1.set_matrix(K1)
ctx1.set_vector(rhs)

ctx1.analyze()
ctx1.factorize()
x1 = ctx1.solve()
#print(f"rank{rank}: x1 =", x1)

ctx1.set_vector(2*rhs)
x1_times_two = ctx1.solve()
#print(f"rank{rank}: 2*x1 =", x1_times_two)

assert np.allclose(x1, x1_times_two/2)

ctx1.set_vector(np.column_stack((rhs,2*rhs)))
x1,x1_times_two = np.unstack(ctx1.solve(), axis=1)
#print(f"rank{rank}: x1 =", x1, "2*x1 =", x1_times_two)


ctx2.set_matrix(K2)
ctx2.set_vector(rhs)

ctx2.analyze()
ctx2.factorize()
x2 = ctx2.solve()
#print(f"rank{rank}: x2 =", x2)

ctx2.set_vector(2*rhs)
x2_times_two = ctx2.solve()
#print(f"rank{rank}: 2*x2 =", x2_times_two)

assert np.allclose(x2, x2_times_two/2)

ctx3 = gf.MumpsContext("unsymmetric","complex")
ctx3.set_ICNTL(4, 0)  # silence MUMPS output

K3 = gf.Spmat('empty', 3, 3)
K3.to_complex()
rhs_cplx = np.array([1+1j, 2+1j, 1+2j])

K3.add(0, 0, 1+0j)
K3.add(1, 1, 2+0j)
K3.add(2, 2, 2+0j)
K3.add(1, 2, 0.2+0j)
K3.add(2, 1, 0.2+0j)

ctx3.set_matrix(K3)
ctx3.set_vector(rhs_cplx)

ctx3.analyze()
ctx3.factorize()
x3 = ctx3.solve()
#print(f"rank{rank}: x3 =", x3)
#print(f"rank{rank}: rhs_cplx =", rhs_cplx)
#print(f"rank{rank}: K3.mult(x3) =", K3.mult(x3))
assert np.allclose(rhs_cplx, K3.mult(x3))

ctx3.set_vector(2*rhs_cplx)
x3_times_two = ctx3.solve()
#print(f"rank{rank}: 2*x3 =", x3_times_two)

assert np.allclose(x3, x3_times_two/2)


if gf.util_mpi_parallelism_level() >= 2:
  ctx1 = gf.MumpsContext("symmetric")
  ctx2 = gf.MumpsContext("unsymmetric")
  ctx1.set_ICNTL(4, 0) # silence MUMPS output
  ctx2.set_ICNTL(4, 0) # silence MUMPS output

  K1 = gf.Spmat('empty', 3, 3)
  K2 = gf.Spmat('empty', 3, 3)
  rhs = np.zeros(3)
  K1.add(2, 2, 2./nprocs)
  K2.add(1, 1, 2./nprocs)
  if rank == 0:
    rhs = np.ones(3)
    K1.add(0, 0, 1)
    K1.add(1, 1, 2)

    K2.add(0, 0, 1)
    K2.add(1, 2, 0.2)
  if rank == nprocs-1:
    K1.add(1, 2, 0.2)
    K1.add(2, 1, 0.2)

    K2.add(2, 2, 3)
    K2.add(2, 1, -0.2)

  ctx1.set_distributed_matrix(K1)
  ctx1.set_vector(rhs)

  x1 = ctx1.analyze_factorize_and_solve()
  #print(f"rank{rank}: x1 =", x1)

  ctx2.set_distributed_matrix(K2)
  ctx2.set_vector(rhs)

  x2 = ctx2.analyze_factorize_and_solve()
  #print(f"rank{rank}: x2 =", x2)

