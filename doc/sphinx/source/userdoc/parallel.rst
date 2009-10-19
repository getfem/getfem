.. $Id$

.. include:: ../replaces.txt

.. highlightlang:: c++

.. _ud-parallel:

Parallelization of |gf|
=======================

Of course, each different problem should require a different parallelization
adapted to its specificities. You may build your own parallelization using the
mesh regions to parallelize assembly procedures.

Nevertheless, the brick system offers a generic parallelization based on MPI
(communication between processes), `METIS
<http://glaros.dtc.umn.edu/gkhome/metis/metis/overview>`_ (partition of the mesh)
and `MUMPS <http://graal.ens-lyon.fr/MUMPS>`_ (parallel sparse direct solver). One
has to compile |gf| with the option ``-D GETFEM_PARA_LEVEL=2`` to use it. With
this option, each mesh used is implicitely partitionned (using METIS) into a
number of regions corresponding to the number of processors and the assembly
procedures are parallelized. This means that the tangent matrix and the constraint
matrix assembled in the model_state variable are distributed. The choice made (for
the moment) is not to distribute the vectors. So that the right hand side vectors
in the model_state variable are communicated to each processor (the sum of each
contribution is made at the end of the assembly and each processor has the
complete vector). Note that you have to think to the fact that the matrices stored
by the bricks are all distributed.

Concerning the constraints, it is preferable to avoid the
``getfem::ELIMINATED_CONSTRAINTS`` option for a better parallelization (i.e. not
to use the constraint matrix).

A model of parallelized program is :file:`tests/elastostatic.cc`.

The following functions are also implicitely parallelized using the option ``-D
GETFEM_PARA_LEVEL=2``:

* computation of norms (``asm_L2_norm``, ``asm_H1_norm``, ``asm_H2_norm`` ..., in
  :file:`getfem/getfem_assembling.h`),

* ``asm_mean_value`` (in :file:`getfem/getfem_assembling.h`),

* ``error_estimate`` (in :file:`getfem/getfem_error_estimate.h`).

This means that these functions have to be called on each processor.

Parallelization of getfem is still considered a "work in progress"...
