.. $Id$

.. include:: ../replaces.txt

.. highlight:: c++

.. _ud-parallel:

MPI Parallelization of |gf|
===========================

Of course, each different problem should require a different
parallelization adapted to its specificities in order to
obtain a good load balancing. You may build your own parallelization
using the mesh regions to parallelize assembly procedures.

Nevertheless, the brick system offers a generic parallelization based on `Open MPI <https://www.open-mpi.org>`_
(communication between processes),
`METIS <http://glaros.dtc.umn.edu/gkhome/metis/metis/overview>`_
(partition of the mesh)
and `MUMPS <http://graal.ens-lyon.fr/MUMPS>`_ (parallel sparse direct solver).
It is available with the compiler option ``-D GETFEM_PARA_LEVEL=2``
and the library itself has to be compiled with the
option ``--enable-paralevel=2`` of the configure script.
Initial MPI parallelization of |gf| has been designed with the help of Nicolas Renon from CALMIP, Toulouse.


When the configure script is run with the option ``--enable-paralevel=2``,
it searches for MPI, METIS and parallel MUMPS libraries.
If the python interface is built, it searches also for MPI4PY library.
In that case, the python interface can
be used to drive the parallel version of getfem (the other interfaces has
not been parallelized for the moment). See demo_parallel_laplacian.py in
the interface/test/python directory.

With the option ``-D GETFEM_PARA_LEVEL=2``, each mesh used is implicitly
partitionned (using METIS) into a
number of regions corresponding to the number of processors and the assembly
procedures are parallelized. This means that the tangent matrix and the
constraint matrix assembled in the model_state variable are distributed.
The choice made (for the moment) is not to distribute the vectors.
So that the right hand side vectors in the model_state variable
are communicated to each processor (the sum of each contribution is made
at the end of the assembly and each processor has the complete vector).
Note that you have to think to the fact that the matrices stored by the
bricks are all distributed.

A model of C++ parallelized program is :file:`tests/elastostatic.cc`.
To run it in parallel you have to launch for instance::

  mpirun -n 4 elastostatic elastostatic.param

For a python interfaced program, the call reads::

  mpirun -n 4 python demo_parallel_laplacian.py

If you do not perform a `make install`, do not forget to first set the shell variable PYTHONPATH to the python-getfem library with for instance::

  export PYTHONPATH=my_getfem_directory/interface/src/python


State of progress of |gf| MPI parallelization
---------------------------------------------

Parallelization of getfem is still considered a "work in progress". A certain number of procedure are still remaining sequential. Of course, a good test to see if the parallelization of your program is correct is to verify that the result of the computation is indeed independent of the number of process.

- Assembly procedures

  Most of assembly procedures (in :file:`getfem/getfem_assembling.h`) have a parameter corresponding to the region in which the assembly is to be computed. They are not parallelized themselves but aimed to be called with a different region in each process to distribute the job. Note that the file :file:`getfem/getfem_config.h` contains a procedures called MPI_SUM_SPARSE_MATRIX allowing to gather the contributions of a distributed sparse matrix.

  The following assembly procedures are implicitly parallelized using the option ``-D GETFEM_PARA_LEVEL=2``:

  * computation of norms (``asm_L2_norm``, ``asm_H1_norm``, ``asm_H2_norm`` ..., in :file:`getfem/getfem_assembling.h`),

  * ``asm_mean_value`` (in :file:`getfem/getfem_assembling.h`),

  * ``error_estimate`` (in :file:`getfem/getfem_error_estimate.h`).

  This means in particular that these functions have to be called on each
  processor.

- Mesh_fem object

  The dof numbering of the getfem::mesh_fem object remains sequential and is executed on each process.  The parallelization is to be done. This could affect the efficiency of the parallelization for very large and/or evoluting meshes.


- Model object and bricks

  The model system is globally parallelized, which mainly means that the
  assembly procedures of standard bricks use a METIS partition of the
  meshes to distribute the assembly. The tangent/stiffness matrices
  remain distibuted and the standard solve call the parallel version
  of MUMPS (which accept distributed matrices).

  For the moment, the procedure ``actualize_sizes()`` of the model
  object remains sequential and is executed on each process.
  The parallelization is to be done.

  Some specificities:

  * The explicit matrix brick: the given matrix is considered to be
    distributed. If it is not, only add it on the master process (otherwise,
    the contribution will be multiplied by the number of processes).

  * The explicit rhs brick: the given vector is not considered to be
    distributed. Only the given vector on the master process is taken into
    account.

  * Constraint brick: The given matrix and rhs are not considered to be
    distributed. Only the given matrix and vector on the master process are
    taken into account.

  * Concerning contact bricks, only integral contact bricks are fully
    parallelized for the moment. Nodal contact bricks work in parallel
    but all the computation is done on the master process.



