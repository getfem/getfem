.. $Id$

.. include:: ../replaces.txt

.. highlightlang:: c++

.. _ud-linalg:

Linear algebra procedures
=========================

The linear algebra library used by |gf| is |gmm| which is now a separate library.
Please see the `GMM++ user documentation
<http://getfem.org/gmm.html>`_.

Note that |gf| includes (since release 1.7) its own version of |sLU| 3.0 (see
`SuperLU web site <http://crd.lbl.gov/~xiaoye/SuperLU>`_) hence a direct sparse
solver is available out of the box. Note that an option of the ``./configure``
file allows to disable the included version of |sLU| in order to use a
pre-installed version.

A small interface to |mumps| is also provided (see `MUMPS web1
<http://graal.ens-lyon.fr/MUMPS>`_ or `MUMPS web2
<http://www.enseeiht.fr/apo/MUMPS>`_). See the file
:file:`gmm/gmm_MUMPS_interface.h`. In order to use |mumps|, you have to indicates
some options to the configure shell::

  --with-mumps-include-dir=" -I /path/to/MUMPS/include "
  --with-mumps=" F90 libraries and libs of MUMPS to be linked "

alternatively, the option ``--enable-mumps`` will search for an installed MUMPS library. Note that if both the sequential and the parallel version is installed on your system (especially on Debian and Ubuntu), the default version will be the parallel one. To select the sequential one it is necessary to add the option ``--with-mumps="-lsmumps_seq -ldmumps_seq -lcmumps_seq -lzmumps_seq"``.

For instance if you want to use the sequential version of |mumps| with double and
complex double::

  --with-mumps-include-dir=" -I /path/to/MUMPS/include "
  --with-mumps=" ...F90libs...  -L /path/to/MUMPS/lib -ldmumps -lzmumps -lpord
              -L /path/to/MUMPS/libseq -lmpiseq "

where ``...F90libs...`` are the libraries of the fortran compiler used to compile
|mumps| (these are highly dependant on the fortran 90 compiler used, the
``./configure`` script should detect the options relative to the default fortran
90 compiler on your machine and display it -- for example, with the intel
``ifort`` compiler, it is ``-L/opt/icc8.0/lib -lifport -lifcoremt -limf -lm
-lcxa -lunwind -lpthread``)
