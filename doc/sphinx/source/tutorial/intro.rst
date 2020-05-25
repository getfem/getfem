.. $Id: install.rst 4738 2014-07-27 12:25:54Z renard $

.. include:: ../replaces.txt

.. highlight:: c++

.. _ud-intro-tut:

Introduction
============

This tutorial is intended to present the main aspects of the |gf| library  on commented simple examples. |gf| allows numerical modeling of complex problems including coupled ones. It is relatively simple to add arbitrary coupling terms thanks to a weak form description language. Expressions in weak form are then compiled in the form of a list of optimized instructions that are executed during assembly.

The main presented functionalities are the following

  - Approximation of multiphysics problems and the use of the generic assembly.
  - The use of fixed-size variables.
  - Contact problems.
  - Fictitious domain functionalities and transformations allowing inter-mesh or inter-region assembly.
  - Continuation / bifurcation problems.

This tutorial does not of course takes the place of the various documentations. The user documentation (that of the C++ library) is the reference one concerning the description of implemented methods. Interfaces documentations describe the use of the different functionalities from interfaces but do not repeat the description of methods. 'Developpers's guide' documentation describes more internal concerning the organization of the library.

|Gf| is a free collaborative project hosted by the `Savannah <http://savannah.gnu.org>`_ site (see https://savannah.nongnu.org/projects/getfem). New contributors are welcome to all aspects of the project.



C++, Python, Scilab or Matlab ?
-------------------------------



|gf| is basically a library written in C++. However, it has an interface that is available in Python, Scilab and Matlab  versions and allows to use all the advanced features. It is recommended to start by using the interface with scripting languages, which makes programming easier. The three versions of the interface differ only by small syntax elements, except for the graphics post-processing (Scilab and Matlab interfaces offer the possibility of integrated post-processing while it is necessary to use dedicated software such as Paraview, Mayavi or gmsh when using the Python interface or directly the library in C++). The Python interface is compiled by default along with the C++ library.

Use |gf| interface is a good strategy even for complex applications, the performance is comparable to the direct use of the C++ library and it is more flexible to use. However, only the Python interface allows for the parallelization (with MPI). The possible addition of functionality to the interface being a relatively simple operation.

The first example of the tutorial (thermo-electric coupling elastico) allows to see the difference in the use of the C++ library and one of the versions of the interface.


Where are demo files ?
----------------------

A large number of demonstration programs can be found

* for the C++ examples in directory::

        tests/

* for the Python interface in the directory::

        interface/tests/python

* for the Scilab interface in the directory::

        interface/scr/scilab/demos

* and for the Matlab interface in the directory::

        interface/tests/matlab

