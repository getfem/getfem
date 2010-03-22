.. $Id$

.. include:: ../replaces.txt

.. highlightlang:: c++

.. _dp-intro:

Introduction
============

The |gf| project focuses on the development of a generic finite element library. 
The goal is to provide a finite element framework which allows to easily build 
numerical code for the modelisation of system described by partial differential 
equations (p.d.e.). A special attention is paid to the flexibility of the use of 
the library in the sense that the switch from a method offered by the library to 
another is made as easy as possible.

The major point allowing this, compared to traditional finite element codes, is 
the complete separation between the description of p.d.e. models and finite 
element methods. Moreover, a separation is made between integration methods 
(exact or approximated), geometric transformations (linear or not) and finite 
element methods of arbitrary degrees described on a reference element. |gf| can 
be used to build very general finite elements codes, where the finite elements, 
integration methods, dimension of the meshes, are just some parameters that can 
be changed very easily, thus allowing a large spectrum of experimentations. 
Numerous examples are available in the :file:`tests` directory of the 
distribution.

The goal is also to make the addition of new finite element method as simple as 
possible. For standard method, a description of the finite element shape 
functions and the type of connection of degrees of freedom on the reference 
element is sufficient. Extensions are provided for Hermite elements, piecewise 
polynomial, non-polynomial, vectorial elements and XFem. Examples of predefined 
available methods are :math:`P_k` on simplices in arbitrary degrees and 
dimensions, :math:`Q_k` on parallelepipeds, :math:`P_1`, :math:`P_2` with bubble 
functions, Hermite elements, elements with hierarchic basis (for multigrid 
methods for instance), discontinuous :math:`P_k` or :math:`Q_k`, XFem, Argyris, 
HCT, Raviart-Thomas.

The library also includes the usual tools for finite elements such as assembly 
procedures for classical PDEs, interpolation methods, computation of norms, mesh 
operations, boundary conditions, post-processing tools such as extraction of 
slices from a mesh ...

|gf| has no meshing capabilities (apart regular meshes, and a not exploitable 
attempt), hence, in many situations, it is necessary to import meshes. Imports 
formats currently known by getfem are `GiD`_, `Gmsh`_ and `EMC2`_ mesh files. 
However, given a mesh, it is possible to refine it automatically.

The aim of the |gf| project is not to provide a ready to use finite element code 
allowing for instance structural mechanics computations with a graphic interface. 
It is basically a library allowing the build of C++ finite element codes. 
However, the matlab and python interfaces allows to easily build application 
coupling the definition of the problem, the finite element methods selection and 
the graphical post-processing.

The future of the project is to continue to develop the finite element framework, 
focusing on the following points.

* Background consolidation of the existing modules (with a reflection on the 
  optimal representation of meshes, degrees of freedom, finite element methods 
  ...).
* Developpement of innovating methods.
* Reflection on the optimal way to represent complex p.d.e. models with the 
  maximum of flexibility and reusability. The brick system is a first step in 
  this direction.

The vocation of |gf| is to remain a free open source project. The advantage given 
by the fact to be an open source project is that by proposing a free use, one 
profits from the experiments of the users who by their tests and the difficulties 
or bug which they meet make progress the robustness of the algorithms. One also 
profits from the possible contributions of the users who can find interest to 
develop new functinalities within the proposed framework. That allows 
constructive exchanges which clarify the weak points and the strong points of the 
project.

Figure :ref:`dp-fig-diagram` describes the diagram of the different modules of 
the |gf| library. The current state and perspective for each module is described 
in section :ref:`dp-libdesc`.

.. _dp-fig-diagram:
.. figure:: images/diagram.png
   :align: center
   :scale: 80

   Diagram of |gf| library

.. include:: ../license.txt
