.. $Id$

.. include:: ../replaces.txt

.. highlight:: c++

.. _ud-intro:

Introduction
============

The |gf| project focuses on the development of a generic and
efficient |c++| library for finite element methods elementary
computations. The goal is to provide a library allowing the
computation of any elementary matrix (even for mixed finite element
methods) on the largest class of methods and elements, and for
arbitrary dimension (i.e. not only 2D and 3D problems).

It offers a complete separation between integration methods (exact or
approximated), geometric transformations (linear or not) and finite
element methods of arbitrary degrees. It can really relieve a more
integrated finite element code of technical difficulties of
elementary computations.

Examples of available finite element method are : Pk on simplices in
arbitrary degrees and dimensions, Qk on parallelepipeds, P1, P2 with
bubble functions, Hermite elements, elements with hierarchic basis
(for multigrid methods for instance), discontinuous Pk or Qk, XFem,
Argyris, HCT, Raviart-Thomas, etc.

The addition of a new finite element method is straightforward. Its
description on the reference element must be provided (in most of the
cases, this is the description of the basis functions, and nothing
more). Extensions are provided for Hermite elements, piecewise
polynomial, non-polynomial and vectorial elements, XFem.

The library also includes the usual tools for finite elements such as
assembly procedures for classical PDEs, interpolation methods,
computation of norms, mesh operations, boundary conditions,
post-processing tools such as extraction of slices from a mesh, etc.

|gf| can be used to build very general finite elements codes, where
the finite elements, integration methods, dimension of the meshes,
are just some parameters that can be changed very easily, thus
allowing a large spectrum of experimentations. Numerous examples are
available in the ``tests`` directory of the distribution.

|gf| has only a (very) experimental meshing procedure (and produces regular meshes), hence it is generally
necessary to import meshes. Imports formats currently known by |gf|
are |gid|, |gmsh| and *emc2* mesh files. However, given a mesh, it
is possible to refine it automatically.

.. include:: ../license.txt
