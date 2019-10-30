.. include:: ../replaces.txt

.. _howtos:

How-tos
=======

Import gmsh mesh
----------------

If we have in the file `quad.geo` a parameterized mesh, as this:

.. literalinclude:: code_samples/quad.geo
   :language: c
   :linenos:

then, when we run::

  $ gmsh -2 quad.geo

the file `quad.msh` is created and contains the encoding of the mesh and its
regions. We can import that file (*quad.msh*) to getfem::

  import getfem as gf

  m = gf.Mesh('import','gmsh','quad.msh')
  print m.regions()

with the second command we can see the *regions ids*. When we import the mesh,
we might be warned with the following::

  Level 3 Warning in getfem_import.cc, line 137:
    All regions must have different number!

this means that the parametrization of the mesh in |gmsh| *.geo file* must
assign a **different** number to each region, the problem exists because in
|gmsh| can coexist, for example, "Physical Surface (200)" and "Physical Line
(200)", as they are different "types of regions" in |gmsh|, that which does
not occur in |gf| since there is only one "type of region".
