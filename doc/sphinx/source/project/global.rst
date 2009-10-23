.. $Id$

.. include:: ../replaces.txt

.. highlightlang:: c++

.. _dp-global:

Global perspectives of structuration, consolidation and growth
==============================================================

intro to the main modifications to be done ...

Modifications to be done are of three kind:

* Background consolidation of the existing modules (with a reflection on the
  optimal representation of meshes, degrees of freedom, finite element methods,
  etc.).

* Developpement of innovating methods.

* Reflection on the optimal way to represent complex p.d.e. models with the
  maximum of flexibility and reusability. The brick system is a first step in
  this direction. It should be replaced soon by a more elaborated system.


Namespace changes
-----------------

After the elimination of the small namespaces ``linkmsg`` and ``ftool`` in
release 3.0, it remains now four namespaces in the |gf| project.

* ``gmm`` (Generic Matrix Methods) : for the linear algebra procedures.

* ``dal`` (Dynamic Array Library) : some basic algorithms including the
  definition of some containers (``dal::dynamic_array``, ``dal::dynamic_tas``,
  ``dal::tree_sorted_array``, ``dal::bit_vector``).

* ``bgeot`` (Basic GEOmetric Tool) : some basic algorithms including the
  definition of geometric objects (convex structure, convex, convex of reference,
  basic mesh).

* ``getfem`` : the main namespace of |gf|.

It is clear that the separation into these remaining four namespaces is mainly
historical. The separate ``gmm`` namespace for |gmm| is clearly justified. The
contour of nemaspaces ``dal`` and ``bgeot`` is more vague. Historically, those
two namespaces had their own justifications.

In the very begining of |gf| (the first files was written in 1995) the S.T.L. was
not available and the containers defined in the ``dal`` namespace was used
everywhere. Now, in |gf|, the S.T.L. containers are mainly used. The remaining
uses of ``dal`` containers are eather historical or due to the specificities of
these containers. It is however clear that this is not the aim of the |gf|
project to developp new container concept. So, the use of the ``dal`` containers
has to be as much as possible reduced.

Now, concerning ``bgeot``, it was containing some other geometrical object at the
begining and was originally designed to be a self-consistent library of geometric
concepts. It slowly derived to be like it is now, a collection of algorithms and
object definition more or less related to geometry (rtree, kdtree, ftool,
polynomials ...).

The conclusion of this is that ``dal`` and ``bgeot`` namespaces can be
advantageously merged to the ``getfem`` namespace, reducing to the minimum the
use of the ``dal`` containers. This should be done preserving the backward
compatibility. An intermediary study would be to see if the ``dal`` cannot be
directly derived from S.T.L. containers preserving the used specificities.


Basic types used
----------------

Basic type of integer, real ... used. to be done.
