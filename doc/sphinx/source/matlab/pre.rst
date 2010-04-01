.. $Id$

.. include:: ../replaces.txt

.. highlightlang:: matlab

.. _mlab-pre:

Preliminary
===========

This is just a short summary of the terms employed in this manual. If you are not
familiar with finite elements, this should be useful (but in any case, you should
definitively read the :ref:`dp`).

The :envvar:`mesh` is composed of :envvar:`convexes`. What we call convexes can be
simple line segments, prisms, tetrahedrons, curved triangles, of even something
which is not convex (in the geometrical sense). They all have an associated
:envvar:`reference convex`: for segments, this will be the :math:`[0,1]` segment,
for triangles this will be the canonical triangle :math:`(0,0)-(0,1)-(1,0)`, etc.
All convexes of the mesh are constructed from the reference convex through a
:envvar:`geometric transformation`. In simple cases (when the convexes are
simplices for example), this transformation will be linear (hence it is easily
inverted, which can be a great advantage). In order to define the geometric
transformation, one defines :envvar:`geometrical nodes` on the reference convex.
The geometrical transformation maps these nodes to the :envvar:`mesh nodes`.

On the mesh, one defines a set of basis functions: the :envvar:`FEM`. A FEM is
associated at each convex. The basis functions are also attached to some
geometrical points (which can be arbitrarily chosen). These points are similar to
the mesh nodes, but **they don't have to be the same** (this only happens on very
simple cases, such as a classical :math:`P_1` fem on a triangular mesh). The set
of all basis functions on the mesh forms the basis of a vector space, on which the
PDE will be solved. These basis functions (and their associated geometrical point)
are the :envvar:`degrees of freedom` (contracted to :envvar:`dof`). The FEM is
said to be :envvar:`Lagrangian` when each of its basis functions is equal to one
at its attached geometrical point, and is null at the geometrical points of others
basis functions. This is an important property as it is very easy to
:envvar:`interpolate` an arbitrary function on the finite elements space.

The finite elements method involves evaluation of integrals of these basis
functions (or product of basis functions etc.) on convexes (and faces of
convexes). In simple cases (polynomial basis functions and linear geometrical
transformation), one can evaluate analytically these integrals. In other cases,
one has to approximate it using :envvar:`quadrature formulas`. Hence, at each
convex is attached an :envvar:`integration method` along with the FEM. If you have
to use an approximate integration method, always choose carefully its order (i.e.
highest degree of the polynomials who are exactly integrated with the method): the
degree of the FEM, of the polynomial degree of the geometrical transformation, and
the nature of the elementary matrix have to be taken into account. If you are
unsure about the appropriate degree, always prefer a high order integration method
(which will slow down the assembly) to a low order one which will produce a
useless linear-system.

The process of construction of a global linear system from integrals of basis
functions on each convex is the :envvar:`assembly`.

A mesh, with a set of FEM attached to its convexes is called a :envvar:`mesh_fem`
object in |gf|.

A mesh, with a set of integration methods attached to its convexes is called a
:envvar:`mesh_im` object in |gf|.

A |mf| can be used to approximate scalar fields (heat, pression, ...), or vector
fields (displacement, electric field, ...). A |mim| will be used to perform
numerical integrations on these fields. Most of the finite elements implemented in
|gf| are scalar (however, :math:`TR_0` and edges elements are also available). Of
course, these scalar FEMs can be used to approximate each component of a vector
field. This is done by setting the :math:`Qdim` of the |mf| to the dimension of
the vector field (i.e. :math:`Qdim=1` :math:`\Rightarrow` scalar field,
:math:`Qdim=2` :math:`\Rightarrow` 2D vector field etc.).

When solving a PDE, one often has to use more than one FEM. The most important one
will be of course the one on which is defined the solution of the PDE. But most
PDEs involve various coefficients, for example:

.. math::

   \nabla\cdot(\lambda(x)\nabla u) = f(x).

Hence one has to define a FEM for the main unknown :math:`u`, but also for the
data :math:`\lambda(x)` and :math:`f(x)` if they are not constant. In order to
interpolate easily these coefficients in their finite element space, one often
choose a Lagrangian FEM.

The convexes, mesh nodes, and dof are all numbered. We sometimes refer to the
number associated to a convex as its :envvar:`convex id` (contracted to
:envvar:`cvid`). Mesh node numbers are also called :envvar:`point id` (contracted
to :envvar:`pid`). Faces of convexes do not have a global numbering, but only a
local number in each convex. Hence functions which need or return a list of faces
will always use a two-rows matrix, the first one containing convex ids, and the
second one containing local face number.

While the dof are always numbered consecutively, **this is not always the case for
point ids and convex ids**, especially if you have removed points or convexes from
the mesh. To ensure that they form a continuous sequence (starting from 1), you
have to call::

  >> gf_mesh_set(m,'optimize structure')
