.. $Id: model.rst 3655 2010-07-17 20:42:08Z renard $

.. include:: ../replaces.txt

.. highlightlang:: c++

.. index:: models, model bricks, Nitsche's method

.. _ud-model-Nitsche:


Nitsche's method for dirichlet and contact boundary conditions
--------------------------------------------------------------

Nitsche's method is a very attractive methods since it allows to take into
account A dirichlet type boundary condition or a contact with friction boundary condition in a weak way without the use of a Lagrange multiplier. \gf provides a generic implementation of Nitche's method. The advantage of Nitsche's method, which is to tranform a Dirichlet boundary condition into weak terms similarly as a Neumann boundary condition, is paid by the fact that the implementation is equation dependent. This method needs the use of an approximation of the corresponding Neumann term. Thus, in order to add a boundary condition with Nitsche's method on a variable of a model, the corresponding brick as to have access to an approximation of the Neumann term of all the partial differential terms applied to this variables. In the following, considering a variable :math:`u`, we will denote by

.. math::
  G(u)

the sum of all the Neumann terms on its variables. For instance, if ther is a Laplace term (:math:`\Delta u`) applied on :math:`u`, the Neumann term will be :math:`G(u) = \Frac{\partial u}{\partial n}` where :math:`n` is the outward unit normal on the considered boundary. If :math:`u` represents the displacement of a deformable body, the Neumann term will be :math:`G(u) = \sigma(u)n`, where :math:`\sigma(u)` is the stress tensor depending on the consitutive law. Of course, :math:`G(u)` may (and generally will) depend on some parameters or some other variables.

In order to propose a generic implementation in which the brick proposing Nitsche's method are not dependent on the partial differential terms applied to the concerned variables, each brick adding a partial differential term is asked to give the expression of the corresponding Neumann term. Of course, it makes the building of a brick a little bit more complicated. So, this mechanism is not mandatory to build a new brick, but of course, it is not possible to use a Nitsche's method brick with a brick which do not propose the corresponding Neumann term. An internal mechanism of \gf controls this.
 





Rule : contrarily to the other ... it has to be added after the P.D.E terms on the variables ... -> Warning à mettre






   


