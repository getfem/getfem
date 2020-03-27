.. $Id$

.. include:: ../replaces.txt

.. highlightlang:: python

.. _ud-hho:

Tools for HHO (Hybrid High-Order) methods
=========================================



HHO method are hybrid methods in the sense that they have both degrees of freedom located on the element of a mesh and on the faces of the elements which represent separated approximations. HHO method are primal methods in the sense that both the degree of freedom in the element and on the faces represent the main unknown of the problem (no lagrange multipliers is introduced). The interest of these methods, first developped in  [Di-Er2015]_, [Di-Er2017]_ is their accuracy and their great robustness, in particular with respect to the element shapes and their locking-free properties. Moreover, they can be extended without difficulty to the approximation of nonlinear problems (see [AB-ER-PI2018]_ for hyper-elasticity, [AB-ER-PI2019]_ for plasticity and [ca-ch-er2019]_ for contact problems).

HHO methods can be applied to arbitrary shape elements. However, the implementation in |gf| is for the moment limited to standard elements : simplices, quadrilaterals, hexahedrons, ... Moreover this implementation is still experimental and not pretending to optimality. For the moment, there is no tool to make an automatic condensation of internal dofs.

HHO elements
------------

HHO elements are composite ones having a polynomial approximation space for the interior of the element and a polynomial approximation for each face of the element. Moreover, this is a discontinous approximation, in the sens that no continuity is prescribed between the approximation inside the element and the approximation on the faces, neither than between the approximations on two different faces of the element. However, when two neighbor elements share a face, the approximation on this face is shared by the two elements. |gf| provide a specific method simply called ``FEM_HHO(fem_int, fem_face1, fem_face2, ...)`` which allows to build an hybrid method from standard finite element spaces. For instance, on a triangle, a possible HHO method can be obtained with::

  getfem::pfem pf = getfem::fem_descriptor("HHO(FEM_SIMPLEX_IPK(2,2), FEM_SIMPLEX_CIPK(1,2))");

The first argument to ``FEM_HHO(...)`` is the fem for the interior of the element. It has to be a discontinuous FEM. The method ``FEM_SIMPLEX_IPK(2,2)`` is a discontinous method having its degrees of freedom in the strict interior of the element, which ensure that no dof identification will be done. The second argument is the fem for the faces (if only one method is given, it will be applied to all faces, but it is also possible to give a different method for each face). Their is no verification on the fact that the given method are of discontinuous type (In fact, a method like ``FEM_HHO(FEM_PK(2,2), FEM_PK(1,2))`` will have no difference with ``FEM_PK(2,2)`` since the degree of freedom on the faces will be identified with the interior ones).

For the moment, the fursnished element for interior and faces are
- ``FEM_SIMPLEX_IPK(n,k)`` : interior PK element of degree k for the simplices in dimension n (equivalent to ``FEM_PK_DISCONTINUOUS(n,k,0.1)``).
- ``FEM_QUAD_IPK(n,k)`` : interior PK element of degree k for the quadrilaterals in dimension n.
- ``FEM_PRISM_IPK(n,k)`` : interior PK element of degree k for the prisms in dimension n.
- ``FEM_SIMPLEX_CIPK(n,k)`` : interior PK element on simplices which is additionnaly connectable. Designed to be use on HHO element face. 
- ``FEM_QUAD_CIPK(k)`` : interior PK element on a quadrilateral which is additionnaly connectable. Designed to be use on HHO element face. 

Reconstruction operators
------------------------

For a variable ``u``, we will note :math:`u_{T}` its value in the interior of the element :math:`T` and :math:`u_{\partial T}` its value on the boundary of :math:`T` (corresponding to the two different approximations). The reconstruction operators are implemeted in |gf| as elementary transformations, as described in the section :ref:`ud-gasm-high-elem-trans`.

Reconstructed gradient
++++++++++++++++++++++

The first reconstruction operator is the reconstructed gradient. Given a certain polynomial space :math:`V_G`, the reconstructed gradient :math:`G(u)` will be the solution to the local problem

.. math::
   \int_T G(u):\tau dx = \int_T \nabla u_T : \tau dx + \int_{\partial T} (u_{\partial T} - u_{T}).(\tau n_T) d\Gamma, ~~~ \forall \tau \in V_G

where :math:`n_T` is the outward unit normal to  :math:`T` on  :math:`\partial T`. Note that the space :math:`V` is a vector-valued one if ``u`` is a scalar field variable (in that case, :math:`G(u):\tau` reduces to :math:`G(u).\tau`) and a matrix-valued one if ``u`` is a vector field variable.

In order to be used, the elementary transformation corresponding to this operator has first to be added to the model by the command:: 

  add_HHO_reconstructed_gradient(model, transname);

where ``transname`` is an arbitrary name which will designate the transformation in GWFL (the generic weak form language). Then, it will be possible to refer to the reconstructed gradient of a variable ``u`` into GWFL as ``Elementary_transformation(u, HHO_grad, Gu)``, if ``transname="HHO_grad"``. The third parameter of the transformation ``Gu`` should be a fem variable or a data of the model. This variable will not be used on itself but will determine the finite element space of the reconstruction (the space :math:`V_G`).

This is an example of use with the Python interface for a two-dimensional triangule mesh ``m`` ::

  mfu   = gf.MeshFem(m, 1)
  mfgu  = gf.MeshFem(m, N)
  mfu.set_fem(gf.Fem('FEM_HHO(FEM_SIMPLEX_IPK(2,2),FEM_SIMPLEX_CIPK(1,2))'))
  mfgu.set_fem(gf.Fem('FEM_PK(2,2)'))

  md = gf.Model('real')
  md.add_fem_variable('u', mfu)
  md.add_fem_data('Gu', mfgu)

  md.add_HHO_reconstructed_gradient('HHO_Grad')
  md.add_macro('HHO_Grad_u', 'Elementary_transformation(u, HHO_Grad, Gu)')
  md.add_macro('HHO_Grad_Test_u', 'Elementary_transformation(Test_u, HHO_Grad, Gu)')

The macro definitions allowing to use the gradient of the variable inside weak formulations as usual. For instance, the addition of a weak term for the Laplace equation can then be simply written::

  md.add_linear_term(mim, 'HHO_Grad_u.HHO_Grad_Test_u')

Two complete examples of use are given in the test programs :file:`interface/tests/demo_laplacian_HHO.py` and :file:`interface/tests/demo_elasticity_HHO.py`.

Reconstructed symmetrized gradient
++++++++++++++++++++++++++++++++++

The symmetrized gradient is only for vector field variables and additionally when the vector field dimension is the same as the domain dimension. This is usually the case for instance for elasticity problems. With the same notation as in the previous section, the reconstructed gradient :math:`G^s(u)` will be the solution to the local problem

.. math::
   \int_T G^s(u):\tau dx = \int_T \nabla^s u_T : \tau dx + \int_{\partial T} (u_{\partial T} - u_{T}).(\tau^s n_T) d\Gamma, ~~~ \forall \tau \in V_G

where :math:`\nabla^s u_T = (\nabla u_T + (\nabla u_T)^T)/2` and :math:`\tau^s = (\tau + \tau^T)/2`.

The elementary transformation corresponding to this operator can be added to the model by the command:: 

  add_HHO_reconstructed_symmetrized_gradient(model, transname);

and then be used into GWFL as ``Elementary_transformation(u, HHO_sym_grad, Gu)``, if ``transname="HHO_sym_grad"``, with ``Gu`` still determining the reconstruction space.

Reconstructed variable
++++++++++++++++++++++

A recontruction of higher order can be done using both the approximation on the interior and the approximation on the faces. The recontructed variable :math:`D(u)` will be the solution to the local Neumann problem on a chosen space :math:`V_D`

.. math::
   \int_T \nabla D(u). \nabla v dx = \int_T \nabla u_T . \nabla v dx + \int_{\partial T} (u_{\partial T} - u_{T}).(\nabla v n_T) d\Gamma, ~~~ \forall v \in V_D

with the additional constraint

.. math::

   \int_T D(u) dx = \int_T u_T dx

The corresponding elementary transformation can be added to the model by the command:: 

  add_HHO_reconstructed_value(model, transname);

and used into GWFL as ``Elementary_transformation(u, HHO_val, ud)``, if ``transname="HHO_val"``, with ``ud`` determining the reconstruction space.

Reconstructed variable with symmetrized gradient
++++++++++++++++++++++++++++++++++++++++++++++++

A variant of the recontruction of a variable is the one using a symmetrized gradient. It can be used only for vector field variables and additionally when the vector field dimension is the same as the domain dimension. The recontructed variable :math:`D(u)` will be the solution to the local Neumann problem on a chosen space :math:`V_D`

.. math::
   \int_T \nabla^s D(u). \nabla^s v dx = \int_T \nabla^s u_T . \nabla^s v dx + \int_{\partial T} (u_{\partial T} - u_{T}).(\nabla^s v n_T) d\Gamma, ~~~ \forall v \in V_D

with the additional constraints

.. math::

  & \int_T D(u) dx = \int_T u_T dx
   
   &\int_T \mbox{Skew}(\nabla D(u)) dx = \int_{\partial T} (n_T \otimes u_{\partial T} - u_{\partial T} \otimes n_T)/2 d\Gamma

where :math:`\mbox{Skew}(\nabla D(u)) = (\nabla D(u) - (\nabla D(u))^T)/2`.

The corresponding elementary transformation can be added to the model by the command:: 

  add_HHO_reconstructed_value(model, transname);

and used into GWFL as ``Elementary_transformation(u, HHO_val, ud)``, if ``transname="HHO_val"``, with ``ud`` determining the reconstruction space.


Stabilization operators
-----------------------

The stabilization operators is an operator that measure in a sense the discontinuity of the approximation. A stabilization is obtained by a penalization term using this operator. The stabilization operator :math:`S(u)` is defined on the boundary space :math:`V_{\partial T}` of the element, with the formula

.. math::
   S(u) = \Pi_{\partial T}(u_{\partial T} - D(u) - \Pi_{T}(u_T - D(u)))

where :math:`D(u)` is the reconstruction operator on a polynomial space one degree higher that the finite element space used for the variable, :math:`\Pi_{\partial T}` is the :math:`L^2` projection onto the space of the face approximations and  :math:`\Pi_{T}` the :math:`L^2` projection onto the space of the interior of the element.

For vector field variables having the same dimension as the domain, there exists also a stabilization operator using the symmetrized gradient, which is defined by

.. math::
   S^s(u) = \Pi_{\partial T}(u_{\partial T} - D^s(u) - \Pi_{T}(u_T - D^s(u)))

The corresponding elementary transformations can be added to the model by the two commands::
  
  add_HHO_stabilization(model, transname);
  add_HHO_symmetrized_stabilization(model, transname);

and used into GWFL as ``Elementary_transformation(u, HHO_stab)``, if ``transname="HHO_stab"``. A third argument is optional to specify the target (HHO) space (the default is one of the variable itself). An example of use is also given in the test programs :file:`interface/tests/demo_laplacian_HHO.py` and :file:`interface/tests/demo_elasticity_HHO.py`.
