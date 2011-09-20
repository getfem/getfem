.. $Id: model.rst 3655 2010-07-17 20:42:08Z renard $

.. include:: ../replaces.txt

.. highlightlang:: c++

.. index:: models, model bricks

.. _ud-model-contact-friction:



Contact with Coulomb friction brick
-----------------------------------

The aim of this brick is to take into account a contact condition with or without friction of an elastic structure on a rigid foundation or between two elastic structures. This brick is restricted to small deformation approximation of contact.

Approximation of contact
++++++++++++++++++++++++

For small deformation problems submitted
a simple (compared to large deformation !) expression of the contact with friction condition is usually used where the tangential displacement do not influence the normal one. This is an approximation in the sense that if an obstacle is not perfectly flat, the tengential displacement of course influence the point where the contact holds. This will not be the case in small deformation where the contact condition can be considered to be described on the reference configuration.

There is mainly to largely used discretizations of the contact with friction condition in this framework: a direct nodal contact condition (usually prescribed on the displacement finite element nodes) or a weak nodal contact condition (usually prescribed on the multiplier finite element nodes). The two discretization leads to similar system. However, the interpretation of quantities is not the same.

More details can be found for instance in [KI-OD1988]_ and [KH-PO-RE2006]_

Direct nodal contact condition
++++++++++++++++++++++++++++++

A nodal contact condition consists in a certain number of contact nodes :math:`a_i`, :math:`i=1..N_c` on which a contact with (or without) friction condition is applied. The contact condition reads

.. math::

  u_N(a_i)-\text{gap}_i \le 0, ~~ \lambda_N^i \le 0,  ~~ (u_N(a_i)-\text{gap}_i) \lambda_N^i = 0,

where :math:`\lambda_N^i` is the equivalent nodal contact force on :math:`a_i` and :math:`u_N(a_i)` is the normal relative displacement between the elastic solid and an obstacle or between two elastic solids. The term :math:`\text{gap}_i` represents the normal gap between the two solids in the reference configuration. The friction condition reads

.. math::

  \|\lambda_T^i\| \le -{\mathscr F} \lambda_N^i,

  \lambda_T^i = {\mathscr F} \lambda_N^i \frac{\dot{u}_T}{\|\dot{u}_T\|} ~~~ \text{ when } \dot{u}_T \ne 0,

where :math:`\dot{u}_T` is the relative slip velocity, :math:`{\mathscr F}` is the friction coefficient and :math:`\lambda_T^i` the equivalent nodal friction force on :math:`a_i`. The friction condition can be summarized by the inclusion

.. math::

  \lambda_T^i \in {\mathscr F} \lambda_N^i \text{Dir}(\dot{u}_T),

where :math:`\text{Dir}(\dot{u}_T)` is the multivalued map being the sub-differential of :math:`x \mapsto \|x_T\|` (i.e. :math:`\text{Dir}(x) = \frac{x}{\|x\|}` when :math:`x \ne 0` and :math:`\text{Dir}(0)` the closed unit ball). For two dimensional cases, :math:`\text{Dir}(\dot{u}_T)` reduces to :math:`\text{Sign}(\dot{u}_T)` where :math:`\text{Sign}` is the multivalued sign map.

A complete linearized elasticity problem with contact with friction reads as

Given an augmentation parameter :math:`r`, the contact and friction conditions can be equivalently expressed in term of projection as

.. math::

  \frac{1}{r}(\lambda_N^i - P_{]-\infty, 0]}(\lambda_N^i - r (u_N(a_i) - \text{gap}_i))) = 0,

  \frac{1}{r}(\lambda_T^i - P_{{\mathscr B}(-{\mathscr F}\lambda_N^i)}(\lambda_T^i - r \dot{u}_T(a_i))) = 0,

where :math:`P_K` is the projection on the convex :math:`K` and :math:`{\mathscr B}(-{\mathscr F}\lambda_N^i)` is the ball of center :math:`0` and radius :math:`-{\mathscr F}\lambda_N^i`.
These expressions will be used to perform a semi-smooth Newton method.

Suppose now that you approximate a linerized elasticity problem submitted to contact with friction. Then, if :math:`U` is the vector of the unknown for the displacement you will be able to express the matrices :math:`B_N` and :math:`B_T` such that

.. math::

  u_N(a_i) = (B_N U)_i,
 
  (\dot{u}_T(a_i))_k = (B_T \dot{U})_{(d-1)(i-1)+k},

where :math:`d` is the dimension of the domain and :math:`k = 1..d-1`. The expression of the elasticity problem with contact with friction can be written as

.. math::
 
  K U = L + B_N^T \lambda_N + B_T^T \lambda_T,

  \frac{1}{r}(\lambda_N^i - P_{]-\infty, 0]}(\lambda_N^i - \alpha_i r ((B_N U)_i - \text{gap}_i))) = 0, ~~ i = 1..N_c,

  \frac{1}{r}(\lambda_T^i - P_{{\mathscr B}(-{\mathscr F}\lambda_N^i)}(\lambda_T^i - \alpha_i r (B_T U - B_T U^{0})_i)) = 0, ~~ i = 1..N_c,

where :math:`\alpha_i` is a parameter which can be added for the homogenization of the augmentation parameter, :math:`(B_T U)_i` denotes here the sub-vector of indices from :math:`(d-1)(i-1)+1` to :math:`(d-1)i` for the sake of simplicity and the sliding velocity :math:`B_T \dot{U}` have been discretized into :math:`\frac{(B_T U - B_T U^{0})}{\Delta t}` whith :math:`U^{0}` the displacement at the previous time step. Note that of course another discretization of the sliding velocity is possible and that the time step :math:`\Delta t` do not appear in the expression of the friction condition since it does not influence the direction of the sliding velocity.


In that case, the homogeneization coefficient :math:`\alpha_i` can be taken

.. math::

  \alpha_i = \frac{\int_{\Gamma_c} \varphi_i d\Gamma}{\ell}

where :math:`\Gamma_c` is the contact boundary, :math:`\varphi_i` is the displacement shape function corresponding to the node :math:`a_i` and :math:`\ell` is a caracteristic lenght, for instance the radius of the domain. In this way, the augmentation parameter :math:`r` can be expressed in :math:`N/m^2` and chosen closed to the Young modulus of the elastic body. Note that the solution is not very sensitiv to the value of the augmentation parameter.


Weak nodal contact condition
++++++++++++++++++++++++++++

The direct nodal condition may have some drawback : locking phenomena, overconstraint. It is in fact often more stable and for the same accuracy to use multiplier of reduced order compared to the displacement (the direct nodal contact condition corresponds more or less to a multiplier described on the same finite element method than the displacement).

Let :math:`\varphi_i` be the shapes functions of the finite element describing the displacement and :math:`\psi_i` be the shape functions of a finite element describing a multiplier on the contact boundary :math:`\Gamma_c`. It is assumed that the set of admissible multiplier describing the normal stress will be

.. math::

  \Lambda_N^h = \{ \mu^h_N = \sum \mu^j_N \psi_j : \mu^h_N(a_i) \le 0, ~i = 1..N_c \}

where :math:`a_i`, :math:`~~i=1..N_c` are the finite element nodes corresponding to the multiplier. The discrete contact condition is now expressed in a weak form by

.. math::

  \int_{\Gamma_c} (\mu_N^h - \lambda_N^h) (u_N - \text{gap}) d\Gamma \ge 0 ~~ \forall \mu_N^h \in \Lambda_N^h. 

In that case, the component :math:`\lambda_N^i` is a contact stress (:math:`N/m^2`) and the matrix :math:`B_N` can be written

.. math::

  (B_N)_{ij} = \int_{\Gamma_c} \psi_i \varphi_j d\Gamma.

The matrix :math:`B_T` can also be written in a similar way. The friction condition can be written in a weak form

.. math::

  \int_{\Gamma_c} (\mu_T^h - \lambda_T^h) \dot{u}_T d\Gamma \ge 0 ~~ \forall \mu_T^h \in \Lambda_T^h({\mathscr F}\lambda_N^h),

where :math:`\Lambda_T^h({\mathscr F}\lambda_N^h)` is the disrete set of admissible friction stress.

Finally, the expression of the direct nodal contact condition are recovered 

.. math::
 
  K U = L + B_N^T \lambda_N + B_T^T \lambda_T,

  \frac{1}{r}(\lambda_N^i - P_{]-\infty, 0]}(\lambda_N^i - \alpha_i r ((B_N U)_i - \text{gap}_i))) = 0, ~~ i = 1..N_c,

  \frac{1}{r}(\lambda_T^i - P_{{\mathscr B}(-{\mathscr F}\lambda_N^i)}(\lambda_T^i - \alpha_i r (B_T U - B_T U^{0})_i)) = 0, ~~ i = 1..N_c,

except that now :math:`\lambda_N^i` and :math:`\lambda_T^i` are force densities, and a good value for :math:`\alpha_i` is now

.. math::

  \alpha_i = \frac{1}{\ell \int_{\Gamma_c}\psi_i},

where :math:`\psi_i` is the shape function of the multiplier for the node :math:`a_i`. In that case, the augmentation parameter :math:`r` can still be chosen close to the Young modulus of the elastic body.


Note that without additional stabilization technique (see [HI-RE2010]_) an inf-sup condition have to be satisfied between the finite element of the displacement and the one for the multipliers. This means in particular that the finite element for the multiplier have to be "less rich" than the one for the displacement.




Add a contact with or without friction to a model
+++++++++++++++++++++++++++++++++++++++++++++++++

Frictionless contact brick
++++++++++++++++++++++++++

In order to add a frictionless contact brick you call the model object method::

     getfem::add_basic_contact_brick
          (md, varname_u, multname_n, dataname_r, BN,dataname_gap , dataname_alpha, symmetrized );

This function adds a frictionless contact brick on ``varname_u`` thanks to a multiplier variable ``multname_n``. If we take :math:`U` is the vector of degrees of freedom on which the unilateral constraint is applied, the matrix :math:`B_N` have to be such that this condition is defined by :math:`B_N U \le 0`. The constraint is prescribed thank to a multiplier ``multname_n`` whose dimension should be equal to the number of lines of :math:`B_N`. The variable ``dataname_r`` is the name of the augmentation parameter :math:`r` should be chosen in a range of acceptabe values. ``dataname_gap`` is an optional parameter representing the initial gap. It can be a single value or a vector of value. ``dataname_alpha`` is an optional homogenization parameter for the augmentation parameter. The parameter ``symmetrized`` indicates that the symmetry of the tangent matrix will be kept or not. 
Note that is possible to change the basic contact matrix :math:`BN` by use::

     getfem::contact_brick_set_BN(md, indbrick);

Hughes stabilized frictionless contact condition
++++++++++++++++++++++++++++++++++++++++++++++++

In order to add a Hughes stabilized frictionless contact brick you call the model object method::

      getfem::add_Hughes_stab_basic_contact_brick
          (md, varname_u, multname_n, dataname_r, BN, DN, dataname_gap, dataname_alpha, symmetrized);

This function adds a Hughes stabilized frictionless contact brick on ``varname_u`` thanks to a multiplier variable ``multname_n``. If we take :math:`U` is the vector of degrees of freedom on which the unilateral constraint is applied, and :math:`\lambda` the multiplier Vector of contact force. Then Hughes stabilized frictionless contact condition is defined by the matrix :math:`BN` and :math:`DN` have to be such that this condition is defined by :math:`B_N U - D_N \lambda \le 0`. Where :math:`DN` is the masse matrix relative to stabilzed term. The variable ``dataname_r`` is the name of the augmentation parameter :math:`r` should be chosen in a range of acceptabe values. ``dataname_gap`` is an optional parameter representing the initial gap. It can be a single value or a vector of value. ``dataname_alpha`` is an optional homogenization parameter for the augmentation parameter. The parameter ``symmetrized`` indicates that the symmetry of the tangent matrix will be kept or not. 
Note that the matrix :math:`DN` is a sum of the basic contact term and the Hughes stabilised term. You can change it with::

      getfem::contact_brick_set_DN(md, indbrick);

