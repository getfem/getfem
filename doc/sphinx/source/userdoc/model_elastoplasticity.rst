.. $Id: model.rst 3655 2010-07-17 20:42:08Z renard $

.. include:: ../replaces.txt

.. highlightlang:: c++

.. index:: models, model bricks

.. _ud-model-elastoplasticity:



Elasto-plasticity brick
-----------------------

The aim of this brick is to add a nonlinear elasto-plastic term to a model.
This brick is restricted to small deformations on isotropic materials and for a quasistatic evolutive model.


Some recalls on elasticity
++++++++++++++++++++++++++


The phenomenon of elasticity refers to the fact that the material on which one applies constraints, it returns to its original form when the constraints are removed.

In order to model such a problem one has to consider:

- the second order small strain tensor :math:`\varepsilon`:

.. math::

   \varepsilon(u) = \frac{1}{2}(\nabla u + \nabla u^t)

where :math:`u` represents the displacements field of the solid.

- the second order symmetric stress tensor :math:`\sigma`

- the isotropic case which implies that:

.. math::

   \sigma = \lambda (Tr \ \varepsilon ) I + 2 \mu \varepsilon

where :math:`\lambda` and :math:`\mu` are the Lame coefficients.

- the fact that the evolutive model is considered as being quasistatic, which means that the model is subjected to the static equilibrium:

.. math::

   - div(\sigma) = f

where :math:`f` represents the volumic forces field applied on the solid by the external environment.

- the elastic constitutive law:

.. math::

   \sigma_{ij} = \sum_{kl} A_{ijkl} \varepsilon_{kl}

Finally, the problem to be solved is:

.. math::

   \text{Find } u \text{ and } \sigma \text{ in } \Omega \text{ such that :}
   \left\{
   \begin{array}{l}
	\sigma = \lambda (Tr \ \varepsilon ) I + 2 \mu \varepsilon \\
   	- div (\sigma) = f \\
   \text{+ boundary conditions} \\
   \end{array}  
   \right. 




Perfect elasto-plasticity problem
+++++++++++++++++++++++++++++++++

Contrary to the elastic phenomenon, the plasticity of a material is characterised by the onset of permanent deformations whithin the solid, resulted from the action of the constraints to which it is subjected.

Generally, these deformations appear beyond a certain stress threshold, noted :math:`s`. In fact, under this threshold the behavior of the material is linear and reversible, which means elastic. On the contrary, when :math:`\sigma` > s, permanent deformations appear and the behavior is not linear at all.

The permanent deformation is defined as the deformation :math:`\varepsilon_p` measured after the solid is unloaded so that the elastic deformation part :math:`\varepsilon_e` is retrieved.

The ideal case where the yield threshold is a material constant independent of the values reached by the plastic deformation, is called perfect elasto-plasticity.



Naturally, one can write the time derivative of the strain tensor as the sum of an elastic and a plastic part:

.. math::

   \dot{\varepsilon}(u) = \dot{\varepsilon_e}(u) + \dot{\varepsilon_p}(u)

Knowing that while :math:`\sigma` < s the material behaves elastically, one can write that:

.. math::

   \dot{\varepsilon_e} = C \dot{\sigma}

where :math:`C` is the elastic compliance tensor.


Then, one can consider :math:`K` = { :math:`\sigma \ / \ \varphi(\sigma) \leq 0` }, the convex defined by the set of admissible plastic constraints, where :math:`\varphi` corresponds to the objective function of plasticity which has to to be defined previously.
Thus, the plastic part :math:`\dot{\varepsilon_p}(u)` of the strain time derivative can be defined as the external normal of this convex on u:

.. math::

   \dot{\varepsilon} \in C \dot{\sigma} + \partial_{\sigma}I_K(\sigma)

where :math:`\partial_{\sigma}I_K(\sigma)` corresponds to the normal cone to :math:`K` on :math:`\sigma` .

This formulation of plasticity is known in the literature as the closest point projection method.

Finally, one has to solve the following problem:

.. math::

   \text{Find } u \text{ and } \sigma \text{ in } \Omega \text{ such that :}
   \left\{
   \begin{array}{l}
   \dot{\varepsilon}(u) \in C \dot{\sigma} + \partial_{\sigma}I_K(\sigma) \\
   - div (\sigma) = f \\
   \text{+ boundary  conditions} \\
   \end{array}
   \right.



Time discretisation
+++++++++++++++++++


One can perform a time discretisation with an implicit Euler scheme (unconditionally stable):

.. math::

   \frac{\varepsilon^{n+1} - \varepsilon^n}{\delta t} - C \frac{\sigma^{n+1} - \sigma^n}{\delta t} \in \partial_{\sigma} I_K(\sigma^{n+1})

.. math::

   \Leftrightarrow \varepsilon^{n+1} - \varepsilon^n - C \sigma^{n+1} + C \sigma^n \in \partial_{\sigma} I_K(\sigma^{n+1})



Weak formulation
++++++++++++++++


The weak problem associated to the above described elasto-plasticity problem is:

.. math::

   \left\{
   \begin{array}{l}
   \text{Find } u^{n+1} \in V = \left\{ v \in C^0_M(\Omega) \right\} \text{ such that :} \\
    \\
   \int_{\Omega} \sigma^{n+1}(\varepsilon(u^{n+1})) : \varepsilon(v) dx = \int_{\Omega} f v \ dx \ \ \ \forall v \in V \\
   \end{array}
   \right.


################

**Property:**

.. math::
   
   \begin{array}{ll}
   \alpha \in \partial_{\sigma}I_K(\sigma) & \Leftrightarrow \sigma = P_K(\sigma + \alpha) \ \ \ \ \forall \alpha \\
   & \Leftrightarrow (\tau - \sigma):(\alpha) \leq 0 \ \ \ \forall \tau \in K \\
   \end{array}

where :math:`P_K` represents the projection operator on :math:`K` associated with the usual scalar product. 

###############


Thus, according to this property, one has:

.. math::

   \begin{array}{l}
   \varepsilon^{n+1} - \varepsilon^n - C \sigma^{n+1} + C \sigma^n \in \partial_{\sigma} I_K(\sigma^{n+1}) \\
   \Leftrightarrow C ( \underbrace{A \varepsilon^{n+1} - A \varepsilon^n + \sigma^n}_\beta - \sigma^{n+1}  ):( \tau - \sigma^{n+1} ) \leq 0 \ \ \ \forall \tau \in K \\
   \Leftrightarrow C ( \beta - \sigma^{n+1} ):( \tau - \sigma^{n+1} ) \leq 0 \ \ \ \forall \tau \in K \\
   \end{array}

with :math:`C = A^{-1}` where :math:`A` is the fourth order, symetric and real, elastic stiffness tensor. :math:`A` is diagonalizable and invertible.

Thus, one has the orthogonality in the sense of the scalar product of elasticity associated with the fourth order tensor :math:`C` and so one can write:

.. math::
   
   \sigma^{n+1} = P_K^C(A \varepsilon^{n+1} - A \varepsilon^n + \sigma^n)

where :math:`P_K^C` represents the projection operator on :math:`K` associated with the elastic scalar product defined above.

Finally, the following weak problem has to be solved:

.. math::

   \left\{
   \begin{array}{l}
   \text{Find } u^{n+1} \in V \text{ such that :} \\
    \\
   \int_\Omega P_K^C(A \varepsilon^{n+1} - A \varepsilon^n + \sigma^n): \varepsilon(v) \ dx = \int_\Omega f v \ dx \ \ \ \forall v \in V \\
   \end{array}
   \right.

In the last equations, the notation :math:`\varepsilon^n = \varepsilon(u^n)` was used for the sake of simplicity.




Elastic projection operator derivative
++++++++++++++++++++++++++++++++++++++


In order to apply a Newton algorithm and thus to obtain the solution of the problem in terms of the evolution time parameter :math:`n`, one has to determine the derivative of :math:`P_K^C(\tau) \ \forall \tau` with respect to :math:`u^{n+1}` , which will be denoted as :math:`\nabla P_K^C` .

By definition, all tensors :math:`\tau` could be decomposed as the sum of a spherical and a deviatoric part as follows:

.. math::

   \tau = \tau^S + \tau^D := \tau_m I + \tau^D, \ \tau_m = \frac{1}{N} Tr(\tau)

where :math:`N` is the dimension of the considered problem.


Moreover, one could prove that: :math:`P_K^C \equiv P_K` , admitted here.

Thus, finding :math:`\nabla P_K^C` is equivalent to finding :math:`\nabla P_K` and it is known that:

.. math::
   
   P_K(\tau) = \tau_m I + inf(|\tau^D|, s) \frac{\tau^D}{|\tau^D|}, \ \forall \tau, \ |\tau^D| > 0

where :math:`|\tau| = (\tau : \tau)^{1/2}` .

In following, three cases have to be considered.


Classical linear elasticity: :math:`|\tau^D| < s` 
##################################################


Here, :math:`\sigma` is whithin the convex :math:`K` and so the projection can be written as following:

.. math::

   P_K(\tau) = \tau_m I + |\tau^D| \frac{\tau^D}{|\tau^D|} = \tau

In that case, :math:`P_K(\tau)` is differentiable with:

.. math::

   <\nabla P_K(\tau), \tau^*> = \tau^*

Thus, :math:`\nabla P_K(\tau) = I_S` , where :math:`I_S` represents here the fourth order identity tensor.


   

Plastic scheme: :math:`|\tau^D| > s` 
#####################################


Here, :math:`\sigma` is out of the convex :math:`K` and the projection can be written as following:

.. math::

   P_K(\tau) = \tau_m I + s \frac{\tau^D}{|\tau^D|}

which is differentiable.

Moreover, knowing that:

.. math::

   h(x) = |x| \ \Rightarrow \ <h'(x), y> = \frac{(x.y)}{|x|}

and that:

.. math::

   g(x) = R \frac{x}{|x|}, \ R \in \Re \ \Rightarrow \ <g'(x), y> = \frac{R}{|x|} [y - (n^*.y) n^*]

with :math:`n^* = \frac{x}{|x|}` and the operator `.` representing the usual scalar product. 

Thus, knowing that the application :math:`\tau \rightarrow \tau^D` is linear, one has:

.. math::

   <\nabla P_K(\tau), \tau^*> = \tau_m^* I + \frac{s}{|\tau|}[{\tau^D}^* - (n:{\tau^D}^*)n]

with :math:`n = \frac{\tau^D}{|\tau^D|}` .

Then, introducing the operator:

.. math::

   I^D : \tau \rightarrow \tau^D

and the relations:

.. math::

   \begin{array}{c}
   (u \otimes v)w = (v.w)u \\
   Tr(\tau) I = (I \otimes I)\tau \\   
   \end{array}

the derivative of the projection becomes:

.. math::

   <\nabla P_K(\tau), \tau^*> = \frac{1}{N}(I \otimes I)\tau^* + \frac{s}{|\tau^D|}[I_S - n \otimes n]I^D {\tau^D}^*

Thus, :math:`\nabla P_K(\tau) = \frac{1}{N}(I \otimes I) + \frac{s}{|\tau^D|}[I_S - n \otimes n]I^D` .




Elastic threshold case: :math:`|\tau^D| = s`
#############################################


Here, one has:

.. math::

   P_K(\tau) = \tau_m I + \tau^D = \tau

In that case, :math:`P_K(\tau)` is not differentiable and its derivative depends on the direction considered:

.. math::

   <\nabla P_K(\tau), \tau^* > = 
   \left\{
   \begin{array}{l l}
   I_S \tau^* & \text{if } \tau^* \in \text{ tangent cone to } K \\
   (I_S - n \otimes n)\tau^* & \text{otherwise} \\
   \end{array}
   \right.




Assembly of Newton's terms
++++++++++++++++++++++++++


In order to apply a Newton algorithm, a tangent matrix and a right hand side vector have to be calculated.


In this problem, the tangent matrix corresponds to the term:

.. math::

   T \equiv \int_\Omega \nabla P_K(A \varepsilon^{n+1} - A \varepsilon^n + \sigma^n) : A \varepsilon(u^*) : \varepsilon(v) \ dx


and the right hand side vector corresponds to:

.. math::

   R \equiv \int_\Omega P_K(A \varepsilon^{n+1} - A \varepsilon^n + \sigma^n) : \varepsilon(v) \ dx - \int_\Omega fv \ dx

Of course, one should add some boundary conditions using appropriate bricks.




Discrete assembly of the terms
++++++++++++++++++++++++++++++


If one denotes:

.. math::

   u \in V_h \ : \ u = \sum_{i = 1}^{N_u} u_i \varphi_i

where :math:`u_i \in \Re` and :math:`\varphi_i : \Omega \rightarrow \Re^N` with :math:`N` the dimension of the problem,

and:

.. math::

   \sigma \in W_h \ : \ \sigma = \sum_{i = 1}^{N_\sigma} \sigma_i \psi_i

where :math:`\sigma_i \in M_{3,3}(\Re)` and :math:`\psi_i : \Omega \rightarrow \Re`,

one has:

.. math::

   R_i = \sum_T \int_T P_K(A \varepsilon^{n+1} - A \varepsilon^n + \sigma^n) : \varepsilon(\varphi_i) \ dx - \sum_T \int_T f_i \ \varphi_i \ dx  


.. math::

   R_i \simeq \sum_T \int_T \sum_{k = 1}^{N_\sigma} [P_K(A \varepsilon^{n+1}(a_{i_k}) - A \varepsilon^n(a_{i_k}) + \sigma^n(a_{i_k})) \psi_{i_k}] : \varepsilon(\varphi_i) \ dx - \sum_T \int_T f_i \ \varphi_i \ dx

where :math:`a_{i_k}` are the nodes of :math:`W_h` on the element T,

and also:

.. math::

   \frac{\partial R}{\partial u_i}[h] \simeq \sum_T \int_T \sum_{k = 1}^{N_\sigma} [\nabla P_K(A \varepsilon^{n+1}(a_{i_k}) - A \varepsilon^n(a_{i_k}) + \sigma^n(a_{i_k})) \psi_{i_k}] : A \varepsilon(h) : \varepsilon(\varphi_i) \ dx

where :math:`h \in V_h` .


In order to compute such a projection, one chooses to interpolate :math:`\varepsilon^{n}` and :math:`\varepsilon^{n+1}` directly on :math:`\sigma` dofs to make sure that the sum will be correctly computed, then to compute the projection on each dofs of :math:`\sigma` and finally to interpolate the result on the dofs of :math:`u` for the integration and the assembly.



Add an elasto-plasticity brick to a model
+++++++++++++++++++++++++++++++++++++++++

The function adding this brick to a model is: ::

      getfem::add_elastoplasticity_brick
          (md, mim, ACP, varname, datalambda, datamu, datathreshold, datasigma, region);

where:
      - ``varname`` represents the main unknown on which the brick is added (u). It should be composed of 2 iterates for the time scheme needed for the Newton algorithm used.
      - ``datalambda`` and ``datamu`` are the data corresponding to the Lame coefficients.
      - ``datathreshold`` represents the plastic threshold of the studied material.
      - ``datasigma`` represents the stress constraint values supported by the material. It should be composed of 2 iterates for the time scheme needed for the Newton algorithm used. Note that the finite element method on which ``datasigma`` is defined should be able to represent the derivative of ``varname``.
      - ``ACP`` corresponds to the type of projection to be used. It has an `abstract_constraints_projection` type and for the moment, only exists the `VM_projection` corresponding to the Von Mises one.


Be careful: ``datalambda``, ``datamu`` and ``datathreshold`` could be constants or described on the same finite element method.

This function assembles the tangent matrix and the right hand side vector which will be solved using a Newton algorithm.


Other useful functions
++++++++++++++++++++++

The function: ::

      getfem::elastoplasticity_next_iter
          (md, mim, varname, ACP, datalambda, datamu, datathreshold, datasigma);

computes the new stress constraint values supported by the material after a load or an unload (once a solve has been done earlier) and upload the variables ``varname`` and ``datasigma`` as follows:

.. math::
   
   u^{n+1} \Rightarrow u^n \ \ \ \ \ and \ \ \ \ \ \sigma^{n+1} \Rightarrow \sigma^n

Then, :math:`u^n` and :math:`\sigma^n` contains the new values computed and one can restart the process.



########################


The function: ::

      getfem::compute_elastoplasticity_Von_Mises_or_Tresca
          (md, datasigma, mf_vm, VM, tresca=false);

computes the Von Mises (or Tresca if ``tresca`` = true) criterion on the stress tensor stored in ``datasigma`` . The stress is evaluated on the `mesh_fem` ``mf_vm`` and stored into the vector ``VM``.
Of course, this function can be used if and only if the previous function ``elastoplasticity_next_iter`` has been called earlier.



##########################


The function: ::

      getfem::compute_plastic_part
          (md, mim, mf_pl, varname, ACP, datalambda, datamu, datathreshold, datasigma, Plast);

computes on ``mf_pl`` the plastic part of the material, that could appear after a load and an unload, into the vector ``Plast``. 

Note that ``datasigma`` should be the vector containing the new stress constraint values, i.e. after a load or an unload of the material.





The program ``tests/plasticity.cc`` can be taken as a model of use of this brick.


    
