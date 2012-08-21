.. $Id: model.rst 3655 2010-07-17 20:42:08Z renard $

.. include:: ../replaces.txt

.. highlightlang:: c++

.. index:: models, model bricks

.. _ud-model-elastoplasticity:



Elasto-plasticity brick
-----------------------

The aim of this brick is to add a nonlinear elasto-plastic term to a model.
This brick is restricted to small deformation on isotropic materials and for a quasistatic evolutive model.


Some recalls on elasticity
++++++++++++++++++++++++++


The phenomenon of elasticity leads on the fact that the material on which one applie constraints returns to its original form when the constraints disappears.

In order to modelize such a problem one have to consider :

- the two ordered small strain tensor :math:`\varepsilon` :

.. math::

   \varepsilon(u) = \frac{1}{2}(\nabla u + \nabla u^t)

where u represents the displacement of the solid.

- the two ordered symetric stress tensor :math:`\sigma`

- the isotropic case which implies that :

.. math::

   \sigma = \lambda (Tr \ \varepsilon ) I + 2 \mu \varepsilon

where :math:`\lambda` and :math:`\mu` are the Lame coefficients.

- the fact that the evolutive model is considered as being quasistatic, which means that the model is subjected to the static fondamentale law :

.. math::

   - div(\sigma) = f

where f represents the volumic force applied on the solid by the external environment.

- the elastic constitutive law :

.. math::

   \sigma_{ij} = \sum_{kl} A_{ijkl} \varepsilon_{kl}

Finally the problem to be solved is :

.. math::

   Find \ u \ and \ \sigma \ in \ \Omega \ such \ as  \ :
   \left\{
   \begin{array}{l}
	\sigma = \lambda (Tr \ \varepsilon ) I + 2 \mu \varepsilon \\
   	- div (\sigma) = f \\
   	+ \ boundary \ conditions \\
   \end{array}  
   \right. 




Perfect elasto-plasticity problem
+++++++++++++++++++++++++++++++++

Contrary to the elastic phenomenon, the plasticity of a material is caracterised by the apparition of permanent deformations whithin the solid, resulted from the action of the constraints to which it is subjected.

Generally, these deformations appear beyond a certain stress threshold, noted ``s``. In fact, under this threshold the behavior of the material is linear and revertible, which means elastic, contrariwise, when :math:`\sigma` > s, permanent deformations appear and the behavior is not linear at all.

The permanent deformation is defined as the deformation :math:`\varepsilon_p` measured after an elastic unload and which retrieve the elastic deformation :math:`\varepsilon_e`.

One consider the ideal case where the flow threshold is a material constant independent of the values reached by the plastic deformation, it is the perfect elasto-plasticity.



Naturally, one write the time derivative of the strain tensor as a sum between an elastic part and a plastic one :

.. math::

   \dot{\varepsilon}(u) = \dot{\varepsilon_e}(u) + \dot{\varepsilon_p}(u)

Knowing that while :math:`\sigma` < s the material behaves elastically, one can write that :

.. math::

   \dot{\varepsilon_e} = C \dot{\sigma}

where C is the elastic compliance tensor.


Then, one note K = { :math:`\sigma \ / \ \varphi(\sigma) \leq 0` }, the convex defined by the set of admissible plastic constraints, where :math:`\varphi` corresponds to the objective function of plasticity which have to to be defined previously.
Thus, one can write the plastic part :math:`\dot{\varepsilon_p}(u)` as being the external normal of this convex on u :

.. math::

   \dot{\varepsilon} \in C \dot{\sigma} + \partial_{\sigma}I_K(\sigma)

where :math:`\partial_{\sigma}I_K(\sigma)` corresponds to the normal cone to K on :math:`\sigma` .

Finally one have to solve the following problem :

.. math::

   Find \ u \ and \ \sigma \ in \ \Omega \ such \ as  \ :
   \left\{
   \begin{array}{l}
   \dot{\varepsilon}(u) \in C \dot{\sigma} + \partial_{\sigma}I_K(\sigma) \\
   - div (\sigma) = f \\
   + \ boundary  \ condition \\
   \end{array}
   \right.



Time discretisation
+++++++++++++++++++


One perform a time discretisation with an implicit Euler scheme (unconditionally stable) :

.. math::

   \frac{\varepsilon^{n+1} - \varepsilon^n}{\delta t} - C \frac{\sigma^{n+1} - \sigma^n}{\delta t} \in \partial_{\sigma} I_K(\sigma^{n+1})

.. math::

   \Leftrightarrow \varepsilon^{n+1} - \varepsilon^n - C \sigma^{n+1} + C \sigma^n \in \partial_{\sigma} I_K(\sigma^{n+1})



Weak formulation
++++++++++++++++


The weak problem associated to our elasto-plasticity problem is :

.. math::

   \left\{
   \begin{array}{l}
   Find \ u^{n+1} \in V = \left\{ v \in C^0_M(\Omega) \right\} \ such \ as \ : \\
    \\
   \int_{\Omega} \sigma^{n+1}(\varepsilon(u^{n+1})) : \varepsilon(v) dx = \int_{\Omega} f v \ dx \ \ \ \forall v \in V \\
   \end{array}
   \right.


################

**Property :**

.. math::
   
   \begin{array}{ll}
   \alpha \in \partial_{\sigma}I_K(\sigma) & \Leftrightarrow \sigma = P_K(\sigma + \alpha) \ \ \ \ \forall \alpha \\
   & \Leftrightarrow (\tau - \sigma):(\alpha) \leq 0 \ \ \ \forall \tau \in K \\
   \end{array}

where :math:`P_K` represents the projection operator on K associated with the usual scalar product. 

###############


Thus, according to this property, one have :

.. math::

   \begin{array}{l}
   \varepsilon^{n+1} - \varepsilon^n - C \sigma^{n+1} + C \sigma^n \in \partial_{\sigma} I_K(\sigma^{n+1}) \\
   \Leftrightarrow C ( \underbrace{A \varepsilon^{n+1} - A \varepsilon^n + \sigma^n}_\beta - \sigma^{n+1}  ):( \tau - \sigma^{n+1} ) \leq 0 \ \ \ \forall \tau \in K \\
   \Leftrightarrow C ( \beta - \sigma^{n+1} ):( \tau - \sigma^{n+1} ) \leq 0 \ \ \ \forall \tau \in K \\
   \end{array}

with :math:`C = A^{-1}` where A is the four ordered, symetric and real, elastic stiffness tensor. A is diagonalizable and invertible.

Thus, one have the orthogonality in the sense of the scalar product of elasticity associated with the four ordered tensor C and so one can write :

.. math::
   
   \sigma^{n+1} = P_K^C(A \varepsilon^{n+1} - A \varepsilon^n + \sigma^n)

where :math:`P_K^C` represents the projection operator on K associated with the elastic scalar product defined above.

Finally, one have to solve the following weak problem :

.. math::

   \left\{
   \begin{array}{l}
   Find \ u^{n+1} \in V \ such \ as \ : \\
    \\
   \int_\Omega P_K^C(A \varepsilon^{n+1} - A \varepsilon^n + \sigma^n): \varepsilon(v) \ dx = \int_\Omega f v \ dx \ \ \ \forall v \in V \\
   \end{array}
   \right.

where one noted :math:`\varepsilon^n = \varepsilon(u^n)` .




Elastic projection operator derivative
++++++++++++++++++++++++++++++++++++++


In order to apply farther a Newton algorithm and thus to obtain the solution of the problem in terms of the evolution time parameter n, one have to determine the derivative, in regards to :math:`u^{n+1}` , of :math:`P_K^C(\tau) \ \forall \tau` ,that one will note :math:`\nabla P_K^C` .

By definition, all tensors :math:`\tau` could be decomposed as a sum between its spherical part and its deviatoric part as follows :

.. math::

   \tau = \tau^S + \tau^D := \tau_m I + \tau^D, \ \tau_m = \frac{1}{N} Tr(\tau)

where N is the dimension of the considered problem.


Moreover, one could prove that : :math:`P_K^C \equiv P_K` , admitted here.

Thus, to find :math:`\nabla P_K^C` is equivalent than to find :math:`\nabla P_K` , and one know that :

.. math::
   
   P_K(\tau) = \tau_m I + inf(|\tau^D|, s) \frac{\tau^D}{|\tau^D|}, \ \forall \tau, \ |\tau^D| > 0

where :math:`|\tau| = (\tau : \tau)^{1/2}` .

Three cases have to be considered.


Classical linear elasticity : :math:`|\tau^D| < s` 
##################################################


Here, :math:`\sigma` is whithin the convex K and so the projection can be written as follows :

.. math::

   P_K(\tau) = \tau_m I + |\tau^D| \frac{\tau^D}{|\tau^D|} = \tau

In that case, :math:`P_K(\tau)` is differentiable and one have :

.. math::

   <\nabla P_K(\tau), \tau^*> = \tau^*

Thus, :math:`\nabla P_K(\tau) = I_S` , where :math:`I_S` represents here the four ordered identity tensor.


   

Plastic scheme : :math:`|\tau^D| > s` 
#####################################


Here, :math:`\sigma` is out of the convex K and the projection can be write as follows :

.. math::

   P_K(\tau) = \tau_m I + s \frac{\tau^D}{|\tau^D|}

which is differentiable.

Moreover, knowing that :

.. math::

   h(x) = |x| \ \Rightarrow \ <h'(x), y> = \frac{(x.y)}{|x|}

and that :

.. math::

   g(x) = R \frac{x}{|x|}, \ R \in \Re \ \Rightarrow \ <g'(x), y> = \frac{R}{|x|} [y - (n^*.y) n^*]

with :math:`n^* = \frac{x}{|x|}` and the operator `.` representing the used scalar product. 

Thus knowing that the application :math:`\tau \rightarrow \tau^D` is linear, one have that :

.. math::

   <\nabla P_K(\tau), \tau^*> = \tau_m^* I + \frac{s}{|\tau|}[{\tau^D}^* - (n:{\tau^D}^*)n]

with :math:`n = \frac{\tau^D}{|\tau^D|}` .

Then, introducing the operator :

.. math::

   I^D : \tau \rightarrow \tau^D

and the relations :

.. math::

   \begin{array}{c}
   (u \oplus v)w = (v.w)u \\
   Tr(\tau) I = (I_S \otimes I_S)\tau \\   
   \end{array}

the derivative of the projection becomes :

.. math::

   <\nabla P_K(\tau), \tau^*> = \frac{1}{N}(I_S \otimes I_S)\tau^* + \frac{s}{|\tau^D|}[I_S - n \otimes n]I^D {\tau^D}^*

Thus, :math:`\nabla P_K(\tau) = \frac{1}{N}(I_S \otimes I_S) + \frac{s}{|\tau^D|}[I_S - n \otimes n]I^D` .




Elastic threshold case : :math:`|\tau^D| = s`
#############################################


Here, one have :

.. math::

   P_K(\tau) = \tau_m I + \tau^D = \tau

In that case, :math:`P_K(\tau)` is not differentiable and its derivative depends on the direction considered :

.. math::

   <\nabla P_K(\tau), \tau^* > = 
   \left\{
   \begin{array}{l l}
   I_S \tau^* & if \ \tau^* \in \ tangent \ cone \ to \ K \\
   (I_S - n \otimes n)\tau^* & otherwise \\
   \end{array}
   \right.




Assembly of Newton's terms
++++++++++++++++++++++++++



In order to apply a Newton algorithm, one have to give a tangent matrix and a right hand side vector.


In this problem, the tangent matrix corresponds to the term :

.. math::

   T \equiv \int_\Omega \nabla P_K(A \varepsilon^{n+1} - A \varepsilon^n + \sigma^n) : A \varepsilon(u^*) : \varepsilon(v) \ dx


and the right hand side vector corresponds to : 

.. math::

   R \equiv \int_\Omega P_K(A \varepsilon^{n+1} - A \varepsilon^n + \sigma^n) : \varepsilon(v) \ dx - \int_\Omega fv \ dx

Of course one should add some boundary condition thanks to appropriated bricks.




Discrete assembly of the terms
++++++++++++++++++++++++++++++


If one denotes:

.. math::

   u \in V_h \ : \ u = \sum_{i = 1}^{N_u} u_i \varphi_i

where :math:`u_i \in \Re` and :math:`\varphi_i : \Omega \rightarrow \Re^N` with N the dimension of the problem;

.. math::

   \sigma \in W_h \ : \ \sigma = \sum_{i = 1}^{N_\sigma} \sigma_i \psi_i

where :math:`\sigma_i \in M_{3,3}(\Re)` and :math:`\psi_i : \Omega \rightarrow \Re` ;

one has:

.. math::

   R_i = \sum_T \int_T P_K(A \varepsilon^{n+1} - A \varepsilon^n + \sigma^n) : \varepsilon(\varphi_i) \ dx - \sum_T \int_T f_i \ \varphi_i \ dx  


.. math::

   R_i \simeq \sum_T \int_T \sum_{k = 1}^{N_\sigma} [P_K(A \varepsilon^{n+1}(a_{i_k}) - A \varepsilon^n(a_{i_k}) + \sigma^n(a_{i_k})) \psi_{i_k}] : \varepsilon(\varphi_i) \ dx - \sum_T \int_T f_i \ \varphi_i \ dx

where :math:`a_{i_k}` are the nodes of :math:`W_h` on the element T;

and also:

.. math::

   \frac{\partial R}{\partial u_i}[h] \simeq \sum_T \int_T \sum_{k = 1}^{N_\sigma} [\nabla P_K(A \varepsilon^{n+1}(a_{i_k}) - A \varepsilon^n(a_{i_k}) + \sigma^n(a_{i_k})) \psi_{i_k}] : A \varepsilon(h) : \varepsilon(\varphi_i) \ dx

where :math:`h \in V_h` .


In order to compute such a projection, one chooses to interpolate :math:`\varepsilon^{n}` and :math:`\varepsilon^{n+1}` directly on :math:`\sigma` dofs to make sure that the sum will be correctly computed, then to compute the projection on each dofs of :math:`\sigma` and finally to interpolate the result on the dofs of u for the integration and the assembly.



Add an elastoplaticity brick to a model
+++++++++++++++++++++++++++++++++++++++

The function adding this brick to a model is ::

      getfem::add_elastoplasticity_brick
          (md, mim, ACP, varname, datalambda, datamu, datathreshold, datasigma, region);

where :
      - ``varname`` represents the main unknown on which the brick is added (u). It should be composed of 2 iterates for the time scheme needed for the Newton algorithm used.
      - ``datalambda`` and ``datamu`` are the data corresponding to the Lame coefficients.
      - ``datathreshold`` represents the plastic threshold of the studied material.
      - ``datasigma`` represents the stress constraint values supported by the material. It should be composed of 2 iterates for the time scheme needed for the Newton algorithm used. Note that the finite element method on which ``datasigma`` is defined should be able to represent the derivative of ``varname``.
      - ``ACP`` corresponds to the type of projection to be used. It has an `abstract_constraints_projection` type and for the moment, only exists the `VM_projection` corresponding to the Von Mises one.


Be careful : ``datalambda``, ``datamu`` and ``datathreshold`` could be constants or described on the same finite element method.

This function assembles the tangent matrix and the right hand side vector which will be solve using a Newton algorithm.


Other useful functions
++++++++++++++++++++++

The function : ::

      getfem::elastoplasticity_next_iter
          (md, mim, varname, ACP, datalambda, datamu, datathreshold, datasigma);

compute the new stress constraint values supported by the material after a load or an unload (once a solve has been done earlier) and upload the variables ``varname`` and ``datasigma`` as follows :

.. math::
   
   u^{n+1} \Rightarrow u^n \ \ \ \ \ and \ \ \ \ \ \sigma^{n+1} \Rightarrow \sigma^n

Then, :math:`u^n` and :math:`\sigma^n` contains the new values computed and one can restart the process.



########################


The function : ::

      getfem::compute_elastoplasticity_Von_Mises_or_Tresca
          (md, datasigma, mf_vm, VM, tresca=false);

compute the Von Mises (or Tresca if ``tresca`` = true) criterion on the stress tensor stored in ``datasigma`` . The stress is evaluated on the `mesh_fem` ``mf_vm`` and stored into the vector ``VM``.
Of course, this function can be used if and only if the previous function ``elastoplasticity_next_iter`` has been done earlier.



##########################


The function : ::

      getfem::compute_plastic_part
          (md, mim, mf_pl, varname, ACP, datalambda, datamu, datathreshold, datasigma, Plast);

compute on ``mf_pl`` the plastic part of the material, that could appears after a load and an unload, into the vector ``Plast``. 

Note that ``datasigma`` should be the vector containing the new stress constraint values, ie after a load or an unload of the material.





The program ``tests/plasticity.cc`` can be taken as a model of use of this brick.


    
