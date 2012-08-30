.. $Id: model.rst 3655 2010-07-17 20:42:08Z renard $

.. include:: ../replaces.txt

.. highlightlang:: c++

.. index:: models, model bricks

.. _ud-model-continuation:



Numerical continuation and bifurcation
--------------------------------------

Let an algebraic problem coming from a discretization of a FEM-model can be
written in the form

.. math::

   F(U) = 0.

In what follows, we shall suppose that the model depends on an additional scalar
parameter :math:`\lambda` so that :math:`F(U) = F(U, \lambda)`. 

Numerical continuation
++++++++++++++++++++++

A numerical continuation method traces solution branches of the system

.. math::

   F(U, \lambda) = 0, \quad F:\mathbb{R}^{N} \times \mathbb{R} \to \mathbb{R}^{N}.

In |gf|, the (approximate) *Moore-Penrose* (also called *Gauss-Newton*)
continuation is implemented (see, for instance, [dh-go-ku2003]_).

Since this method does not make an explicit difference between the state
variable :math:`U` and the parameter :math:`\lambda`, we shall denote
:math:`Y := (U, \lambda)` for brevity. Nevertheless, to avoid bad scaling of the
values of the continuation parameter, we shall use the following weighted scalar
product and norm:

.. math::

   \langle Y, \tilde{Y} \rangle_{w} := \kappa \langle U, \tilde{U} \rangle + \lambda \tilde{\lambda},\quad \lVert Y \rVert_{w} := \sqrt{\kappa \lVert U \rVert^{2} + \lambda^{2}},\qquad Y = (U, \lambda),\, \tilde{Y} = (\tilde{U}, \tilde{\lambda}).

Here, :math:`\kappa` should be chosen so that
:math:`\kappa \langle U, \tilde{U} \rangle` is proportional to the scalar
product of the corresponding space variables in :math:`L^{2}`. One can take, for
example, :math:`\kappa = h^{d}`, where :math:`h` is the mesh size and :math:`d`
stands for the dimension of the problem. Alternatively, :math:`\kappa` can be
chosen as the reciprocal of the total number of degrees of freedom.

The Moore-Penrose continuation consists in computing a sequence of consecutive
points :math:`Y_{j}` on a chosen solution branch and the corresponding unit
tangent vectors :math:`T_{j}`:

.. math::

   F(Y_{j}) = 0,\quad \nabla F(Y_{j}) T_{j} = 0,\quad \lVert T_{j} \rVert_{w} = 1,\quad j = 0, 1,\dotsc.

To describe the technique, let us suppose that we have a couple
:math:`(Y_{j}, T_{j})` satisfying the relations above at our disposal. The next
couple is calculated in two steps -- *prediction* and *correction*.

In the prediction, an initial approximation of :math:`(Y_{j+1}, T_{j+1})` is
given by

.. math::

   Y_{j+1}^{0} := Y_{j} + h_{j} T_{j},\quad T_{j+1}^{0} := T_{j},

where :math:`h_{j}` is a step size. Its choice will be discussed later on.

In the correction, one computes a sequence
:math:`\{(Y_{j+1}^{l}, T_{j+1}^{l})\}`, where
:math:`T_{j+1}^{l} := \tilde{T}_{j+1}^{l} / \lVert \tilde{T}_{j+1}^{l} \rVert_{w}`
and the couple :math:`(Y_{j+1}^{l}, \tilde{T}_{j+1}^{l})` is given by one
iteration of the Newton method applied to the equation :math:`F^{l}(Y, T) = 0`
with

.. math::

   F^{l}(Y, T) := \begin{pmatrix}F(Y)\\ (T_{j+1}^{l-1})^{\top}(Y - Y_{j+1}^{l-1})\\ \nabla F(Y_{j+1}^{l-1})T\\ \langle T_{j+1}^{l-1}, T \rangle_{w} - \langle T_{j+1}^{l-1}, T_{j+1}^{l-1} \rangle_{w}\end{pmatrix}.

.. _ud_fig_correction:
.. figure:: images/getfemusercorrection.png
   :align: center

   Correction

A new couple :math:`(Y_{j+1}, T_{j+1})` is set to
:math:`(Y_{j+1}^{l}, T_{j+1}^{l})` iff
:math:`\lVert F(Y_{j+1}^{l})\rVert \leq \varepsilon`,
:math:`\lVert Y_{j+1}^{l} - Y_{j+1}^{l-1}\rVert_{w} \leq \varepsilon'`, and the
cosine of the angle between :math:`T_{j+1}^{l}` and :math:`T_{j}` is greater or
equal to :math:`c_{\mathrm{min}}`. Let us note that the partial gradient
:math:`\nabla_{U} F` is assembled analytically whereas
:math:`\nabla_{\lambda} F` is evaluated by forward finite differences with an
increment equal to 1e-8. 

Finally, the step size :math:`h_{j+1}` in the next prediction depends on how
this Newton correction is successful. Denoting the number of iterations needed
by :math:`l_{\mathrm{it}}`, it is selected as

.. math::

   h_{j+1} := \begin{cases}\max\{h_{\mathrm{dec}} h_{j}, h_{\mathrm{min}}\}& \text{if no new couple was accepted},\\ \min\{h_{\mathrm{inc}} h_{j}, h_{\mathrm{max}}\}& \text{if a new couple was accepted and } l_{\mathrm{it}} < l_{\mathrm{thr}},\\ h_{j}& \text{otherwise},\end{cases}

where :math:`0 < h_{\mathrm{dec}} < 1 < h_{\mathrm{inc}}`,
:math:`0 < l_{\mathrm{thr}}` as well as
:math:`0 < h_{\mathrm{min}} < h_{\mathrm{max}}` are given constants. At the
beginning, one sets :math:`h_{1} := h_{\mathrm{init}}` for some
:math:`h_{\mathrm{min}} \leq h_{\mathrm{init}} \leq h_{\mathrm{max}}`.

In |gf|, the Moore-Penrose continuation is implemented for two ways of
parametrisation of the model:

1. The parameter :math:`\lambda` is directly a scalar datum that the model 
   depends on.

2. The model is parametrised by the scalar parameter :math:`\lambda` *via* a
   vector datum :math:`P` that the model depends on. In this case, one takes the
   linear path

   .. math::

      \lambda \mapsto P(\lambda) := (1 - \lambda)P^{0} + \lambda P^{1},

   where :math:`P^{0}` and :math:`P^{1}` are given values of :math:`P`, and one
   traces the solution set of the problem

   .. math::

      F(U, P(\lambda)) = 0.

Numerical bifurcation
+++++++++++++++++++++

A point :math:`\bar{Y}` is called a *bifurcation point* of the equation
:math:`F(Y) = 0` if :math:`F(\bar{Y}) = 0` and two or more distinct solution
branches pass through it. The following result gives a test for bifurcation
points (see, e.g., [georg2001]_):

Let :math:`s \mapsto Y(s)` be a parametrisation of a solution branch and
:math:`\bar{Y} := Y(\bar{s})` a bifurcation point. Moreover, let
:math:`T^{\top} \dot{Y}(\bar{s}) > 0` and
:math:`B \notin \mathrm{im}(J(\bar{Y}))`,
:math:`C \notin \mathrm{im}(J(\bar{Y})^{\top})` with

   .. math::

      J(Y) := \begin{pmatrix}\nabla F(Y)\\ T^{\top}\end{pmatrix}.

Define :math:`\tau(Y)` via

   .. math::

      \begin{pmatrix}J(Y)& B\\ C^{\top}& 0\end{pmatrix} \begin{pmatrix}V(Y)\\ \tau(Y)\end{pmatrix} = \begin{pmatrix}0\\ 1\end{pmatrix}.

Then :math:`\tau(Y(s))` changes sign at :math:`s = \bar{s}`.

Obviously, if one takes the vectors :math:`B` and :math:`C` randomly, it is
highly possible that they satisfy the two conditions above. Consequently, by
taking the vectors :math:`Y` and :math:`T` supplied by the correction at each 
continuation step and monitoring the sign of :math:`\tau`, a numerical
continuation method is able to detect bifurcation points.

Once a bifurcation point :math:`\bar{Y}` is detected by the sign change in the
test function :math:`\tau`, i.e., :math:`\tau(Y_{j}) \tau(Y_{j+1}) < 0`, it can
be approximated more precisely by the predictor-corrector steps described above
with a special step-length adaptation (see Sect. 8.1 in [all-ge1997]_). In
particular, one can take the subsequent step lengths as 

   .. math::

      h_{j+1} := -\frac{\tau(Y_{j+1})}{\tau(Y_{j+1}) - \tau(Y_{j})}h_{j}

until :math:`\lvert h_{j+1} \rvert < h_{\mathrm{min}}`, which corresponds to the
secant method for finding a zero of :math:`s \mapsto \tau(Y(s))`.

Finally, it would be desirable to switch solution branches. To this end, we
shall consider the case of the so-called *simple bifurcation point*, where only
two distinct solution branches intersect.

Let :math:`\tilde{Y}` be an approximation of :math:`\bar{Y}` that we are given
and :math:`V(\tilde{Y})` be the first part of the solution of the augmented
system for computing the test function :math:`\tau(\tilde{Y})`. As proposed in
[georg2001]_, to obtain a point on the bifurcating (new) branch, one can take
:math:`V(\tilde{Y})` as a predictor direction and do one continuation step
starting with :math:`(\tilde{Y}, V(\tilde{Y}))`. After this continuation step
has been performed successfully and a point on the new branch has been
recovered, one can proceed with the usual predictor-corrector steps to trace
this branch.

Approximation of solution branches of a model
+++++++++++++++++++++++++++++++++++++++++++++

The numerical continuation is defined in ``getfem/getfem_continuation.h``. In
order to use it, one has to do the initialisation first::

  getfem::cont_struct_getfem_model S(model, parameter_name[, initdata_name, finaldata_name, currentdata_name],
                           	     sfac, ls, bifurcations, h_init, h_max, h_min, h_inc, h_dec, maxit, thrit,
				     maxres, maxdiff, mincos, maxres_solve, noisy);
  getfem::init_Moore_Penrose_continuation(S, U, lambda, T_U, T_lambda, h);

where ``parameter_name`` is the name of the model datum representing
:math:`\lambda`, ``sfac`` represents the scale factor :math:`\kappa`, ``ls`` is
the name of the solver to be used for the linear systems incorporated in the
process (e.g., ``getfem::default_linear_solver<getfem::model_real_sparse_matrix, getfem::model_real_plain_vector>(model)``), and the boolean value of
``bifurcations`` determines whether the tools for detection and treatment of
bifurcation points have to be used. The real numbers ``h_init``, ``h_max``,
``h_min``, ``h_inc``, ``h_dec`` denote :math:`h_{\mathrm{init}}`,
:math:`h_{\mathrm{max}}`, :math:`h_{\mathrm{min}}`, :math:`h_{\mathrm{inc}}`,
and :math:`h_{\mathrm{dec}}`, the integers ``maxit`` and ``thrit`` are the
maximum number of iterations allowed in the correction and
:math:`l_{\mathrm{thr}}`, respectively,  ``maxres``, ``maxdiff``, ``mincos``,
and ``maxres_solve`` denote :math:`\varepsilon`, :math:`\varepsilon'`,
:math:`c_{\mathrm{min}}`, and the target residual value for the linear systems
to be solved. Finally, the non-negative integer ``noisy`` determines how
detailed information has to be displayed in the course of the continuation
process (the larger value the more details). Under the optional data names
``initdata_name`` and ``finaldata_name``, :math:`P^{0}` and :math:`P^{1}`
should be stored in the case of the parametrisation by a vector datum,
respectively. Under ``currentdata_name``, the values of :math:`P(\lambda)` are
stored then, that is, actual values of the datum the model depends on. Further,
``U`` should be a solution for the value of parameter :math:`\lambda` equal to
``lambda`` so that :math:`Y_{0}=` (\ ``U``\ ,\ ``lambda``\ ). In accordance with
the sign of the initial value ``T_lambda``, an initial unit tangent
:math:`T_{0}` corresponding to :math:`Y_{0}` is computed and returned in
``T_U``, ``T_lambda``. Moreover, ``h`` is set to the initial step size 
``h_init``.

Consequently, one step of the continuation can be called by ::

  getfem::Moore_Penrose_continuation(S, U, lambda, T_U, T_lambda, h);

After each call, a new point on the solution curve and the corresponding tangent
are returned in the variables ``U``, ``lambda`` and ``T_U``, ``T_lambda``. Step
size to the next prediction is returned in ``h``. It the option ``bifurcations``
has been chosen, the test function for bifurcations is evaluated at the end of
each continuation step. Furthermore, if a bifurcation point is detected, the
procedure for numerical bifurcation is performed and the approximation of the
branching point as well as tangents to both bifurcating branches are saved in
the continuation structure ``S``. From there, they can easily be recovered with
member functions of ``S`` so that one can initialise the continuation to trace
either of the branches next time.

For a complete example of use, see the test programs
``tests/test_continuation.cc``, ``interface/tests/matlab/demo_continuation.m``
or ``interface/src/scilab/demos/demo_continuation.sce``.