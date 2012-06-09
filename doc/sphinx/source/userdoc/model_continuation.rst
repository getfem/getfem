.. $Id: model.rst 3655 2010-07-17 20:42:08Z renard $

.. include:: ../replaces.txt

.. highlightlang:: c++

.. index:: models, model bricks

.. _ud-model-continuation:



Numerical continuation
----------------------

Let an algebraic problem coming from a discretization of a FEM-model can be
written in the form

.. math::

   F(U) = 0.

In what follows, we shall suppose that the model depends on an additional scalar
parameter :math:`\lambda` so that :math:`F(U) = F(U, \lambda)`. Our aim is to
find the solution set of the problem

.. math::

   F(U, \lambda) = 0.

This can be done by using *continuation* (*path-following*) techniques. In |gf|,
the (approximate) *Moore-Penrose* (also called *Gauss-Newton*) continuation is
implemented (see, for instance, [dh-go-ku2003]_).

Since this method does not make an explicit difference between the variable
:math:`U` and the parameter :math:`\lambda`, we shall denote
:math:`Y := (U, \lambda)` for brevity. Nevertheless, to avoid bad scaling of the
values of the continuation parameter, we shall use the following weighted scalar
product and norm:

.. math::

   \langle Y, \tilde{Y} \rangle := \kappa \langle U, \tilde{U} \rangle + \lambda \tilde{\lambda},\quad \lVert Y \rVert := \sqrt{\kappa \lVert U \rVert^{2} + \lambda^{2}},\qquad Y = (U, \lambda),\, \tilde{Y} = (\tilde{U}, \tilde{\lambda}).

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

   F(Y_{j}) = 0,\quad \nabla F(Y_{j}) T_{j} = 0,\quad \lVert T_{j} \rVert = 1,\quad j = 0, 1,\dotsc.

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
:math:`T_{j+1}^{l} := \tilde{T}_{j+1}^{l} / \lVert \tilde{T}_{j+1}^{l} \rVert`
and the couple :math:`(Y_{j+1}^{l}, \tilde{T}_{j+1}^{l})` is given by one
iteration of the Newton method applied to the equation :math:`F^{l}(Y, T) = 0`
with

.. math::

   F^{l}(Y, T) := \begin{pmatrix}F(Y)\\ \langle T_{j+1}^{l-1}, Y - Y_{j+1}^{l-1} \rangle\\ \nabla F(Y_{j+1}^{l-1})T\\ \langle T_{j+1}^{l-1}, T \rangle - \langle T_{j+1}^{l-1}, T_{j+1}^{l-1} \rangle\end{pmatrix}.

.. _ud_fig_correction:
.. figure:: images/getfemusercorrection.png
   :align: center

   Correction

A new couple :math:`(Y_{j+1}, T_{j+1})` is set to
:math:`(Y_{j+1}^{l}, T_{j+1}^{l})` iff
:math:`\lVert F(Y_{j+1}^{l})\rVert \leq \varepsilon`,
:math:`\lVert Y_{j+1}^{l} - Y_{j+1}^{l-1}\rVert \leq \varepsilon'` and
:math:`\langle T_{j+1}^{l}, T_{j} \rangle \geq \vartheta_{\mathrm{min}}`. Let us
note that the partial gradient :math:`\nabla_{U} F` is assembled analytically
whereas :math:`\nabla_{\lambda} F` is evaluated by a finite difference with an
increment equal to :math:`\epsilon`. 

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

1. The parameter :math:`\lambda` is directly some scalar datum that the model
   depends on.

2. The model is parametrised by the scalar parameter :math:`\lambda` *via* some
   vector datum :math:`V` that the model depends on. In this case, one takes the
   linear path

   .. math::

      \lambda \mapsto V(\lambda) := (1 - \lambda)V^{0} + \lambda V^{1},

   where :math:`V^{0}` and :math:`V^{1}` are given values of :math:`V`, and one
   traces the solution set of the problem

   .. math::

      F(U, V(\lambda)) = 0.

In order to use the continuation, one has to do its initialisation first::

  getfem::S_getfem_model s(model, parameter_name[, initdata_name, finaldata_name, currentdata_name],
                           ls, sfac, maxit, thrit, maxres, maxdiff, minang, h_init, h_max, h_min,
			   h_inc, h_dec, eps, maxres_solve, noisy);
  getfem::init_Moore_Penrose_continuation(s, U, lambda, T_U, T_lambda, h);

where ``parameter_name`` is the name of the model datum representing
:math:`\lambda`, ``ls`` is the name of the solver to be used for the linear
systems incorporated in the process (e.g.
``getfem::default_linear_solver<getfem::model_real_sparse_matrix, getfem::model_real_plain_vector>(model)``),
``sfac`` represents the scale factor :math:`\kappa`, the integers ``maxit`` and
``thrit`` stand for the maximal number of iterations allowed in the correction
and :math:`l_{\mathrm{thr}}`, respectively, the real numbers ``maxres``,
``maxdiff``, ``minang``, ``h_init``, ``h_max``, ``h_min``, ``h_inc``, ``h_dec``,
``eps``, ``maxres_solve`` denote :math:`\varepsilon`, :math:`\varepsilon'`,
:math:`\vartheta_{\mathrm{min}}`, :math:`h_{\mathrm{init}}`,
:math:`h_{\mathrm{max}}`, :math:`h_{\mathrm{min}}`, :math:`h_{\mathrm{inc}}`,
:math:`h_{\mathrm{dec}}`, :math:`\epsilon` and the target residual value for the
linear systems to be solved, and the non-negative integer ``noisy`` determines
how detailed information has to be displayed in the course of the continuation
process (the larger value the more details). Under the optional data names
``initdata_name`` and ``finaldata_name``, :math:`V^{0}` and :math:`V^{1}`
should be stored in the case of the parametrisation by vector datum,
respectively. Under ``currentdata_name``, the values of :math:`V(\lambda)` are
stored then, that is, actual values of the datum the model depends on. Further,
``U`` should be a solution for the value of parameter :math:`\lambda` equal to
``lambda`` so that :math:`Y_{0}=` (\ ``U``\ ,\ ``lambda``\ ). In accordance with
the sign of the initial value ``T_lambda``, an initial unit tangent
:math:`T_{0}` corresponding to :math:`Y_{0}` is computed and returned in
``T_U``, ``T_lambda``. Moreover, ``h`` is set to the initial step size
``h_init``.

Consequently, one step of the continuation can be called by ::

  getfem::Moore_Penrose_continuation(s, U, lambda, T_U, T_lambda, h);

After each call, a new point on the solution curve and the corresponding tangent
are returned in the variables ``U``, ``lambda`` and ``T_U``, ``T_lambda``. Step
size to the next prediction is returned in ``h``.

For a complete example of use, see the test programs
``tests/test_continuation.cc``, ``interface/tests/matlab/demo_continuation.m``
or ``interface/src/scilab/demos/demo_continuation.sce``.