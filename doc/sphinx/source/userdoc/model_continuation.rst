.. $Id: model.rst 3655 2010-07-17 20:42:08Z renard $

.. include:: ../replaces.txt

.. highlightlang:: c++

.. index:: models, model bricks

.. _ud-model-continuation:



Numerical continuation
----------------------

Let an algebraic problem coming from discretization of a FEM-model can be
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
implemented (see, for instance, [dh-go-ku2003]_): 

Denoting :math:`Y := (U, \lambda)`, this method consists in computing a sequence
of consecutive points :math:`Y_{j}` on a chosen solution branch and the
corresponding unit tangent vectors :math:`T_{j}`:

.. math::

   F(Y_{j}) = 0,\quad \nabla F(Y_{j}) T_{j} = 0,\quad \lVert T_{j} \rVert = 1,\quad j = 0, 1,\dotsc.

More precisely, let us suppose that we have a couple :math:`(Y_{j}, T_{j})`
satisfying the relations above at our disposal. The next couple is calculated in
two steps -- *prediction* and *correction*.

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

   F^{l}(Y, T) := \begin{pmatrix}F(Y)\\ (T_{j+1}^{l-1})^{T}(Y - Y_{j+1}^{l-1})\\ \nabla F(Y_{j+1}^{l-1})T\\ (T_{j+1}^{l-1})^{T} T - (T_{j+1}^{l-1})^{T}T_{j+1}^{l-1})\end{pmatrix}.

A new couple :math:`(Y_{j+1}, T_{j+1})` is set to
:math:`(Y_{j+1}^{l}, T_{j+1}^{l})` iff
:math:`\lVert F(Y_{j+1}^{l})\rVert \leq \varepsilon`,
:math:`\lVert Y_{j+1}^{l} - Y_{j+1}^{l-1}\rVert \leq \varepsilon'` and
:math:`(T_{j+1}^{l})^{T} T_{j} \geq \vartheta_{\mathrm{min}}`. Let us note that
the partial gradient :math:`\nabla_{U} F` is assembled analytically whereas
:math:`\nabla_{\lambda} F` is evaluated by a finite difference with the
increment equal to :math:`\epsilon`. 

Finally, the step size :math:`h_{j+1}` in the next prediction depends on how
this Newton correction is successfull. Denoting the number of iterations needed
by :math:`l_{\mathrm{it}}`, it is selected as

.. math::

   h_{j+1} := \begin{cases}\max\{h_{\mathrm{dec}} h_{j}, h_{\mathrm{min}}\}& \text{if no new couple was accepted},\\ \min\{h_{\mathrm{inc}} h_{j}, h_{\mathrm{max}}\}& \text{if a new couple was accepted and } l_{\mathrm{it}} < l_{\mathrm{thr}},\\ h_{j}& \text{otherwise},\end{cases}

where :math:`0 < h_{\mathrm{dec}} < 1 < h_{\mathrm{inc}}`,
:math:`0 < l_{\mathrm{thr}}` as well as
:math:`0 < h_{\mathrm{min}} < h_{\mathrm{max}}` are given constants. At the
beginning, one sets :math:`h_{1} := h_{\mathrm{init}}` for some
:math:`h_{\mathrm{min}} \leq h_{\mathrm{init}} \leq h_{\mathrm{max}}`.

In order to apply the Moore-Penrose continuation on a model defined in |gf|, one
has to do the initialization first::

  getfem::S_getfem_model s(model, parameter_name, ls, maxit, thrit, maxres, maxdiff, minang, h_init, h_max, h_min, h_inc, h_dec, eps, maxres_solve, noisy);
  getfem::init_Moore_Penrose_continuation(s, U, lambda, T_U, T_lambda, h);

where ``parameter_name`` is the name of the variable representing the parameter
in the model, ``ls`` is the name of a solver to be used for the linear systems
incorporated in the process (e.g.
``getfem::default_linear_solver<getfem::model_real_sparse_matrix, getfem::model_real_plain_vector>(model)``),
the integers ``maxit`` and ``thrit`` stand for the maximal number of iterations
allowed in the correction and :math:`l_{\mathrm{thr}}`, respectively, the reals
``maxres``, ``maxdiff``, ``minang``, ``h_init``, ``h_max``, ``h_min``,
``h_inc``, ``h_dec``, ``eps``, ``maxres_solve`` denote :math:`\varepsilon`,
:math:`\varepsilon'`, :math:`\vartheta_{\mathrm{min}}`,
:math:`h_{\mathrm{init}}`, :math:`h_{\mathrm{max}}`, :math:`h_{\mathrm{min}}`,
:math:`h_{\mathrm{inc}}`, :math:`h_{\mathrm{dec}}`, :math:`\epsilon` and the
target residual value for the linear systems to be solved, and the non-negative
integer ``noisy`` determines how detailed information has to be displayed in the
course of the continuation process (the larger value the more details). Further,
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

For a complete example of use, you can see the test programs
``tests/test_continuation.cc``, ``interface/tests/matlab/demo_continuation.m``
or ``interface/src/scilab/demos/demo_continuation.sce``.