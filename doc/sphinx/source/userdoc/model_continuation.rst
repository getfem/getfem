.. $Id: model.rst 3655 2010-07-17 20:42:08Z renard $

.. include:: ../replaces.txt

.. highlightlang:: c++

.. index:: models, model bricks

.. _ud-model-continuation:



Numerical continuation and bifurcation
--------------------------------------

Let an algebraic problem coming from discretisation of a FEM-model can be
written in the form

.. math::

   F(U) = 0.

In what follows, we shall suppose that the model depends on an additional scalar
parameter :math:`\lambda` so that :math:`F(U) = F(U, \lambda)`. 

Numerical continuation
++++++++++++++++++++++

Methods of numerical continuation serve for tracing solutions of the system

.. math::

   F(U, \lambda) = 0, \quad F\colon \mathbb{R}^{N} \times \mathbb{R} \to \mathbb{R}^{N}.

In |gf|, a continuation technique for piecewise :math:`C^{1}` (:math:`PC^{1}`)
solution curves is implemented (see [li-re]_ for more details). Since it does
not make an explicit difference between the state variable :math:`U` and the
parameter :math:`\lambda`, we shall denote :math:`Y := (U, \lambda)` for
brevity. Nevertheless, to avoid bad scaling when calculating tangents, for
example, we shall use the following weighted scalar product and norm:

.. math::

   \langle Y, \tilde{Y} \rangle_{w} := \kappa \langle U, \tilde{U} \rangle + \lambda \tilde{\lambda},\quad \lVert Y \rVert_{w} := \sqrt{\kappa \lVert U \rVert^{2} + \lambda^{2}},\qquad Y = (U, \lambda),\, \tilde{Y} = (\tilde{U}, \tilde{\lambda}).

Here, :math:`\kappa` should be chosen so that
:math:`\kappa \langle U, \tilde{U} \rangle` is proportional to the scalar
product of the corresponding space variables, usually in :math:`L^{2}`. One can
take, for example, :math:`\kappa = h^{d}`, where :math:`h` is the mesh size and
:math:`d` stands for the dimension of the underlying problem. Alternatively,
:math:`\kappa` can be chosen as :math:`1/N` for simplicity.

The idea of the continuation strategy is to continue smooth pieces of solution
curves by a classical predictor-corrector method and to join the smooth pieces
continuously.

The particular predictor-corrector method employed is a slight modification of
the *inexact Moore-Penrose* continuation implemented in MATCONT [dh-go-ku2003]_.
It computes a sequence of consecutive points :math:`Y_{j}` lying approximately
on a solution curve and a sequence of the corresponding unit tangent vectors
:math:`T_{j}`:

.. math::

   \lVert F(Y_{j}) \rVert \leq \varepsilon,\quad F'(Y_{j}; T_{j}) = 0,\quad \lVert T_{j} \rVert_{w} = 1,\quad j = 0, 1,\dotsc.

To describe it, let us suppose that we have a couple :math:`(Y_{j}, T_{j})`
satisfying the relations above at our disposal. In the *prediction*, an initial
approximation of :math:`(Y_{j+1}, T_{j+1})` is taken as

.. math::

   Y_{j+1}^{0} := Y_{j} + h_{j} T_{j},\quad T_{j+1}^{0} := T_{j},

where :math:`h_{j}` is a step size. Its choice will be discussed later on.

In the *correction*, one computes a sequence
:math:`\{(Y_{j+1}^{l}, T_{j+1}^{l})\}`, where
:math:`T_{j+1}^{l} := \tilde{T}_{j+1}^{l} / \lVert \tilde{T}_{j+1}^{l} \rVert_{w}`
and the couple :math:`(Y_{j+1}^{l}, \tilde{T}_{j+1}^{l})` is given by one
iteration of the Newton method applied to the equation :math:`F^{l}(Y, T) = 0`
with

.. math::

   F^{l}(Y, T) := \begin{pmatrix}F(Y)\\ (T_{j+1}^{l-1})^{\top}(Y - Y_{j+1}^{l-1})\\ \nabla F(Y_{j+1}^{l-1})T\\ \langle T_{j+1}^{l-1}, T \rangle_{w} - \langle T_{j+1}^{l-1}, T_{j+1}^{l-1} \rangle_{w}\end{pmatrix}

and the initial approximation :math:`(Y_{j+1}^{l-1}, T_{j+1}^{l-1})`. Due to
potential non-differentiability of :math:`F`, a piecewise-smooth variant of the
Newton method is used (7.2.14 Algorithm in [fa-pa2003]_).

.. _ud_fig_correction:
.. figure:: images/getfemusercorrection.png
   :align: center

   Correction.

A couple :math:`(Y_{j+1}^{l}, T_{j+1}^{l})` is accepted for
:math:`(Y_{j+1}, T_{j+1})` if
:math:`\lVert F(Y_{j+1}^{l})\rVert \leq \varepsilon`,
:math:`\lVert Y_{j+1}^{l} - Y_{j+1}^{l-1}\rVert_{w} \leq \varepsilon'`, and the
cosine of the angle between :math:`T_{j+1}^{l}` and :math:`T_{j}` is greater or
equal to :math:`c_{\mathrm{min}}`. Let us note that the partial gradient of
:math:`F` (or of one of its selection functions in the case of
non-differentiability) with respect to :math:`U` is assembled analytically
whereas the partial gradient with respect to :math:`\lambda` is evaluated by
forward finite differences with the increment equal to 1e-8.

The step size :math:`h_{j+1}` in the next prediction depends on how the Newton
correction has been successful. Denoting the number of iterations needed by
:math:`l_{\mathrm{it}}`, it is selected as

.. math::

   h_{j+1} := \begin{cases}\max\{h_{\mathrm{dec}} h_{j}, h_{\mathrm{min}}\}& \text{if no new couple was accepted},\\ \min\{h_{\mathrm{inc}} h_{j}, h_{\mathrm{max}}\}& \text{if a new couple was accepted and } l_{\mathrm{it}} < l_{\mathrm{thr}},\\ h_{j}& \text{otherwise},\end{cases}

where :math:`0 < h_{\mathrm{dec}} < 1 < h_{\mathrm{inc}}`,
:math:`0 < l_{\mathrm{thr}}` and
:math:`0 < h_{\mathrm{min}} < h_{\mathrm{max}}` are given constants. At the
beginning, one sets :math:`h_{1} := h_{\mathrm{init}}` for some
:math:`h_{\mathrm{min}} \leq h_{\mathrm{init}} \leq h_{\mathrm{max}}`.

Now, let us suppose that we have approximated a piece of a solution curve
corresponding to one sub-domain of smooth behaviour of :math:`F` and we want to
recover a piece corresponding to another sub-domain of smooth behaviour. Let
:math:`(Y_{j},T_{j})` be the last computed couple.

.. _ud_fig_transition:
.. figure:: images/getfemusertransition.png
   :align: center

   Transition between smooth pieces of a solution curve.

To approximate the tangent to the other smooth piece, we first take a point
:math:`Y_{j} + h T_{j}` with :math:`h` a bit greater than
:math:`h_{\mathrm{min}}` so that this point belongs to the interior of the other
sub-domain of smooth behaviour. Then, we find :math:`\tilde{T}` such that

.. math::

   \nabla F(Y_{j} + h T_{j}) \tilde{T} = 0,\quad \lVert \tilde{T} \rVert_{w} = 1

and it remains to determine an appropriate direction of this vector. This can be
done on the basis of the following observations:  First, there exists
:math:`r \in \{\pm 1\}` such that :math:`Y_{j} - r \tilde{h} \tilde{T}` remains
in the same sub-domain as :math:`Y_{j}` for any :math:`\tilde{h}` positive.
This can be characterised by the fact that
:math:`\frac{\lvert T_{-}^{\top} \tilde{T}\rvert}{\lVert T_{-} \rVert \lVert \tilde{T} \rVert}`
is significantly smaller than 1 for :math:`T_{-}` with
:math:`\nabla F(Y_{j} - r \tilde{h} \tilde{T}) T_{-} = 0`. Second,
:math:`Y_{j} + r \tilde{h} \tilde{T}` appears in the other sub-domain for
:math:`\tilde{h}` larger than some positive threshold and for such values,
:math:`\frac{\lvert T_{+}^{\top} \tilde{T}\rvert}{\lVert T_{+} \rVert \lVert \tilde{T} \rVert}`
is close to 1 for :math:`T_{+}` with
:math:`\nabla F(Y_{j} + r \tilde{h} \tilde{T}) T_{+} = 0`.

This suggests the following procedure for selecting the desired direction of
:math:`\tilde{T}`: Increase the values of :math:`\tilde{h}` successively from
:math:`h_{\mathrm{min}}` and when you arrive at :math:`\tilde{h}` and
:math:`r \in \{\pm 1\}` such that

.. math::

   \frac{\lvert T^{\top} \tilde{T}\rvert}{\lVert T \rVert \lVert \tilde{T} \rVert} \approx 1\quad \text{if}\ \nabla F(Y_{j} + r \tilde{h} \tilde{T}) T = 0,

take :math:`r \tilde{T}` as an approximation of the tangent to the other smooth
piece.

After finding a new tangent, say :math:`\tilde{T}_{j}`, we restart the
predictor-corrector with :math:`(Y_{j}, \tilde{T}_{j})`.

In |gf|, the continuation is implemented for two ways of parametrisation of the
model:

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

Approximation of solution curves of a model
+++++++++++++++++++++++++++++++++++++++++++++

The numerical continuation is defined in ``getfem/getfem_continuation.h``. In
order to use it, one has to do the following initialisation first::

  getfem::cont_struct_getfem_model S(model, parameter_name[, initdata_name, finaldata_name, currentdata_name],
                           	     sfac, ls, bifurcations, h_init, h_max, h_min, h_inc, h_dec, maxit, thrit,
				     maxres, maxdiff, mincos, maxres_solve, noisy, non-smooth);
  getfem::init_Moore_Penrose_continuation(S, U, lambda, T_U, T_lambda, h);

where ``parameter_name`` is the name of the model datum representing
:math:`\lambda`, ``sfac`` represents the scale factor :math:`\kappa`, ``ls`` is
the name of the solver to be used for the linear systems incorporated in the
process (e.g., ``getfem::default_linear_solver<getfem::model_real_sparse_matrix, getfem::model_real_plain_vector>(model)``), and the boolean value of
``bifurcations`` determines whether the tools for detection and treatment of
bifurcation points have to be used. The real numbers ``h_init``, ``h_max``,
``h_min``, ``h_inc``, ``h_dec`` denote :math:`h_{\mathrm{init}}`,
:math:`h_{\mathrm{max}}`, :math:`h_{\mathrm{min}}`, :math:`h_{\mathrm{inc}}`,
and :math:`h_{\mathrm{dec}}`, the integer ``maxit`` is the maximum number of
iterations allowed in the correction and ``thrit``, ``maxres``, ``maxdiff``,
``mincos``, and ``maxres_solve`` denote :math:`l_{\mathrm{thr}}`,
:math:`\varepsilon`, :math:`\varepsilon'`, :math:`c_{\mathrm{min}}`, and the
target residual value for the linear systems to be solved, respectively. The
non-negative integer ``noisy`` determines how detailed information has to be
displayed in the course of the continuation process (the larger value the more
details) and the boolean value of ``non-smooth`` determines whether only the
predictor-corrector for smooth curves has to be used or the technique for
searching for other smooth pieces of solution curves has to be employed as well.
Under the optional data names ``initdata_name`` and ``finaldata_name``,
:math:`P^{0}` and :math:`P^{1}` should be stored, respectively, in the case of
the parametrisation by a vector datum. Under ``currentdata_name``, the values of
:math:`P(\lambda)` are stored then, that is, actual values of the datum the
model depends on. Further, ``U`` should be a solution for the value of parameter
:math:`\lambda` equal to ``lambda`` so that
:math:`Y_{0}=` (\ ``U``\ ,\ ``lambda``\ ). In accordance with the sign of the
initial value ``T_lambda``, an initial unit tangent :math:`T_{0}` corresponding
to :math:`Y_{0}` is computed and returned in ``T_U``, ``T_lambda``. Moreover,
``h`` is set to the initial step size ``h_init``.

Consequently, one step of the continuation can be called by ::

  getfem::Moore_Penrose_continuation(S, U, lambda, T_U, T_lambda, h);

After each call, a new point on a solution curve and the corresponding tangent
are returned in the variables ``U``, ``lambda`` and ``T_U``, ``T_lambda``. The
step size for the next prediction is returned in ``h``. If the option
``bifurcations`` has been chosen, the test function for bifurcations is
evaluated at the end of each continuation step. Furthermore, if a bifurcation
point is detected, the procedure for numerical bifurcation is performed and an
approximation of the branching point as well as tangents to both bifurcating
branches are saved in the continuation structure ``S``. From there, they can
easily be recovered with member functions of ``S`` so that one can initialise
the continuation to trace either of the branches next time.

Complete examples of use are shown in the test programs
``tests/test_continuation.cc``, ``interface/tests/matlab/demo_continuation.m``
and ``interface/src/scilab/demos/demo_continuation.sce`` with the
predictor-corrector only and in
``interface/src/scilab/demos/demo_continuation_block.sce`` or
``interface/src/scilab/demos/demo_continuation_vee.sce`` with the adaptation to
non-smooth problems described above.