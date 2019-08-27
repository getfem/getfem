.. $Id$

.. include:: ../replaces.txt

.. highlightlang:: c++

.. _dp-appendixa:

Appendix A. Some basic computations between reference and real elements
=======================================================================

Volume integral
---------------

One has

.. math::

   \int_T f(x)\ dx = \int_{\widehat{T}} \widehat{f}(\widehat{x})
   |\mbox{vol}\left(
   \frac{\partial\tau(\widehat{x})}{\partial \widehat{x}_0};
   \frac{\partial\tau(\widehat{x})}{\partial \widehat{x}_1};
   \ldots;
   \frac{\partial\tau(\widehat{x})}{\partial \widehat{x}_{P-1}}
   \right)|\ d\widehat{x}.


Denoting :math:`J_{\tau}(\widehat{x})` the jacobian

.. math::

   \fbox{$ J_{\tau}(\widehat{x}) :=
   |\mbox{vol}\left(
   \frac{\partial\tau(\widehat{x})}{\partial \widehat{x}_0};
   \frac{\partial\tau(\widehat{x})}{\partial \widehat{x}_1};
   \ldots;
   \frac{\partial\tau(\widehat{x})}{\partial \widehat{x}_{P-1}}
   \right)| =
   (\mbox{det}(K(\widehat{x})^T K(\widehat{x})))^{1/2}$,}

one finally has

.. math::

   \fbox{$\int_T f(x)\ dx = \int_{\widehat{T}} \widehat{f}(\widehat{x}) J_{\tau}(\widehat{x})\ d\widehat{x}$.}

When :math:`P = N`, the expression of the jacobian reduces to :math:`J_{\tau}(\widehat{x})
= |\mbox{det}(K(\widehat{x}))|`.


Surface integral
----------------

With :math:`\Gamma` a part of the boundary of :math:`T` a real element and
:math:`\widehat{\Gamma}` the corresponding boundary on the reference element :math:`\widehat{T}`,
one has

.. math::

   \fbox{$\int_{\Gamma} f(x)\ d\sigma =
   \int_{\widehat{\Gamma}}\widehat{f}(\widehat{x}) \|B(\widehat{x})\widehat{n}\| J_{\tau}(\widehat{x})\ d\widehat{\sigma}$,}

where :math:`\widehat{n}` is the unit normal to :math:`\widehat{T}` on :math:`\widehat{\Gamma}`. In a same
way

.. math::

   \fbox{$\int_{\Gamma} F(x)\cdot n\ d\sigma =
   \int_{\widehat{\Gamma}} \widehat{F}(\widehat{x})\cdot(B(\widehat{x})\cdot\widehat{n}) J_{\tau}(\widehat{x})\ d\widehat{\sigma}$,}

For :math:`n` the unit normal to :math:`T` on :math:`\Gamma`.


Derivative computation
----------------------

One has

.. math::

   \nabla f(x) = B(\widehat{x})\widehat{\nabla} \widehat{f}(\widehat{x}).


Second derivative computation
-----------------------------

Denoting

.. math::

   \nabla^2 f =
   \left[\frac{\partial^2 f}{\partial x_i \partial x_j}\right]_{ij},

the :math:`N \times N` matrix and

.. math::

   \widehat{X}(\widehat{x}) =
   \sum_{k = 0}^{N-1}\widehat{\nabla}^2\tau_k(\widehat{x})\frac{\partial f}{\partial x_k}(x) =
   \sum_{k = 0}^{N-1}\sum_{i = 0}^{P-1}
   \widehat{\nabla}^2\tau_k(\widehat{x})B_{ki}\frac{\partial \widehat{f}}{\partial \widehat{x}_i}(\widehat{x}),

the :math:`P \times P` matrix, then

.. math::

   \widehat{\nabla}^2 \widehat{f}(\widehat{x}) = \widehat{X}(\widehat{x}) + K(\widehat{x})^T \nabla^2 f(x) K(\widehat{x}),

and thus

.. math::

   \nabla^2 f(x) = B(\widehat{x})(\widehat{\nabla}^2 \widehat{f}(\widehat{x}) - \widehat{X}(\widehat{x})) B(\widehat{x})^T.

In order to have uniform methods for the computation of elementary matrices, the
Hessian is computed as a column vector :math:`H f` whose components are
:math:`\frac{\partial^2 f}{\partial x^2_0}, \frac{\partial^2 f}{\partial
x_1\partial x_0},\ldots, \frac{\partial^2 f}{\partial x^2_{N-1}}`. Then, with
:math:`B_2` the :math:`P^2 \times P` matrix defined as

.. math::

   \left[B_2(\widehat{x})\right]_{ij} =
   \sum_{k = 0}^{N-1}
   \frac{\partial^2 \tau_k(\widehat{x})}{\partial \widehat{x}_{i / P} \partial \widehat{x}_{i\mbox{ mod }P}}
   B_{kj}(\widehat{x}),

and :math:`B_3` the :math:`N^2 \times P^2` matrix defined as

.. math::

   \left[B_3(\widehat{x})\right]_{ij} =
   B_{i / N, j / P}(\widehat{x}) B_{i\mbox{ mod }N, j\mbox{ mod }P}(\widehat{x}),

one has

.. math::

   \fbox{$H f(x) = B_3(\widehat{x})
   \left(\widehat{H}\ \widehat{f}(\widehat{x}) - B_2(\widehat{x})\widehat{\nabla} \widehat{f}(\widehat{x})\right)$.}


Example of elementary matrix
----------------------------

Assume one needs to compute the elementary "matrix":

.. math::

   t(i_0, i_1, \ldots, i_7) =
   \int_{T}\varphi_{i_1}^{i_0}
   \partial_{i_4}\varphi_{i_3}^{i_2}
   \partial^2_{i_7/ P, i_7\mbox{ mod } P}\varphi_{i_6}^{i_5}\ dx,

The computations to be made on the reference elements are

.. math::

   \widehat{t}_0(i_0, i_1, \ldots,i_7) =
   \int_{\widehat{T}}(\widehat{\varphi})_{i_1}^{i_0}
   \partial_{i_4}(\widehat{\varphi})_{i_3}^{i_2}
   \partial^2_{i_7 / P, i_7\mbox{ mod } P}(\widehat{\varphi})_{i_6}^{i_5} J(\widehat{x})\ d\widehat{x},

and

.. math::

   \widehat{t}_1(i_0, i_1, \ldots, i_7) =
   \int_{\widehat{T}}(\widehat{\varphi})_{i_1}^{i_0}
   \partial_{i_4}(\widehat{\varphi})_{i_3}^{i_2}
   \partial_{i_7}(\widehat{\varphi})_{i_6}^{i_5} J(\widehat{x})\ d\widehat{x},

Those two tensor can be computed once on the whole reference element if the
geometric transformation is linear (because :math:`J(\widehat{x})` is constant). If the
geometric transformation is non-linear, what has to be stored is the value on
each integration point. To compute the integral on the real element a certain
number of reductions have to be made:

* Concerning the first term (:math:`\varphi_{i_1}^{i_0}`) nothing.

* Concerning the second term (:math:`\partial_{i_4}\varphi_{i_3}^{i_2}`) a
  reduction with respect to :math:`i_4` with the matrix :math:`B`.

* Concerning the third term (:math:`\partial^2_{i_7 / P, i_7\mbox{ mod }P}
  \varphi_{i_6}^{i_5}`) a reduction of :math:`\widehat{t}_0` with respect to :math:`i_7`
  with the matrix :math:`B_3` and a reduction of :math:`\widehat{t}_1` with respect also
  to :math:`i_7` with the matrix :math:`B_3 B_2`


The reductions are to be made on each integration point if the geometric
transformation is non-linear. Once those reductions are done, an addition of all
the tensor resulting of those reductions is made (with a factor equal to the load
of each integration point if the geometric transformation is non-linear).

If the finite element is non-:math:`\tau`-equivalent, a supplementary reduction of the
resulting tensor with the matrix :math:`M` has to be made.
