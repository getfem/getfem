.. $Id$

.. include:: ../replaces.txt

.. highlight:: matlab

.. _mlab-examples:

Examples
========

.. _mlab-laplacianexample:

A step-by-step basic example
----------------------------

This example shows the basic usage of getfem, on the über-canonical problem above
all others: solving the :envvar:`Laplacian`, :math:`-\Delta u = f` on a square,
with the Dirichlet condition :math:`u = g(x)` on the domain boundary. You can find
the **m-file** of this example under the name **demo_step_by_step.m** in the
directory ``interface/tests/matlab-octave/`` of the |gf| distribution.

The first step is to **create a mesh**. It is possible to create simple structured meshes or unstructured meshes for simple geometries (see ``gf_mesh('generate', mesher_object mo, scalar h)``) or to rely on an external mesher (see ``gf_mesh('import', string
FORMAT, string FILENAME))``).  For this example, we
just consider a regular **cartesian mesh** whose nodes are
:math:`\{x_{i=0\ldots10,j=0..10}=(i/10,j/10)\}`::

  >> % creation of a simple cartesian mesh
  >> m = gf_mesh('cartesian',[0:.1:1],[0:.1:1]);
  m =
       id: 0
      cid: 0

If you try to look at the value of ``m``, you'll notice that it appears to be a
structure containing two integers. The first one is its identifier, the second one
is its class-id, i.e. an identifier of its type. This small structure is just an
"handle" or "descriptor" to the real object, which is stored in the |gf| memory
and cannot be represented via |octv| and |mlab| data structures. Anyway, you can still inspect the |gf| objects via the command ``gf_workspace('stats')``.

Now we can try to have a **look at the mesh**, with its vertices numbering and the
convexes numbering::

  >> % we enable vertices and convexes labels
  >> gf_plot_mesh(m, 'vertices', 'on', 'convexes', 'on');

As you can see, the mesh is regular, and the numbering of its nodes and convexes
is also regular (this is guaranteed for cartesian meshes, but do not hope a
similar numbering for the degrees of freedom).

The next step is to **create a mesh_fem object**. This one links a mesh with a set
of FEM::

  >> % create a mesh_fem of for a field of dimension 1 (i.e. a scalar field)
  >> mf = gf_mesh_fem(m,1);
  >> gf_mesh_fem_set(mf,'fem',gf_fem('FEM_QK(2,2)'));

The first instruction builds a new |mlab_mf| object, the second argument specifies
that this object will be used to interpolate scalar fields (since the unknown
:math:`u` is a scalar field). The second instruction assigns the :math:`Q_2` FEM
to every convex (each basis function is a polynomial of degree 4, remember that
:math:`P_k` are polynomials of degree :math:`k`, while
:math:`Q_k` are polynomials of degree :math:`2k`). As :math:`Q_2` is a
polynomial FEM, you can view the expression of its basis functions on the
reference element::

  >> gf_fem_get(gf_fem('FEM_QK(2,2)'), 'poly_str');
  ans =
      '1 - 3*x - 3*y + 2*x^2 + 9*x*y + 2*y^2 - 6*x^2*y - 6*x*y^2 + 4*x^2*y^2'
      '4*x - 4*x^2 - 12*x*y + 12*x^2*y + 8*x*y^2 - 8*x^2*y^2'
      '-x + 2*x^2 + 3*x*y - 6*x^2*y - 2*x*y^2 + 4*x^2*y^2'
      '4*y - 12*x*y - 4*y^2 + 8*x^2*y + 12*x*y^2 - 8*x^2*y^2'
      '16*x*y - 16*x^2*y - 16*x*y^2 + 16*x^2*y^2'
      '-4*x*y + 8*x^2*y + 4*x*y^2 - 8*x^2*y^2'
      '-y + 3*x*y + 2*y^2 - 2*x^2*y - 6*x*y^2 + 4*x^2*y^2'
      '-4*x*y + 4*x^2*y + 8*x*y^2 - 8*x^2*y^2'
      'x*y - 2*x^2*y - 2*x*y^2 + 4*x^2*y^2'

It is also possible to make use of the "object oriented" features of |octv| and |mlab|. As
you may have noticed, when a class "foo" is provided by the |gfi|, it is build
with the function ``gf_foo``, and manipulated with the functions ``gf_foo_get``
and ``gf_foo_set``. But you may also create the
object with the ``gfFoo`` constructor , and manipulated with the ``get(..)`` and
``set(..)`` methods. For example, the previous steps could have been::

  >> gfFem('FEM_QK(2,2)');
  gfFem object ID=0 dim=2, target_dim=1, nbdof=9,[EQUIV, POLY, LAGR], est.degree=4
    -> FEM_QK(2,2)
  >> m=gfMesh('cartesian', [0:.1:1], [0:.1:1]);
  gfMesh object ID=0 [16512 bytes], dim=2, nbpts=121, nbcvs=100
  >> mf=gfMeshFem(m,1);
  gfMeshFem object: ID=1 [804 bytes], qdim=1, nbdof=0,
    linked gfMesh object: dim=2, nbpts=121, nbcvs=100
  >> set(mf, 'fem', gfFem('FEM_QK(2,2)'));
  >> mf
  gfMeshFem object: ID=1 [1316 bytes], qdim=1, nbdof=441,
    linked gfMesh object: dim=2, nbpts=121, nbcvs=100

Now, in order to perform numerical integrations on ``mf``, we need to **build a
mesh_im object**::

  >> % assign the same integration method on all convexes
  >> mim = gf_mesh_im(m, gf_integ('IM_EXACT_PARALLELEPIPED(2)'));

The integration method will be used to compute the various integrals on each
element: here we choose to perform exact computations (no :envvar:`quadrature
formula`), which is possible since the geometric transformation of these convexes
from the reference convex is linear (this is true for all simplices, and this is
also true for the parallelepipeds of our regular mesh, but it is not true for
general quadrangles), and the chosen FEM is polynomial. Hence it is possible to
analytically integrate every basis function/product of basis
functions/gradients/etc. There are many alternative FEM methods and integration
methods (see :ref:`ud`).

Note however that in the general case, approximate integration methods are a
better choice than exact integration methods.

Now we have to **find the** ":envvar:`boundary`" **of the domain**, in order to
set a Dirichlet condition. A mesh object has the ability to store some sets of
convexes and convex faces. These sets (called "regions") are accessed via an
integer #id::

  >> % detect the border of the mesh
  >> border = gf_mesh_get(m,'outer faces');
  >> % mark it as boundary #42
  >> gf_mesh_set(m, 'region', 42, border);
  >> gf_plot_mesh(m, 'regions', [42]); % the boundary edges appears in red

Here we find the faces of the convexes which are on the boundary of the mesh (i.e.
the faces which are not shared by two convexes).

 Remark:

   we could have used ``gf_mesh_get(m, 'OuTEr_faCes')``, as the interface is
   case-insensitive, and whitespaces can be replaced by underscores.

The array ``border`` has two rows, on the first row is a convex number, on the
second row is a face number (which is local to the convex, there is no global
numbering of faces). Then this set of faces is assigned to the region number 42.

At this point, we just have to describe the model and run the solver to get the
solution! The ":envvar:`model`" is created with the ``gf_model`` (or ``gfModel``)
constructor. A model is basically an object which build a global linear system
(tangent matrix for non-linear problems) and its associated right hand side.
Typical modifications are insertion of the stiffness matrix for the problem
considered (linear elasticity, laplacian, etc), handling of a set of contraints,
Dirichlet condition, addition of a source term to the right hand side etc. The
global tangent matrix and its right hand side are stored in the ":envvar:`model`"
structure.

Let us build a problem with an easy solution: :math:`u=x(x-1)y(y-1)+x^5`, then we
have :math:`\Delta u=2(x^2+y^2)-2(x+y)+20x^3` (the FEM won't be able to catch the
exact solution since we use a :math:`Q^2` method).

We start with an empty real model::

  >> % empty real model
  >> md = gf_model('real');

(a model is either ``'real'`` or ``'complex'``). And we declare that ``u`` is an
unknown of the system on the finite element method `mf` by::

  >> % declare that "u" is an unknown of the system
  >> % on the finite element method `mf`
  >> gf_model_set(md, 'add fem variable', 'u', mf);

Now, we add a "generic elliptic" brick, which handles :math:`-\nabla\cdot(A:\nabla
u) = \ldots` problems, where :math:`A` can be a scalar field, a matrix field, or
an order 4 tensor field. By default, :math:`A=1`. We add it on our main variable
``u`` with::

  >> % add generic elliptic brick on "u"
  >> gf_model_set(md, 'add Laplacian brick', mim, 'u');


Next we add a Dirichlet condition on the domain boundary::

  >> % add Dirichlet condition
  >> Uexact = gf_mesh_fem_get(mf, 'eval', {'(x-.5).^2 + (y-.5).^2 + x/5 - y/3'});
  >> gf_model_set(md, 'add initialized fem data', 'DirichletData', mf, Uexact);
  >> gf_model_set(md, 'add Dirichlet condition with multipliers', mim, 'u', mf, 42, 'DirichletData');


The two first lines defines a data of the model which represents the value of the
Dirichlet condition. The third one add a Dirichlet condition to the variable ``u``
on the boundary number ``42``. The dirichlet condition is imposed with lagrange
multipliers. Another possibility is to use a penalization. A |mlab_mf| argument is
also required, as the Dirichlet condition :math:`u=g` is imposed in a weak form
:math:`\int_\Gamma u(x)v(x) = \int_\Gamma g(x)v(x) ~ \forall v` where :math:`v` is
taken in the space of multipliers given by here by ``mf``.


.. topic:: Remark:

   the polynomial expression was interpolated on ``mf``. It is possible only if
   ``mf`` is of Lagrange type. In this first example we use the same |mlab_mf| for
   the unknown and for the data such as ``g``, but in the general case, ``mf``
   won't be Lagrangian and another (Lagrangian) |mf| will be used for the
   description of Dirichlet conditions, source terms etc.

A source term can be added with the following lines::

  >> % add source term
  >> f = gf_mesh_fem_get(mf, 'eval', { '2(x^2+y^2)-2(x+y)+20x^3' });
  >> gf_model_set(md, 'add initialized fem data', 'VolumicData', mf, f);
  >> gf_model_set(md, 'add source term brick', mim, 'u', 'VolumicData');

It only remains now to launch the solver. The linear system is assembled and solve
with the instruction::

  >> % solve the linear system
  >> gf_model_get(md, 'solve');

The model now contains the solution (as well as other things, such as the linear
system which was solved). It is extracted, a display into a |octv| or |mlab| figure::

  >> % extracted solution
  >> u = gf_model_get(md, 'variable', 'u');
  >> % display
  >> gf_plot(mf, u, 'mesh','on');


Another Laplacian with exact solution
-------------------------------------

This is the :file:`tests/matlab-octave/demo_laplacian.m` example.

.. literalinclude:: code_samples/demo_laplacian.m


Linear and non-linear elasticity
--------------------------------

This example uses a mesh that was generated with `GiD`_. The object is meshed
with quadratic tetrahedrons. You can find the ``m-file`` of this example under
the name :file:`demo_tripod.m` in the directory :file:`tests/matlab-octave` of the
toolbox distribution.

.. literalinclude:: code_samples/demo_tripod.m

Here is the final figure, displaying the :envvar:`Von Mises` stress:

.. _malb-fig-tripod-vm:
.. figure:: images/tripodvonmiseswithmesh.png
   :width: 300pt
   :align: center

   deformed tripod


Avoiding the bricks framework
-----------------------------

The model bricks are very convenient, as they hide most of the details of the
assembly of the final linear systems. However it is also possible to stay at a
lower level, and handle the assembly of linear systems, and their resolution,
directly in |octv| or |mlab|. For example, the demonstration :file:`demo_tripod_alt.m` is
very similar to the :file:`demo_tripod.m` except that the assembly is explicit::

  nbd=get(mfd, 'nbdof');
  F = gf_asm('boundary_source', 1, mim, mfu, mfd, repmat([0;-10;0],1,nbd));
  K = gf_asm('linear_elasticity', mim, mfu, mfd, ...
             lambda*ones(1,nbd),mu*ones(1,nbd));

  % handle Dirichlet condition
  [H,R]=gf_asm('dirichlet', 2, mim, mfu, mfd, repmat(eye(3),[1,1,nbd]), zeros(3, nbd));
  [N,U0]=gf_spmat_get(H, 'dirichlet_nullspace', R);
  KK=N'*K*N;
  FF=N'*F;
  % solve ...
  disp('solving...'); t0 = cputime;
  lsolver = 1 % change this to compare the different solvers
  if (lsolver == 1),     % conjugate gradient
    P=gfPrecond('ildlt',KK);
    UU=gf_linsolve('cg',KK,FF,P,'noisy','res',1e-9);
  elseif (lsolver == 2), % superlu
    UU=gf_linsolve('superlu',KK,FF);
  else                   % the matlab "slash" operator
    UU=KK \ FF;
  end;
  disp(sprintf('linear system solved in \%.2f sec', cputime-t0));
  U=(N*UU).'+U0;

In |gfi|, the assembly of vectors, and matrices is done via the ``gf_asm``
function. The Dirichlet condition :math:`u(x) = r(x)` is handled in the weak form
:math:`\int (h(x)u(x)).v(x) = \int r(x).v(x)\quad \forall v` (where :math:`h(x)`
is a :math:`3\times3` matrix field -- here it is constant and equal to the
identity). The reduced system ``KK UU = FF`` is then built via the elimination of
Dirichlet constraints from the original system. Note that it might be more
efficient (and simpler) to deal with Dirichlet condition via a penalization
technique.


Other examples
--------------

* the :file:`demo_refine.m` script shows a simple 2D or 3D bar whose extremity is
  clamped. An adaptative refinement is used to obtain a better approximation in
  the area where the stress is singular (the transition between the clamped area
  and the neumann boundary).

* the :file:`demo_nonlinear_elasticity.m` script shows a 3D bar which is is
  bended and twisted. This is a quasi-static problem as the deformation is
  applied in many steps. At each step, a non-linear (large deformations)
  elasticity problem is solved.

* the :file:`demo_stokes_3D_tank.m` script shows a Stokes (viscous fluid) problem
  in a tank. The :file:`demo_stokes_3D_tank_draw.m` shows how to draw a nice plot
  of the solution, with mesh slices and stream lines. Note that the
  :file:`demo_stokes_3D_tank_alt.m` is the old example, which uses the deprecated
  ``gf_solve`` function.

* the :file:`demo_bilaplacian.m` script is just an adaption of the |gf| example
  :file:`tests/bilaplacian.cc`. Solve the bilaplacian (or a Kirchhoff-Love plate
  model) on a square.

* the :file:`demo_plasticity.m` script is an adaptation of the |gf| example
  :file:`tests/plasticity.cc`: a 2D or 3D bar is bended in many steps, and the
  plasticity of the material is taken into account (plastification occurs when
  the material's Von Mises exceeds a given threshold).

* the :file:`demo_wave2D.m` is a 2D scalar wave equation example (diffraction of
  a plane wave by a cylinder), with high order geometric transformations and high
  order FEMs.


Using Octave/Matlab Object-Oriented features
--------------------------------------------

The basic functions of the |gf| toolbox do not use any advanced |octv| or |mlab| features
(except that the handles to getfem objects are stored in a small structure). But the toolbox comes with a set of objects, which encapsulate
the handles and make them look as real |octv| / |mlab| objects. The aim is not to provide
extra-functionalities, but to have a better integration of the toolbox.

Here is an example of its use::

  >> m=gf_mesh('cartesian',0:.1:1,0:.1:1)
  m =
       id: 0
      cid: 0

  >> m2=gfMesh('cartesian',0:.1:1,0:.1:1)
  gfMesh object ID=1 [17512 bytes], dim=2, nbpts=121, nbcvs=100
  % while \kw{m} is a simple structure, \kw{m2} has been flagged
  % as  an object of class gfMesh.  Since the \texttt{display} method for
  % these  objects  have  been  overloaded,  the  toolbox  displays  some
  % information about the mesh instead of the content of the structure.
  >> gf_mesh_get(m,'nbpts')
  ans =
     121
  % pseudo member access (which calls ##gf_mesh_get(m2,'nbpts'))
  >> m2.nbpts
  ans =
     121

Refer to the OO-commands reference :ref:`mlab-oocmd` for more details.
