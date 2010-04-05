.. $Id$

.. include:: ../replaces.txt

.. highlightlang:: c++

.. index:: models, model bricks

.. _ud-model:

=====================
The model description
=====================

This part is a work in progress for |gf| 4.0. The model description of |gf| allows
to quickly build some fem applications on complex linear or nonlinear PDE coupled
models. The principle is to propose predefined bricks which can be assembled to
describe a complex situation. A brick can describe either an equation (Poisson
equation, linear elasticity ...) or a boundary condition (Dirichlet, Neumann ...)
or any relation between two variables. Once a brick is written, it is possible to
use it in very different situations. This allows a reusability of the produced
code and the possibility of a growing library of bricks. An effort as been made in
order to facilitate as much as possible the definition of a new brick. A brick is
mainly defined by its contribution in the tangent linear system to be solved.

This model description is an evolution of the model bricks of previous versions of
|gf|. Compared to the old system, it is more flexible, more general, allows the
coupling of model (multiphysics) in a easier way and facilitate the writing of new
components. It also facilitate the write of time integration schemes for evolving
PDEs.

The kernel of the model description is contained in the file
:file:`getfem/getfem_models.h`. The two main objects are the |mo| and the |br|.


The model object
----------------

The aim of the |mo| object, defined in file :file:`getfem/getfem_models.h`, is to
globally describe a PDE model. It mainly contains two lists: a list of variables
(related or not to the |mf| objects) and data (also related or not to the |mf|
objects) and a list of bricks. The role of the |mo| object is to coordinate the
module and make them produce a linear system of equations. If the model is
linear, this will simply be the linear system of equation on the corresponding
dofs. If the model is nonlinear, this will be the tangent linear system. There is
two versions of the |mo| object: a real one and complex one.

The declaration of a model object is done by::

  getfem::model md(complex_version = false);

The parameter of the constructor is a boolean which sets if the model deals with
complex number or real numbers. The default is false for a model dealing with real
numbers.

.. _ud-fig-syslin:
.. figure:: images/getfemuserlinearsys.png
   :align: center
   :width: 7cm

   The (tangent) linear system

There are different kinds of variables/data in the model. The variables are the 
unknown of the model. They will be (generally) computed by solving the (tangent) 
linear system build by the model. Generally, the model will have several 
variables. Each variable as a certain size (number of degrees of freedom) and the 
different variables are sorted in alphanumeric order to form the global unknown 
(:math:`U` in Fig. :ref:`ud-fig-syslin`). Each variable will be associated to an 
interval :math:`I = [n_1, n_2]` which will represent the degrees of freedom 
indices correspondig to this variable in the global system. The model stores also 
some data (in the same format than the variables). The difference between data 
and variables is that a data is not an unknown of the model. The value of the 
data should be provided. In some cases (nonlinear models) some variables can be 
considered as some data for certain terms. Variables and data are of two kind. 
They can have a fixed size, or they can depend on a finite element method (be the 
d.o.f. of a finite element method).

For instance, in the situation described in Fig. :ref:`ud-fig-syslin`, there is 
four variables in the model, namely :math:`X, Y, V` and :math:`W`. The role of 
the model object will be to assemble the linear system, i.e. to fill the sub 
matrices corresponding to each variable (:math:`R_{X,X}, R_{Y,Y}, R_{V,V}`, and 
:math:`R_{W,W}`) and the coupling terms between two variables (:math:`R_{X,Y}, 
R_{X,V}, R_{W,V}, \cdots`). This different contributions will be given by the 
different bricks added to the model.

The main usefull methods on a |mo| object are

.. cfunction:: m.is_complex()

   A boolean which says if the model deals with real or complex unknowns and data.

.. cfunction:: add_fixed_size_variable(name, size, niter=1)

   Add a variable of fixed size. ``name`` is a string which designate the
   variable. ``niter`` is the number of copy of the variable (used for time
   integration schemes).

.. cfunction:: add_fixed_size_data(name, size, niter=1)

   Add a data of fixed size. ``name`` is a string which designate the data.
   ``niter`` is the number of copy of the data (used for time integration
   schemes).

.. cfunction:: add_initialized_fixed_size_data(name, V)

   Add a data of fixed size initialized with the given vector ``V``. ``name`` is a
   string which designate the data.

.. cfunction:: add_initialized_scalar_data(name, e)

   Add a data of size 1 initialized with the given scalar value ``e``. ``name`` is
   a string which designate the data.

.. cfunction:: add_fem_variable(name, mf, niter=1)

   Add a variable being the dofs of a finite element method ``mf``. ``name`` is a
   string which designate the variable. ``niter`` is the number of copy of the
   variable (used for time integration schemes).

.. cfunction:: add_fem_data(name, mf, niter=1)

   Add a data being the dofs of a finite element method ``mf``. ``name`` is a
   string which designate the data. ``niter`` is the number of copy of the data
   (used for time integration schemes).

.. cfunction:: add_initialized_fem_data(name, mf, V, niter=1)

   Add a data being the dofs of a finite element method ``mf`` initialized with
   the given vector ``V``. ``name`` is a string which designate the data.
   ``niter`` is the number of copy of the data (used for time integration
   schemes).

.. cfunction:: add_multiplier(name, mf, primal_name, niter=1)

   Add a special variable linked to the finite element method ``mf`` and being a
   multiplier for certain constraints (Dirichlet condition for instance) on a
   primal variable ``primal_name``. The most important is that the degrees of
   freedom will be filtered thanks to a ``partial_mesh_fem`` object in order to
   retain only a set of linearly independent constraints. To ensure this, a call
   to the bricks having a term linking the multiplier and the primal variable is
   done and a special algorithm is called to extract independent constraints. This
   algorithm is optimized for boundary multipliers (see gmm::range_basis). Use it
   with care for volumic multipliers. ``niter`` is the number of copy of the
   variable (used for time integration schemes). Note that for complex terms, only
   the real part is considered to filter the multiplier.

.. cfunction:: real_variable(name, niter=1)

   Gives the access to the vector value of a variable or data. Real version.

.. cfunction:: complex_variable(name, niter=1)

   Gives the access to the vector value of a variable or data. Complex version.

.. cfunction:: mesh_fem_of_variable(name)

   Gives a reference on the |mf| on which the variable is defined. Throw an
   exception if this is not a fem variable.

.. cfunction:: real_tangent_matrix()

   Gives the access to tangent matrix. Real version. A computation of the tangent
   system have to be done first.

.. cfunction:: complex_tangent_matrix()

   Gives the access to tangent matrix. Complex version. A computation of the
   tangent system have to be done first.

.. cfunction:: real_rhs()

   Gives the access to right hand side vector of the linear system. real version.
   A computation of the tangent system have to be done first.

.. cfunction:: complex_rhs()

   Gives the access to right hand side vector of the linear system. Complex
   version. A computation of the tangent system have to be done first.


The |br| object
---------------

A model brick is an object which is supposed to represent a part of a model. It
aims to represent some integral terms in a weak formulation of a pde model. The
model object will contain a list of brick. All the terms described by the brick
will be finally assembled to build the linear system to be solved (the tangent
linear system for a nonlinear problem). For instance if a term :math:`\Delta u` is
present on the pde model (Laplacian of :math:`u`) then the weak formulation will
contain the term :math:`\int_{\Omega}\nabla u\cdot\nabla v\ dx`, where :math:`v`
is the test function corresponding to :math:`u`. Then the role of the
correspponding brick is to assemble the term :math:`\int_{\Omega}\nabla\varphi_i
\cdot\nabla\varphi_j\ dx`, where :math:`\varphi_i` and :math:`\varphi_j` are the
shape functions of the finite element method describing :math:`u`. This term will
be added by the model object to the global linear system on a diagonal block
corresponding to the variable :math:`u`. The only role of the brick is thus to
call the corresponding assembly procedure when the model object ask for it. The
construction of a brick for such a linear term is thus very simple.

Basically, the brick object will derive from the object ``virtual_brick`` defined
in :file:`getfem/getfem_models.h` and should redefine the method
``asm_real_tangent_terms`` or ``asm_complex_tangent_terms`` depending on whether
it is a real term or an intrinsic complex term.


How to build a new brick
------------------------

According to the spirit in which the brick has been designed, a brick should avoid
as much as possible to store additional data. The parameter of a brick should be
contained in the variable and data of the model. For instance, the parameter of a
linear elasticity brick are the elasticity coefficient. This coefficients have to
be some data of the model. When the brick is called by the model obejct, a list of
variables and data is given to the brick. The great majority of the predefined
bricks do not store any data. This allows to instantiate such a bricks only once.

An example of a brick corresponding to the laplacian term is the following (other
examples can be found in the file :file:`getfem_models.cc` which contains the
very standards bricks)::

  struct my_Laplacian_brick: public getfem::virtual_brick {

    void asm_real_tangent_terms(const getfem::model &md, size_type ib,
                                const getfem::model::varnamelist &varl,
                                const getfem::model::varnamelist &datal,
                                const getfem::model::mimlist &mims,
                                getfem::model::real_matlist &matl,
                                getfem::model::real_veclist &vecl,
                                getfem::model::real_veclist &vecl_sym,
                                size_type region, build_version nl) const {
      GMM_ASSERT1(matl.size() == 1,
                  "My Laplacian brick has one and only one term");
      GMM_ASSERT1(mims.size() == 1,
                  "My Laplacian brick need one and only one mesh_im");
      GMM_ASSERT1(varl.size() == 1 && datal.size() == 0,
                  "Wrong number of variables for my Laplacian brick");

      const getfem::mesh_fem &mf_u = md.mesh_fem_of_variable(varl[0]);
      const getfem::mesh_im &mim = *mims[0];

      gmm::clear(matl[0]);
      getfem::asm_stiffness_matrix_for_homogeneous_laplacian
      (matl[0], mim, mf_u, region);
    }

    my_Laplacian_brick(void)
    { set_flags("My Laplacian brick", true /* linear */,
                                      true /* symmetric */,
                                      true /* coercivity */,
                                      true /* real version defined */,
                                      false /* no complex version*/);
    }
  };

The constructor of a brick should call the method ``set_flags``. The first
parameter of this method is a name for the brick (this allows to list the bricks
of a model and facilitate their identification). The other parameters are some
flags, respectively:

* if the brick terms are all linear or not.

* if the brick terms are globally symmetric (conjugated in the complex version) or
  at least do not affect the symmetry. The terms corresponding to two different
  variables and declared symmetric are added twice in the global linear system
  (the term and the transpose of the term).

* if the terms do not affect the coercivity.

* if the terms have a real version or not. If yes, the method
  ``asm_real_tangent_terms`` should be redefined.

* if the terms have a complex version or not. If yes, the method
  ``asm_complex_tangent_terms`` should be redefined.

The method ``asm_real_tangent_terms`` will be called by the model object for the
assembly of the tangent system. The model object gives the whole framework to the
brick to build its terms. The parameter ``md`` of the ``asm_real_tangent_terms``
method is the model that called the brick, ``ib`` being the brick number in the
model. The parameter ``varl`` is an array of variable/data names defined in this
model and needed in the brick. ``mims`` is an array of |mim| pointers. It
corresponds to the integration methods needed to assemble the terms. ``matl`` is
an array of matrices to be computed. ``vecl`` is an array of vectors to be
computed (rhs or residual vectors).  ``vecl_sym`` is an array of vectors to be
computed only for symmetric terms and corresponding to the rhs of the second
variable. A brick can have an arbitrary number of terms. For each term, at least
the corresponding matrix or the corresponding vector has to be filled (or both the
two, but only in the nonlinear case, see the description of the terms below, next
section). ``region`` is a mesh region number indicated that the terms have to be
assembled on a certain region. ``nl`` is for nonlinear bricks only. It says if the
tangent matrix or the residual or both the two are to be computed (for linear
bricks, all is to be computed at each call).

For the very simple Laplacian brick defined above, only one variable is used and
no data and there is only one term. The lines::

      GMM_ASSERT1(matl.size() == 1,
                  "My Laplacian brick has one and only one term");
      GMM_ASSERT1(mims.size() == 1,
                  "My Laplacian brick need one and only one mesh_im");
      GMM_ASSERT1(varl.size() == 1 && datal.size() == 0,
                  "Wrong number of variables for my Laplacian brick");

are not mandatory and just verify that the good number of terms (1), integration
methods (1), variables(1), data(0) are passed to the ``asm_real_tangent_terms``
method.

The lines::

      const getfem::mesh_fem &mf_u = md.mesh_fem_of_variable(varl[0]);
      const getfem::mesh_im &mim = *mims[0];

takes the |mf| object from the variable on which the Laplacian term will be added
and the |mim| object in the list of integrations methods. Finally, the lines::

      gmm::clear(matl[0]);
      getfem::asm_stiffness_matrix_for_homogeneous_laplacian
      (matl[0], mim, mf_u, region);

call a standard assembly procedure for the Laplacian term defined in the file
:file:`getfem/getfem_assembling.h`. The clear method is necessary because
although it is guaranteed that the matrices in ``matl`` have the good sizes they
maybe not cleared before the call of ``asm_real_tangent_terms``.

Note that this simple brick have only one term and is linear. In the case of a
linear birck, either the matrix or the right hand side vector have to be filled
but not both the two. Depending on the declaration of the term. See below the
integration of the birck to the model.

Le us see now a second example of a simple brick which prescribe a Dirichlet
condition thanks to the use of a Lagrange multiplier. The Dirichlet condition is
of the form

.. math::

   u = u_D \text{ on } \Gamma,

where :math:`u` is the variable, :math:`u_D` is a given value and :math:`\Gamma`
is a part on the boundary of the considered domain. The weak terms corresponding
to this condition prescribed with a Lagrange multiplier are

.. math::

   \int_{\Gamma} u \mu\ d\Gamma = \int_{\Gamma} u_D \mu\ d\Gamma, \forall \mu \in M,

where :math:`M` is a appropriate multiplier space. The contributions to the 
global linear system can be viewed in Fig. :ref:`ud-fig-syslinDir`. The matrix 
:math:`B` is the "mass matrix" between the finite element space of the variable 
:math:`u` and the finite element space of the multiplier :math:`\mu`. 
:math:`L_{u}` is the right end side corresponding to the data :math:`u_D`.

.. _ud-fig-syslinDir:
.. figure:: images/getfemuserlinsysDir.png
   :align: center
   :width: 7cm

   Contributions of the simple Dirchlet brick

The brick can be defined as follows::

  struct my_Dirichlet_brick: public getfem::virtual_brick {

    void asm_real_tangent_terms(const getfem::model &md, size_type ib,
                                const getfem::model::varnamelist &varl,
                                const getfem::model::varnamelist &datal,
                                const getfem::model::mimlist &mims,
                                getfem::model::real_matlist &matl,
                                getfem::model::real_veclist &vecl,
                                getfem::model::real_veclist &vecl_sym,
                                size_type region, build_version nl) const {
      GMM_ASSERT1(matl.size() == 1,
                  "My Dirichlet brick has one and only one term");
      GMM_ASSERT1(mims.size() == 1,
                  "My Dirichlet brick need one and only one mesh_im");
      GMM_ASSERT1(varl.size() == 2 && datal.size() == 1,
                  "Wrong number of variables for my Laplacian brick");

      const getfem::mesh_fem &mf_u = md.mesh_fem_of_variable(varl[0]);
      const getfem::mesh_fem &mf_mult = md.mesh_fem_of_variable(varl[1]);
      const getfem::mesh_im &mim = *mims[0];
      const getfem::model_real_plain_vector &A = md.real_variable(datal[ind]);
      const getfem::mesh_fem *mf_data = md.pmesh_fem_of_variable(datal[ind]);

      if (mf_data)
        getfem::asm_source_term(vecl[0], mim, mf_mult, *mf_data, A, region);
      else
        getfem::asm_homogeneous_source_term(vecl[0], mim, mf_mult, A, region);

      gmm::clear(matl[0]);
      getfem::asm_mass_matrix(matl[0], mim, mf_mult, mf_u, region);
    }

    my_Dirichlet_brick(void)
    { set_flags("My Dirichlet brick", true /* linear */,
                                      true /* symmetric */,
                                      false /* coercivity */,
                                      true /* real version defined */,
                                      false /* no complex version */);
    }
  };

This brick has again only one term but define both the matrix and the right hand
side parts. Two variables are concerned, the primal variable on which the
Dirichlet condition is prescribed, and the multiplier variable which should be
defined on a mesh region corresponding to a boundary (it should be added to the
model with the method ``add_multiplier``). The term of the brick will be declared
symmetric (see the next section).

The lines::

      const getfem::model_real_plain_vector &A = md.real_variable(datal[ind]);
      const getfem::mesh_fem *mf_data = md.pmesh_fem_of_variable(datal[ind]);

allow to have the access to the value of the data corresponding to the right hand
side of the Dirichlet condition and to the |mf| on which this data is defined. If
the data is constant (no derscribe on a fem) then ``mf_data`` is a null pointer.

The lines::

      if (mf_data)
        getfem::asm_source_term(vecl[0], mim, mf_mult, *mf_data, A, region);
      else
        getfem::asm_homogeneous_source_term(vecl[0], mim, mf_mult, A, region);

make the assembly of the right hand side. The two versions correspond to a data
defined on a finite element method or constant size data.

( + some example with a nonlinear term ... )


How to add the brick to a model
-------------------------------

In order to add a brick to a model, a certain information have to be passed to the
model:

* A pointer to the brick itself.
* The set of variable names concerned with the terms of the brick.
* The set of data names concerned with the terms of the brick.
* A list of terms description.
* A list of integration methods.
* Eventually the concerned mesh region.

This is done by the call of the |mo| object method::

   md.add_brick(pbr, const getfem::model::varnamelist &varnames,
                     const getfem::model::varnamelist &datanames,
                     const getfem::model::termlist &terms,
                     const getfem::model::mimlist &mims,
                     size_t region);

The method return the index of the brick in the model. The call of this method is
rather complex because it can be adapted to many situations. The construction of a
new brick should be accompagned to the definition of a function that add the new
brick to the model calling this method and more simple to use.

For instance, for the simple Laplacian brick described above, this function can be
defined as folows::

  size_t add_my_Laplacian_brick(getfem::model &md, const getfem::mesh_im &mim,
                                const std::string &varname,
                                size_t region = size_t(-1)) {
    getfem::pbrick pbr = new my_Laplacian_brick;
    getfem::model::termlist tl;

    tl.push_back(getfem::model::term_description(varname, varname, true));
    return md.add_brick(pbr, getfem::model::varnamelist(1, varname),
                        getfem::model::varnamelist(), tl,
                        getfem::model::mimlist(1, &mim), region);
  }

This function will be called by the user of your brick. The type
``getfem::model::varnamelist`` is a ``std::vector<std::string>`` and represent an
array of variable names. The type ``getfem::model::mimlist`` is a
``std::vector<const getfem::mesh_im *>`` and represent an array of pointers to
integration methods. The type ``getfem::model::termlist`` is an array of terms
description. There is two kind of terms. The terms adding only a right hand side
to the linear (tangent) system which have to be added to the list by::

  tl.push_back(getfem::model::term_description(varname));

and the terms having a contribution to the matrix of the linear system which have
to be added to the list by::

  tl.push_back(getfem::model::term_description(varname1, varname2, true/false));

In this case, the matrix term is added in the rows corresponding to the variable
``varname1`` and the columns corresponding to the variable ``varname2``. The
boolean being the third parameter is to declare if the term is symmetric or not.
If it is symmetric and if the two variables are different then the assembly
procedure add the corresponding term AND its transpose. The number of terms is
arbitrary. For each term declared, the brick have to fill the corresponding right
hand side vector (parameter ``vecl`` of ``asm_real_tangent_terms`` above) or/and
the matrix term (parameter ``matl`` of ``asm_real_tangent_terms``) depending on
the declaration of the term. Note that for nonlinear bricks, both the matrix and
the right hand side vectors have to be filled. For linear bricks, if the right
hand side is filled for a term declared to be a matrix term, it is IGNORED.

The variable names and the data names are given in two separate arrays because the
dependence of the brick is not the same in both cases. A linear term have to be
recomputed if the value of a data is changed but not if the value of a variable is
changed.

The function allowing to add the simple Dirichlet brick described above can be
defined as follows::

  size_t add_my_Dirichlet_condition_brick(model &md, const mesh_im &mim,
                                          const std::string &varname,
                                          const std::string &multname,
                                          size_t region,
                                          const std::string &dataname) {
    pbrick pbr = new my_Dirichlet_brick;
    model::termlist tl;
    tl.push_back(model::term_description(multname, varname, true));
    model::varnamelist vl(1, varname);
    vl.push_back(multname);
    model::varnamelist dl;
    if (dataname.size()) dl.push_back(dataname);
    return md.add_brick(pbr, vl, dl, tl, model::mimlist(1, &mim), region);
  }

Again, here, the term is declared symmetric and then the matrix term and its
transpose will be added.


Generic elliptic brick
----------------------

This brick add an elliptic term on a variable of a model.  The shape of the
elliptic term depends both on the variable and a given coefficient. This
corresponds to a term:

.. math::

   -\text{div}(a\nabla u),

where :math:`a` is the coefficient and :math:`u` the variable. The coefficient can
be a scalar, a matrix or an order four tensor. The variable can be vector valued
or not. This means that the brick treats several different situations. If the
coefficient is a scalar or a matrix and the variable is vector valued then the
term is added componentwise. An order four tensor coefficient is allowed for
vector valued variable only.  The coefficient can be constant or described on a
FEM. Of course, when the coefficient is a tensor described on a finite element
method (a tensor field) the corresponding data can be a huge vector. The
components of the matrix/tensor have to be stored with the fortran order
(columnwise) in the data vector corresponding to the coefficient (compatibility
with blas). The symmetry and coercivity of the given matrix/tensor is not verified
(but assumed).

This brick can be added to a model ``md`` thanks to two functions. The first one
is::

  size_type getfem::add_Laplacian_brick(md, mim, varname, region = -1);

that adds an elliptic term relatively to the variable ``varname`` of the model
with a constant coefficient equal to :math:`1` (a Laplacian term). This
corresponds to the Laplace operator. ``mim`` is the integration method which will
be used to compute the term. ``region`` is an optional region number. If it is
ommited, it is assumed that the term will be computed on the whole mesh. The
result of the function is the brick index in the model.

The second function is::

  size_type getfem::add_generic_elliptic_brick(md, mim, varname, dataname, region = -1);

It adds a term with an arbitrary coefficient given by the data ``dataname`` of the
model. This data have to be defined first in the model.

Note that very general equations can be obtained with this brick. For instance,
linear anisotropic elasticity can be obtained with a tensor data. When an order
four tensor is used, the corresponding weak term is the following

.. math::

   \int_{\Omega} \sum_{i,j,k,l} a_{i,j,k,l}\partial_i u_j \partial_k v_l dx

where :math:`a_{i,j,k,l}` is the order four tensor and :math:`\partial_i u_j` is
the partial derivative with respect to the :math:`i^{th}` variable of the
component :math:`j` of the unknown :math:`k`. :math:`v` is the test function.
However, for linear isotropic elasticity, a more adapted brick is available (see
below).

The brick has a working complex version.


Dirichlet condition brick
-------------------------

The aim of the Dirichlet condition brick is to prescribe a Dirichlet condition on
a part of the boundary of the domain for a variable of the model. This means that
the value of this variable is prescribed on the boundary. There is two versions of
this brick. The first version prescribe the Dirichlet thank to a multiplier. The
associated weak form of the term is the following:

.. math::

   \int_{\Gamma} u \mu d\Gamma = \int_{\Gamma} u_D \mu d\Gamma, \forall \mu \in M.

where :math:`u` is the variable, :math:`M` is the space of multipliers, :math:`u`
is the variable and :math:`\Gamma` the Dirichlet boundary. For this version, an
additional variable have to be added to represent the multiplier. It can be done
directly to the model or thanks to the functions below. There are three functions
allowing to add a Dirichlet condition prescribed with a multiplier. The first one
is::

  add_Dirichlet_condition_with_multipliers(md, mim, varname,
                                           multname, region,
                                           dataname = std::string());

adding a Dirichlet condition on ``varname`` thanks to a multiplier variable
``multname`` on the mesh region ``region`` (which should be a boundary). The value
of the variable on that boundary is described by the data ``dataname`` which
should be previously defined in the model. If the data is ommitted, the Dirichlet
condition is assumed to be an homogeneous one (vanishing variable on the
boundary). The data can be constant or described on a FEM. It can also be scalar
or vector valued, depending on the variable. The variable ``multname`` should be
added to the model by the method ``add_multiplier``. The function returns the
brick index in the model. The second function is::

  add_Dirichlet_condition_with_multipliers(md, mim, varname,
                                           mf_mult, region,
                                           dataname = std::string());

The only difference is that ``multname`` is replaced by ``mf_mult`` which means
that only the finite element on which the multiplier will be built is given. The
function adds itself the multiplier variable to the model. The third function is
very similar::

  add_Dirichlet_condition_with_multipliers(md, mim, varname,
                                           degree, region,
                                           dataname = std::string());

The parameter ``mf_mult`` is replaced by a integer ``degree`` indicating that the
multiplier will be build on a classical finite element method of that degree.

Note, that in all the cases, when a variable is added by the method
``add_multiplier`` of the model object, the |mf| will be filtered (thank to a
``partial_mesh_fem_object`` in order to retain only the degrees of freedom having
a non vanishing contribution on the considered boundary.

Finally, the variable name of the multiplier can be obtained thank to the
function::

  mult_varname_Dirichlet(md, ind_brick);

where ``ind_brick`` is the brick index in the model. This function has an
undefined behavior if it applied to another kind of brick.

The second version of the Dirichlet condition brick is the one with penalization.
The function allowing to add this brick is::

  add_Dirichlet_condition_with_penalization(md, mim, varname,
                                            penalization_coeff, region,
                                            dataname = std::string());

The penalization consists in computing the mass matrix of the variable and add it
multiplied by the penalization coefficient to the stiffness matrix. The
penalization coefficient is added as a data of the model and can be changed thanks
to the function::

  change_penalization_coeff(md, ind_brick, penalisation_coeff);


Generalized Dirichlet condition brick
-------------------------------------

The generalized Dirichlet condition is a boundary condition of a vector field u of 
the type

.. math::

   H u  = r

where :math:`H` is a matrix field. The functions adding the corresponding bricks 
are similar to the ones of the standard Dirichlet condition except that they need 
the supplementary parameter `Hname` which gives the name of the data corresponding 
to :math:`H`. This data can be a matrix field described on a scalar fem or a 
constant matrix.::


  add_Dirichlet_condition_with_multipliers(md, mim, varname,
                                           multname, region,
                                           dataname, Hname);


  add_Dirichlet_condition_with_multipliers(md, mim, varname,
                                           mf_mult, region,
                                           dataname, Hname);

  add_Dirichlet_condition_with_multipliers(md, mim, varname,
                                           degree, region,
                                           dataname, Hname);


  add_Dirichlet_condition_with_penalization(md, mim, varname,
                                            penalization_coeff, region,
                                            dataname, Hname);


Source term bricks (and Neumann condition)
------------------------------------------

This brick add a source term, i.e. a term which occurs only in the right hand side
of the linear (tangent) system build by the model. If :math:`f` denotes the value
of the source term, the weak form of such a term is

.. math::

   \int_{\Omega} f v\ dx

where :math:`v` is the test function. The value :math:`f` can be constant or
described on a finite element method.

It can also represent a Neumann condition if it is applied on a boundary of the
domain.

The function to add a source term to a model is::

  add_source_term_brick(md, mim,
                        varname, dataname, region = -1,
                        directdataname = std::string());

where ``md``is the model object, ``mim`` is the integration method, ``varname`` is
the variable of the model for which the source term is added, ``dataname`` is the
name of the data in the model which represents the source term. It has to be
scalar or vector valued depending on the fact that the variable is scalar or
vector valued itself. ``region`` is a mesh region on which the term is added. If
the region corresponds to a boundary, the source term will represent a Neumann
condition. ``directdataname`` is an optional additional data which will directly
be added to the right hand side without assembly.

The brick has a working complex version.

A slightly different brick, especially dedicated to deal with a Neumann condition,
is added by the following function::

  add_normal_source_term_brick(md, mim,
                               varname, dataname, region);

The difference compared to the basic source term brick is that the data should be
a vector field (a matrix field if the variable ``varname`` is itself vector
valued) and a scalar product with the outward unit normal is performed on it.


Predefined solvers
------------------

Of course, for many problems, it will be more convenient to make a specific
solver. Even so, one generic solver is available to test your models quickly. It
can also be taken as an example to build your own solvers. It is defined in
:file:`getfem/getfem_model_solvers.h` and the call is::

  getfem::standard_solve(md, iter);

where ``md`` is the model object and ``iter`` is an iteration object from |gmm|.
See also the next section for an example of use.

Note that |sLU| is used by default on "small" problems. You can also link
|mumps| with |gf| (see section :ref:`ud-linalg`) and used the parallele
version.


Example of a complete Poisson problem
-------------------------------------

The following example is a part of the test program
:file:`tests/laplacian_with_bricks.cc`. Construction of the mesh and finite
element methods are omitted. It is assumed that a mesh is build and two finite
element methods ``mf_u`` and ``mf_rhs`` are build on this mesh. Is is also
assumed that ``NEUMANN_BOUNDARY_NUM`` and ``DIRICHLET_BOUNDARY_NUM`` are two
valid boundary indices on that mesh. The code begins by the definition of three
functions which are interpolated on ``mf_rhs`` in order to build the data for the
source term, the Neumann condition and the Dirichlet condition. Follows the
declaration of the model object, the addition of the bricks and the solving of
the problem::

  using bgeot::base_small_vector;
  // Exact solution. Allows an interpolation for the Dirichlet condition.
  scalar_type sol_u(const base_node &x) { return sin(x[0]+x[1]); }
  // Righ hand side. Allows an interpolation for the source term.
  scalar_type sol_f(const base_node &x) { return 2*sin(x[0]+x[1]); }
  // Gradient of the solution. Allows an interpolation for the Neumann term.
  base_small_vector sol_grad(const base_node &x)
  { return base_small_vector(cos(x[0]+x[1]), cos(x[0]+x[1]); }

  int main(void) {

    // ... definition of a mesh
    // ... definition of a finite element method mf_u
    // ... definition of a finite element method mf_rhs
    // ... definition of a integration method mim
    // ... definition of boundaries NEUMANN_BOUNDARY_NUM
    //                        and DIRICHLET_BOUNDARY_NUM

    // Model object
    getfem::model laplacian_model;

    // Main unknown of the problem
    laplacian_model.add_fem_variable("u", mf_u);

    // Laplacian term on u.
    getfem::add_Laplacian_brick(laplacian_model, mim, "u");

    // Volumic source term.
    std::vector<scalar_type> F(mf_rhs.nb_dof());
    getfem::interpolation_function(mf_rhs, F, sol_f);
    laplacian_model.add_initialized_fem_data("VolumicData", mf_rhs, F);
    getfem::add_source_term_brick(laplacian_model, mim, "u", "VolumicData");

    // Neumann condition.
    gmm::resize(F, mf_rhs.nb_dof()*N);
    getfem::interpolation_function(mf_rhs, F, sol_grad);
    laplacian_model.add_initialized_fem_data("NeumannData", mf_rhs, F);
    getfem::add_normal_source_term_brick
    (laplacian_model, mim, "u", "NeumannData", NEUMANN_BOUNDARY_NUM);

    // Dirichlet condition.
    gmm::resize(F, mf_rhs.nb_dof());
    getfem::interpolation_function(mf_rhs, F, sol_u);
    laplacian_model.add_initialized_fem_data("DirichletData", mf_rhs, F);
    getfem::add_Dirichlet_condition_with_multipliers
    (laplacian_model, mim, "u", mf_u, DIRICHLET_BOUNDARY_NUM, "DirichletData");

    gmm::iteration iter(residual, 1, 40000);
    getfem::standard_solve(laplacian_model, iter);

    std::vector<scalar_type> U(mf_u.nb_dof());
    gmm::copy(laplacian_model.real_variable("u"), U);

    // ... doing something with the solution ...

    return 0;
  }

Note that the brick can be added in an arbitrary order.


Constraint brick
----------------

The constraint brick allows to add an explicit constraint on a variable. Explicit
means that no integration is done. if :math:`U` is a variable then a constraint of
the type

.. math::

   BU = L,

can be added with the two following functions::

  indbrick = getfem::add_constraint_with_penalization(md, varname,
                                                      penalisation_coeff, B, L);
  indbrick = getfem::add_constraint_with_multipliers(md, varname,
                                                     multname, B, L);

In the second case, a (fixed size) variable which will serve as a multiplier
should be first added to the model.

For the penalized version ``B`` should not contain a plain row, otherwise the
whole tangent matrix will be plain. The penalization parameter can be changed
thanks to the function::

  change_penalization_coeff(md, ind_brick, penalisation_coeff);

It is possible to change the constraints at any time thanks to the two following
functions::

  getfem::set_private_data_matrix(md, indbrick, B)
  getfem::set_private_data_rhs(md, indbrick, L)

where ``indbrick`` is the index of the brick in the model.


Other "explicit" bricks
-----------------------

Two (very simple) bricks allow to add some explicit terms to the tangent system.

The function::

  indbrick = getfem::add_explicit_matrix(md, varname1, varname2, B
                                         issymmetric = false,
                                         iscoercive = false);

adds a brick which just adds the matrix ``B`` to the tangent system relatively to
the variables ``varname1`` and ``varname2``. The given matrix should have as many
rows as the dimension of ``varname1`` and as many columns as the dimension of
``varname2``. If the two variables are different and if ``issymmetric`` is set to
true then the transpose of the matrix is also added to the tangent system (default
is false). Set ``iscoercive`` to true if the term does not affect the coercivity
of the tangent system (default is false). The matrix can be changed by the
command::

  getfem::set_private_data_matrix(md, indbrick, B);

The function::

  getfem::add_explicit_rhs(md, varname, L);

add a brick which just add the vector ``L`` to the right hand side of the tangent
system relatively to the variable ``varname``. The given vector should have the
same size as the variable ``varname``. The value of the vector can by changed by
the command::

  getfem::set_private_data_rhs(md, indbrick, L);


Helmholtz brick
---------------

This brick represents the complex or real Helmholtz problem:

.. math::

   \Delta u + k^2 u = \ldots

where :math:`k` the wave number is a real or complex value. For a complex
version, a complex model has to be used (see :file:`tests/helmholtz.cc`).

The function adding a Helmholtz brick to a model is::

  getfem::add_Helmholtz_brick(md, mim, varname, dataname, region);

where ``varname`` is the variable on which the Helmholtz term is added and
``dataname`` should contain the wave number.


Fourier-Robin brick
-------------------

This brick can be used to add boundary conditions of Fourier-Robin type like:

.. math::

   \frac{\partial u}{\partial \nu} = Qu

for scalar problems, or

.. math::

   \sigma\cdot \nu = Qu

for linearized elasticity problems. ``Q`` is a scalar field in the scalar case or
a matrix field in the vectorial case. This brick works for both real or complex
terms in scalar or vectorial problems.

The function adding this brick to a model is::

  add_Fourier_Robin_brick(md, mim, varname, dataname, region);

where ``dataname`` is the data of the model which represents the coefficient
:math:`Q`.

Note that an additional right hand side can be added with a source term brick.


Isotropic linearized elasticity brick
-------------------------------------

This brick represents a term

.. math::

   -div(\sigma) = \ldots

with

.. math::

   \sigma &= \lambda\mbox{tr}(\varepsilon(u))I + 2\mu\varepsilon(u) \\
   \varepsilon(u) &= (\nabla u + \nabla u^T)/2

:math:`\varepsilon(u)` is the small strain tensor, :math:`\sigma` is the stress
tensor, :math:`\lambda` and :math:`\mu` are the Lamé coefficients. This represents
the system of linearized isotropic elasticity. It can also be used with
:math:`\lambda=0` together with the linear incompressible brick to build the
Stokes problem.

The function which adds this brick to a model is::

  ind_brick = getfem::add_isotropic_linearized_elasticity_brick
              (md, mim, varname, dataname_lambda, dataname_mu,
               region = size_type(-1));

where ``dataname_lambda`` and ``dataname_mu`` are the data of the model
representing the Lamé coefficients (constant or described on a finite element
method).

The function::

  getfem::compute_isotropic_linearized_Von_Mises_or_Tresca
    (md, varname, dataname_lambda, dataname_mu, mf_vm, VM, tresca_flag = false);

compute the Von Mises criterion (or Tresca if ``tresca_flag`` is set to true) on
the displacement field stored in ``varname``. The stress is evaluated on the |mf|
``mf_vm`` and stored in the vector ``VM``.

The program :file:`tests/elastostatic.cc` can be taken as a model of use of this
brick.


linear incompressibility (or nearly incompressibility) brick
------------------------------------------------------------

This brick adds a linear incompressibility condition (or a nearly incompressible
condition) in a problem of type:

.. math::

   \mbox{div}(u) = 0,\quad (\mbox{ or } \mbox{div}(u) = \varepsilon p)

This constraint is enforced with Lagrange multipliers representing the pressure,
introduced in a mixed formulation.

The function adding this incompressibility condition is::

  ind_brick = getfem::add_linear_incompressibility
              (md, mim, varname, multname_pressure, region = size_type(-1),
               dataname_penal_coeff = std::string());

where ``varname`` is the variable on which the incompressibility condition is
prescribed, ``multname_pressure`` is a variable which should be described on a
scalar fem representing the multiplier (the pressure) and ``dataname_penal_coeff``
is an optional penalization coefficient (constant or described on a finite element
method) for the nearly incompressible condition.

In nearly incompressible homogeneous linearized elasticity, one has
:math:`\varepsilon = 1 / \lambda` where :math:`\lambda` is one of the Lamé
coefficient and :math:`\varepsilon` the penalization coefficient.

For instance, the following program defines a Stokes problem with a source term
and an homogeneous Dirichlet condition on boundary 0. ``mf_u``, ``mf_data`` and
``mf_p`` have to be valid finite element description on the same mesh. ``mim``
should be a valid integration method on the same mesh::

  typedef std::vector<getfem::scalar_type> plain_vector;
  size_type N = mf_u.linked_mesh().dim();

  getfem::model Stokes_model;

  laplacian_model.add_fem_variable("u", mf_u);

  getfem::scalar_type mu = 1.0;
  Stokes_model.add_initialized_data("lambda", plain_vector(1, 0.0));
  Stokes_model.add_initialized_data("mu", plain_vector(1, mu));

  getfem::add_isotropic_linearized_elasticity_brick(Stokes_model, mim,
                                                    "u", "lambda", "mu");

  laplacian_model.add_fem_variable("p", mf_p);
  getfem::add_linear_incompressibility(Stokes_model, mim, "u", "p");

  plain_vector F(mf_data.nb_dof()*N);
  for (int i = 0; i < mf_data.nb_dof()*N; ++i) F(i) = ...;
  Stokes_model.add_initialized_fem_data("VolumicData", mf_data, F);
  getfem::add_source_term_brick(Stokes_model, mim, "u", "VolumicData");

  getfem::add_Dirichlet_condition_with_multipliers(Stokes_model, mim,
                                                   "u", mf_u, 1);

  gmm::iteration iter(residual, 1, 40000);
  getfem::standard_solve(Stokes_model, iter);

  plain_vector U(mf_u.nb_dof());
  gmm::copy(Stokes_model.real_variable("u"), U);

An example for a nearly incompressibility condition can be found in the program
:file:`tests/elastostatic.cc`.


Mass brick
----------

This brick represents a weak term of the form

.. math::

   \int_{\Omega} \rho u\cdot v\ dx + \ldots

It mainly represents a mass term for transient problems but can also be used for
other applications (it can be used on a boundary). Basically, this brick adds a
mass matrix on the tangent linear system with respect to a certain variable.

The function which adds this brick to a model is::

  ind_brick = getfem::add_mass_brick
              (md, mim, varname, dataname_rho="", region = size_type(-1));

where ``dataname_rho`` is an optional data of the model representing the density
:math:`\rho`. If it is omitted, the density is assumed to be equal to one.

Note that for time integrations scheme, there exist specific bricks for the
discretisation of time derivatives.


The time dispatchers: integration of transient problems
-------------------------------------------------------

The role of time dispatchers is to allow the integration of transient problems 
with some pre-defined time integration schemes. The principle of the time 
dispatchers is to dispatch the terms of a brick on the different time steps of the 
considered time integration scheme. When time derivative terms are present in the 
model (this is generally the case except for quasistatic models), the time 
dispatcher will be associated to a specific brick representing this time 
derivative term (:math:`\partial u / \partial t` or :math:`\partial^2 u / \partial 
t^2` for instance). For this, a number of tools are available in |gf| to help the 
construction of a time dispatcher. Mainly they are the two following:

* The variables can be duplicated to take into account the differents version 
  corresponding to each time iteration. For instance, for simplest time 
  integration scheme, two versions :math:`U^n` and :math:`U^{n+1}` of a variable 
  :math:`U` are stored. The addition of a variable :math:`u` with two versions can 
  be done with the method of the model object::

    model.add_fem_variable("u", mf_u, 2);

  where :math:`2` is here the number of versions. The variable which is actually 
  computed have always the index 0 and will be accessed with 
  ``model.real_variable("u", 0)`` or simply with ``model.real_variable("u")``. It 
  will generally represent the version :math:`U^{n+1}`. The version :math:`U^{n}` 
  (corresponding to the previous time step) will be accessed with 
  ``model.real_variable("u", 1)``. Generally, it will be necessary to set this 
  version with ``model.set_real_variable("u", 1)`` to define the initial condition 
  of the model. At the end of each iteration, the different versions of a variable 
  are automatically shifted (version 0 :math:`\rightarrow` version 1 ...).

* The right hand side of a brick is dispatched into several right hand sides for 
  each time iteration which are stored. To avoid unnecessary computation, the time 
  dispatcher can shift these extra right hand sides at the end of each time 
  iteration.


Theta-method dispatcher
-----------------------

This is the simplest time dispatcher. The use of this dispatcher will be described 
in details. Since the use of the other dispatchers is similar, only their 
specificities will be described later on.

The principle of the :math:`\theta`-method is to dispatch the term :math:`F` into 
:math:`(\theta) F^{n+1} + (1-\theta) F^{n},`

For specific values of :math:`\theta` one obtains some classical schemes: backward 
Euler for :math:`\theta = 1`, forward Euler for :math:`\theta = 0` and 
Crank-Nicholson scheme for :math:`\theta = 1/2` (which is an order two scheme).

For instance, if the dispatcher is applied to a brick representing a linear 
elliptic term :math:`KU` where :math:`K` is the stiffness matrix and :math:`U` the 
unknown, it will be transformed into :math:`(\theta) KU^{n+1} + (1-\theta) 
KU^{n}`.

Since :math:`U^{n+1}` is the real unknown, the effect will be to multiply by 
:math:`\theta` the stiffness matrix and to add to the right hand side the term 
:math:`(1-\theta) KU^{n}`. This means also that :math:`U^{n}` have to be 
initialized (with something like ``gmm::copy(U0, model.real_variable("u",1))``). 
It represents an initial data for the problem. Remember this principle: each time 
you apply a time dispatcher to a brick, the corresponding variables have to have 
the right number of versions (see above) and should be initialized before the 
first time iteration.

You can apply the dispatcher to a brick having only a right hand side (a source 
term for instance). It is not necessary if the term is constant in time.

When a brick represents a constraint (Dirichlet condition, incompressibility ...) 
this is not mandatory to apply the dispatcher. Of course, the result will not be 
exactly the same if you apply or not the dispatcher. If you do not apply it, the 
constraint will be applied to the current variable (:math:`U^{n+1}` for the 
:math:`\theta`-method). If you apply it, the constraint will be in a sense applied 
to :math:`(\theta) U^{n+1} + (1-\theta) U^{n}`. If the constraint is applied 
thanks to a multiplier, this multiplier will need to have different versions and 
will need to have an initial value.

In order to apply the :math:`\theta`-method dispatcher to a set of brick you must 
execute::

  model.add_initialized_scalar_data("theta", theta);
  getfem::add_theta_method_dispatcher(model, transient_bricks, "theta");

where ``transient_bricks`` is a ``dal::bit_vector`` containing the indices of the 
corresponding bricks. The value of :math:`\theta` can be modified from an 
iteration to another.

The global structure of the loop solving the different time steps should be the 
following::

  gmm::iteration solver_iter(residual, 0, 40000);

  // Set here the initial values.

  model.first_iter(); // initialize the iterations.

  for (scalar_type t = 0; t < T; t += dt) {

    solver_iter.init();
    getfem::standard_solve(model, solver_iter); // solve an iteration.

    model.next_iter(); // shift the variables and additional right hand sides.
  }

where ``model.first_iter()`` should be called before the first iteration to 
initialize the right hand side of the time dispatchers. The initial data should be 
set before the call to ``model.first_iter()``. The method ``model.next_iter()`` is 
to be called at the end of each iteration. It calls the dispatcher to shift there 
additional right hand side and it shifts the version of the variables.


Basic first order time derivative brick
+++++++++++++++++++++++++++++++++++++++

A term like :math:`\rho \partial u / \partial t` will be represented in the model 
by :math:`(MU^{n+1} - MU^{n}) / dt`, where :math:`M` is the mass matrix and 
:math:`dt` is the time step. The :math:`\theta`-method is compatible with this. A 
brick is dedicated to represent this term. It can be added to the model by the 
function::

  getfem::add_basic_d_on_dt_brick(model, mim, varname, dataname_dt,
                                  dataname_rho = std::string(),
                                  region = size_type(-1));

where ``varname`` is the name of the variable on which the time derivative is 
applied (should have at least two versions), ``dataname_dt`` is the name of the 
data corresponding to the time step (added by 
``model.add_initialized_scalar_data("dt", dt)`` for instance) which could be 
modified from an iteration to another and ``dataname_rho`` is an optional 
parameter (whose default value is 1) corresponding to the term :math:`\rho` in 
:math:`\rho \partial u / \partial t`.

NOTE that the time dispatcher should not be applied to this brick !

A good model of the use of this brick and the :math:`\theta`-method time 
dispatcher can be found in the test program ``tests/heat_equation.cc``.


Basic second order time derivative brick
++++++++++++++++++++++++++++++++++++++++

This brick represents a second order time derivative like :math:`\rho \partial^2 u 
/ \partial t^2`. The problem with such a term is that the :math:`\theta`-method 
should be applied both on :math:`u` and :math:`\partial u / \partial t` which 
means that :math:`\partial u / \partial t` is a natural unknown of the problem. 
The easiest way is then to add the time derivative of the variable :math:`u` has a 
independent variables of the model (a drawback, of course, is that one has twice 
as much unknowns). This Basic second order time derivative brick does not apply 
this strategy. The time derivative :math:`\partial u / \partial t` is considered 
as a data which is updated at a post-traitment stage (in some cases, this strategy 
cannot be applied if the time derivative appears to be a required unknown of the 
model).

The term :math:`\rho \partial^2 u / \partial t^2` will be represented by 
:math:`(MU^{n+1} - MU^{n}) / (\alpha dt^2) - M V^n / (\alpha dt) ~~~~~~~~(*)`, 
where :math:`M` is the mass matrix, :math:`dt` is the time step, :math:`\alpha` is 
a parameter which is equal to :math:`\theta` for the :math:`\theta`-method and 
:math:`V^n` the time derivative at the previous time step. This means in 
particular that :math:`V` should be added as a data on the model with (at least) 
two versions.

The function adding the brick is::

  getfem::add_basic_d2_on_dt2_brick(model, mim, varname, dataname_V,
             dataname_dt, dataname_alpha, dataname_rho = std::string(),
             region = size_type(-1));

where ``varname`` is the name of the variable on which the second order time 
derivative is applied, ``dataname_V`` is the data representing the time 
derivative, ``dataname_dt`` is the name of the data corresponding to the time step 
(added by ``model.add_initialized_scalar_data("dt", dt)`` for instance) which 
could be modified from an iteration to another, ``dataname_alpha`` is the name of 
the data containing the parameter :math:`\alpha` in (*) and ``dataname_rho`` is an 
optional parameter (whose default value is 1) corresponding to the term 
:math:`\rho` in :math:`\rho \partial^2 u / \partial t^2`.

At the end of each iteration, the data ``dataname_V`` should be updated (before 
the call to ``model.next_iter()`` by the call to::

  getfem::velocity_update_for_order_two_theta_method
      (model, varname, dataname_V, dataname_dt, dataname_alpha);

A good model of the use of this brick and the :math:`\theta`-method time 
dispatcher can be found in the test program ``tests/wave_equation.cc``.


Midpoint dispatcher
-------------------

The principle of the midpoint scheme is to dispacth a term :math:`F(U)` into 
:math:`F((U^{n+1}-U^{n})/2),`

It is different from the Crank-Nicholson scheme (:math:`\theta`-method for 
:math:`\theta=1/2`) only for nonlinear terms.

The real unknown remains :math:`U^{n+1}`. the effect will be to multiply by 
:math:`1/2` the stiffness (or tangent) matrix and to add to a right hand side the 
term :math:`(KU^{n}/2` for a linear matrix term :math:`K`. As for the 
:math:`\theta`-method, the variables have to have two version and the second 
version have to be initialized.

You can apply the dispatcher to a brick having only a right hand side (a source 
term for instance). It is not necessary if the term is constant in time.

NOTE that if the brick depend on a data which is not constant in time, the data 
either have to have to versions (and the mean of the two versions are taken into 
account) or evaluated at the middle of the time step.

When a brick represents a constraint (Dirichlet condition, incompressibility ...) 
this is not mandatory to apply the dispatcher. Of course, the result will not be 
exactly the same if you apply or not the dispatcher. If you do not apply it, the 
constraint will be applied to the current variable :math:`U^{n+1}`. If you apply 
it, the constraint will be applied to :math:`(U^{n+1} + U^{n})/2`. If the 
constraint is applied thanks to a multiplier, this multiplier will need to have 
different versions and will need to have an initial value.

In order to apply the midpoint dispatcher to a set of brick you must execute::

  getfem::add_midpoint_dispatcher(model, transient_bricks);

where ``transient_bricks`` is a ``dal::bit_vector`` containing the indices of the 
corresponding bricks.


Basic first order time derivative brick
+++++++++++++++++++++++++++++++++++++++

The same brick as for the :math:`\theta`-method can be used to represent a first 
order time derivative.


Basic second order time derivative brick
++++++++++++++++++++++++++++++++++++++++

The same brick as for the :math:`\theta`-method can be used to represent a second 
order time derivative. The value of :math:`\alpha` should be :math:`1/2`.


Newmark scheme
--------------

For a system

.. math::

   M\ddot{U} + K(U) = F,

the Newmark scheme of parameter :math:`\beta` and :math:`\gamma` is defined by

.. math::

   M(U^{n+1} - U^{n}) = dt M V^n + dt^2/2( 2\beta(F^{n+1}-K(U^{n+1})) + (1-2\beta)(F^{n}-K(U^{n}))),\\
   M(V^{n+1} - V^{n}) = dt ( 2\gamma(F^{n+1}-K(U^{n+1})) + (1-2\gamma)(F^{n}-K(U^{n}))),

where :math:`V` represents the time derivative of :math:`U`.

The implementation of the Newmark scheme proposed is not optimal and should be 
adapted. It can be optained using the basic second order time derivative brick 
(see :math:`\theta`-method) and the :math:`\theta`-method time dispatcher used 
with :math:`\theta = 2\beta`. Additionaly, one has to use the following function 
which compute the time derivative of the variable as a post-computation::

  getfem::velocity_update_for_Newmark_scheme
      (model, id2dt2, varname, dataname_V, dataname_dt, dataname_alpha);

where ``id2dt2`` is the index of the basic second order time derivative brick (see 
the section on the :math:`\theta`-method for more details and the implementation 
in the test program ``tests/wave_equation.cc``).

This implementation of the Newmark-scheme is not optimal since the latter function 
inverts the mass matrix to compute the time derivative using a conjugate gradient. 
This linear system solve could be avoided by keeping the multiplication of the 
mass matrix with the time derivative as a data, with an addaptation of the time 
derivative brick.


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

  \lambda_T^i = {\mathscr F} \lambda_N^i \frac{\dot{u}_T}{\|\dot{u}_T\|} ~~~ \text{ si } \dot{u}_T \ne 0,

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

The direct nodal condition may have some drawback : locking phenomena, overconstraint. It is in fact often more stable and for the same accuracy to use multiplier of reduced order compared to the displacement (the direct nodal contact condition correspond more or less to a multiplier described on the same finite element method than the displacement).

Let :math:`\varphi_i` be the shapes functions of the finite element describing the displacement and :math:`\psi_i` be the shape functions of a finite element describing a multiplier on the contact boundary :math:`\Gamma_c`. It is assumed that the set of admissible multiplier describing the normal stress will be

.. math::

  \Lambda_N^h = \{ \mu^h_N = \sum \mu^j_N \psi_j : \mu^h_N(a_i) \le 0, ~i = 1..N_c \}

where :math:`a_i`, :math:`~~i=1..N_c` are the finite element nodes corresponding to the multiplier. The discrete contact condition is now expressed in a weak form by

.. math::

  \int_{\Gamma_c} (\mu_N^h - \lambda_N^h) u_N d\Gamma \le 0 ~~ \forall \mu_N^h \in \Lambda_N^h. 

In that case, the component :math:`\lambda_N^i` is a contact stress (:math:`N/m^2`) and the matrix :math:`B_N` can be written

.. math::

  (B_N)_{ij} = \int_{\Gamma_c} \psi_i \varphi_j d\Gamma.

The matrix :math:`B_T` can also be written in a similar way. The friction condition can be written in a weak form

.. math::

  \int_{\Gamma_c} (\mu_T^h - \lambda_T^h) \dot{u}_T d\Gamma \le 0 ~~ \forall \mu_T^h \in \Lambda_T^h({\mathscr F}\lambda_N^h),

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

to be done.
