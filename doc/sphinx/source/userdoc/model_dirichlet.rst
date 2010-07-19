.. $Id: model.rst 3655 2010-07-17 20:42:08Z renard $

.. include:: ../replaces.txt

.. highlightlang:: c++

.. index:: models, model bricks

.. _ud-model-Dirichlet:


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
should be previously defined in the model. If the data is omitted, the Dirichlet
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

The parameter ``mf_mult`` is replaced by an integer ``degree`` indicating that the
multiplier will be built on a classical finite element method of that degree.

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
                                            dataname = std::string(),
					    *mf_mult = 0);

The penalization consists in computing the mass matrix of the variable and add it
multiplied by the penalization coefficient to the stiffness matrix.
The parameter `mf_mult` (a pointer to a ``getfem::mesh_fem`` object) is optional. It allows to weaken the Dirichlet condition for locking situations. In that case, the penalization matrix is of the form :math:`B^TB` where :math:`B` is the "mass matrix" on the boundary between the shape functions of the variable `varname` and the shape function of the multiplier space.
The penalization coefficient is added as a data of the model and can be
changed thanks to the function::

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


  add_generalized_Dirichlet_condition_with_multipliers(md, mim, varname,
                                           multname, region,
                                           dataname, Hname);


  add_generalized_Dirichlet_condition_with_multipliers(md, mim, varname,
                                           mf_mult, region,
                                           dataname, Hname);

  add_generalized_Dirichlet_condition_with_multipliers(md, mim, varname,
                                           degree, region,
                                           dataname, Hname);


  add_generalized_Dirichlet_condition_with_penalization(md, mim, varname,
                                            penalization_coeff, region,
                                            dataname, Hname);
