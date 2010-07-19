.. $Id: model.rst 3655 2010-07-17 20:42:08Z renard $

.. include:: ../replaces.txt

.. highlightlang:: c++

.. index:: models, model bricks

.. _ud-model-poisson:



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
  // Right hand side. Allows an interpolation for the source term.
  scalar_type sol_f(const base_node &x) { return 2*sin(x[0]+x[1]); }
  // Gradient of the solution. Allows an interpolation for the Neumann term.
  base_small_vector sol_grad(const base_node &x)
  { return base_small_vector(cos(x[0]+x[1]), cos(x[0]+x[1]); }

  int main(void) {

    // ... definition of a mesh
    // ... definition of a finite element method mf_u
    // ... definition of a finite element method mf_rhs
    // ... definition of an integration method mim
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

