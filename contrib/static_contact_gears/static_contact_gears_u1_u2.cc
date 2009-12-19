// -*- c++ -*- (enables emacs c++ mode)
//===========================================================================
//
// Copyright (C) 2009-2009 Konstantinos Poulios.
//
// This file is a part of GETFEM++
//
// Getfem++  is  free software;  you  can  redistribute  it  and/or modify it
// under  the  terms  of the  GNU  Lesser General Public License as published
// by  the  Free Software Foundation;  either version 2.1 of the License,  or
// (at your option) any later version.
// This program  is  distributed  in  the  hope  that it will be useful,  but
// WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
// or  FITNESS  FOR  A PARTICULAR PURPOSE.  See the GNU Lesser General Public
// License for more details.
// You  should  have received a copy of the GNU Lesser General Public License
// along  with  this program;  if not, write to the Free Software Foundation,
// Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301, USA.
//
//===========================================================================

#include "getfem/getfem_regular_meshes.h"
#include "getfem/getfem_model_solvers.h"
#include "getfem/getfem_import.h"   /* import functions (load a mesh from file)    */
#include "getfem/getfem_export.h"   /* export functions (save solution in a file)  */
#include "getfem/getfem_Coulomb_friction.h"
#include "gmm/gmm.h"
#include <map>

/* some Getfem++ types that we will be using */
using bgeot::dim_type;
using bgeot::size_type;   /* = unsigned long */
using bgeot::scalar_type; /* = double */
using bgeot::base_node;   /* geometrical nodes(derived from base_small_vector)*/
using bgeot::base_matrix; /* small dense matrix. */

/* definition of some matrix/vector types. These ones are built
 * using the predefined types in Gmm++
 */
typedef getfem::model_real_plain_vector  plain_vector;

/*
  structure for the elastostatic problem
*/
struct elastostatic_problem {

  enum { DIRICHLET_BOUNDARY = 0,
         DIRICHLET_BOUNDARY_1 = 1, DIRICHLET_BOUNDARY_2 = 2,
         CONTACT_BOUNDARY_1 = 3, CONTACT_BOUNDARY_2 = 4};
  getfem::mesh mesh1, mesh2; /* the meshes */
  getfem::mesh_im mim1, mim2;/* the integration methods */
  getfem::mesh_fem mf_u1;    /* 1st mesh_fem, for the elastostatic solution  */
  getfem::mesh_fem mf_u2;    /* 2nd mesh_fem, for the elastostatic solution  */
  getfem::mesh_fem mf_rhs1;  /* 1st mesh_fem for the right hand side         */
  getfem::mesh_fem mf_rhs2;  /* 2nd mesh_fem for the right hand side         */
  scalar_type lambda, mu;    /* elastic coefficients.                        */

  scalar_type residual;      /* max residual for the iterative solvers       */
  scalar_type rot_angle;     /* rotation angle of the pinion gear            */
  scalar_type threshold;     /* threshold distance for contact finding       */

  size_type N;               /* dimension of the problem                     */

  bool solve(plain_vector &, plain_vector &, plain_vector &,
             plain_vector &, plain_vector &, plain_vector &, plain_vector &);
  void init(void);
  elastostatic_problem(void) : mim1(mesh1), mim2(mesh2),
    mf_u1(mesh1), mf_u2(mesh2), mf_rhs1(mesh1), mf_rhs2(mesh2) {}
};


/* Define the problem parameters, import the mesh, set finite element
   and integration methods and detect the boundaries.
*/
void elastostatic_problem::init(void) {

  std::string FEM_TYPE  = "FEM_QK(3,1)";
  std::string INTEGRATION = "IM_HEXAHEDRON(5)";

  /* First step : import the mesh */
  getfem::import_mesh("gmsh:./gear1.msh", mesh1);
  getfem::import_mesh("gmsh:./gear2.msh", mesh2);

  for (size_type ii = 0; ii <= 1; ii++) {
    getfem::mesh &m = ii ? mesh2 : mesh1;
    getfem::mesh_region &cb = m.region(ii ? CONTACT_BOUNDARY_2 : CONTACT_BOUNDARY_1);
    getfem::mesh_region &db = m.region(ii ? DIRICHLET_BOUNDARY_2 : DIRICHLET_BOUNDARY_1);

    for (getfem::mr_visitor i(m.region(113)); !i.finished(); ++i)
      if (i.is_face()) cb.add(i.cv(), i.f());
    for (getfem::mr_visitor i(m.region(133)); !i.finished(); ++i)
      if (i.is_face()) db.add(i.cv(), i.f());
    for (getfem::mr_visitor i(m.region(142)); !i.finished(); ++i)
      if (i.is_face()) db.add(i.cv(), i.f());
    for (getfem::mr_visitor i(m.region(143)); !i.finished(); ++i)
      if (i.is_face()) db.add(i.cv(), i.f());
    for (getfem::mr_visitor i(m.region(173)); !i.finished(); ++i)
      if (i.is_face()) db.add(i.cv(), i.f());
    for (getfem::mr_visitor i(m.region(182)); !i.finished(); ++i)
      if (i.is_face()) db.add(i.cv(), i.f());
    for (getfem::mr_visitor i(m.region(183)); !i.finished(); ++i)
      if (i.is_face()) db.add(i.cv(), i.f());
  }

  N = mesh1.dim(); //=mesh2.dim()

  residual = 1e-6;
  rot_angle = -1.5e-2;
  threshold = 10.;

  lambda = 1.18e+5;
  mu = 0.83e+5;

  mf_u1.set_qdim(dim_type(N));
  mf_u2.set_qdim(dim_type(N));

  /* set the finite element on the mf_u */
  getfem::pfem pf_u = getfem::fem_descriptor(FEM_TYPE);
  mf_u1.set_finite_element(pf_u);
  mf_u2.set_finite_element(pf_u);

  getfem::pintegration_method ppi = getfem::int_method_descriptor(INTEGRATION);
  mim1.set_integration_method(ppi);
  mim2.set_integration_method(ppi);

  /* set the finite element on mf_rhs */
  std::string data_fem_name = "FEM_QK(3,1)";
  mf_rhs1.set_finite_element(mesh1.convex_index(),
                             getfem::fem_descriptor(data_fem_name));
  mf_rhs2.set_finite_element(mesh2.convex_index(),
                             getfem::fem_descriptor(data_fem_name));
}

/*  Construction and solution of the Model.
*/
bool elastostatic_problem::solve(plain_vector &U1, plain_vector &U2, plain_vector &RHS,
                                 plain_vector &Forces1, plain_vector &Forces2,
                                 plain_vector &CForces1, plain_vector &CForces2) {
  size_type nb_dof_rhs1 = mf_rhs1.nb_dof();
  size_type nb_dof_rhs2 = mf_rhs2.nb_dof();

  cout << "mf_rhs1.nb_dof() :" << mf_rhs1.nb_dof() << endl;
  cout << "mf_u1.nb_dof()   :" << mf_u1.nb_dof() << endl;
  cout << "mf_rhs2.nb_dof() :" << mf_rhs2.nb_dof() << endl;
  cout << "mf_u2.nb_dof()   :" << mf_u2.nb_dof() << endl;

  getfem::model md;
  md.add_fem_variable("u1", mf_u1);
  md.add_fem_variable("u2", mf_u2);

  // Linearized elasticity brick.
  md.add_initialized_scalar_data("lambda", lambda);
  md.add_initialized_scalar_data("mu", mu);
  getfem::add_isotropic_linearized_elasticity_brick(md, mim1, "u1", "lambda", "mu");
  getfem::add_isotropic_linearized_elasticity_brick(md, mim2, "u2", "lambda", "mu");

  // Nonlinear elasticity brick.
//  base_vector p(2); p[0] = lambda; p[1] = mu;
//  getfem::SaintVenant_Kirchhoff_hyperelastic_law pl;
//  getfem::mdbrick_nonlinear_elasticity<>  ELAS(pl, mim, mf_u, p);
  
  // Defining the contact condition.
  std::string varname_u1="u1";
  std::string varname_u2="u2";
  std::string dataname_r="r";
  md.add_initialized_scalar_data
    (dataname_r, mu * (3*lambda + 2*mu) / (lambda + mu) );  // r ~= Young modulus
  std::string multname_n;
  getfem::add_frictionless_contact_brick
    (md, mim1, mim2, varname_u1, varname_u2, multname_n, dataname_r,
     CONTACT_BOUNDARY_1, CONTACT_BOUNDARY_2);

  // Defining the DIRICHLET condition.
  plain_vector F(nb_dof_rhs1 * N);
  dal::bit_vector
    cn = mf_rhs1.basic_dof_on_region(mesh1.region(DIRICHLET_BOUNDARY_1));
  for (dal::bv_visitor i(cn); !i.finished(); ++i) {
    base_node node = mf_rhs1.point_of_basic_dof(i);
    F[i*N] = -node[1] * rot_angle;
    F[i*N+1] = node[0] * rot_angle;
  }
  md.add_initialized_fem_data("DirichletData1", mf_rhs1, F);
  getfem::add_Dirichlet_condition_with_multipliers
    (md, mim1, "u1", mf_u1, DIRICHLET_BOUNDARY_1, "DirichletData1");

  gmm::resize(F, nb_dof_rhs2 * N);
  gmm::clear(F);
  md.add_initialized_fem_data("DirichletData2", mf_rhs2, F);
  getfem::add_Dirichlet_condition_with_multipliers
    (md, mim2, "u2", mf_u2, DIRICHLET_BOUNDARY_2, "DirichletData2");

  // Defining the surface pressure term for the NEUMANN boundary.
//  base_vector f(N);
//  for (size_type i = 0; i < nb_dof_rhs; ++i) {
//    gmm::copy(f, gmm::sub_vector(F, gmm::sub_interval(i*N, N)));
//  }
//  getfem::mdbrick_source_term<> NEUMANN(FRICTION, mf_rhs, F, NEUMANN_BOUNDARY_NUM);

  gmm::iteration iter(residual, 1, 40000);

  gmm::default_newton_line_search ls;
  getfem::standard_solve(md, iter, getfem::rselect_linear_solver(md,"superlu"), ls);

  gmm::resize(U1, mf_u1.nb_dof());
  gmm::resize(U2, mf_u2.nb_dof());
  gmm::resize(RHS, md.nb_dof());
  gmm::resize(CForces1, mf_u1.nb_dof());
  gmm::resize(CForces2, mf_u2.nb_dof());

  gmm::copy(md.real_variable("u1"), U1);
  gmm::copy(md.real_variable("u2"), U2);
  gmm::copy(md.real_rhs(), RHS);

  gmm::copy(gmm::sub_vector(RHS, md.interval_of_variable("u1")), Forces1);
  gmm::scale(Forces1, -1.0);
  gmm::copy(gmm::sub_vector(RHS, md.interval_of_variable("u2")), Forces2);
  gmm::scale(Forces2, -1.0);

  gmm::mult(gmm::sub_matrix(md.real_tangent_matrix(),
                            md.interval_of_variable("u1"),
                            md.interval_of_variable(multname_n) ),
            md.real_variable(multname_n),
            CForces1);
  gmm::scale(CForces1, -1.0);
  gmm::mult(gmm::sub_matrix(md.real_tangent_matrix(),
                            md.interval_of_variable("u2"),
                            md.interval_of_variable(multname_n) ),
            md.real_variable(multname_n),
            CForces2);
  gmm::scale(CForces2, -1.0);

  return (iter.converged());
}

/**************************************************************************/
/*  main program.                                                         */
/**************************************************************************/

int main(int argc, char *argv[]) {

  elastostatic_problem p;
  p.init();

  plain_vector U1(p.mf_u1.nb_dof()), U2(p.mf_u2.nb_dof());
  plain_vector RHS;
  plain_vector Forces1(p.mf_u1.nb_dof()), Forces2(p.mf_u2.nb_dof());
  plain_vector CForces1(p.mf_u1.nb_dof()), CForces2(p.mf_u2.nb_dof());
  if (!p.solve(U1, U2, RHS, Forces1, Forces2, CForces1, CForces2))
    cout << "Solve has failed\n";

  p.mesh1.write_to_file("static_contact_gears1.mesh");
  p.mesh2.write_to_file("static_contact_gears2.mesh");
  p.mf_u1.write_to_file("static_contact_gears1.mf", true);
  p.mf_u2.write_to_file("static_contact_gears2.mf", true);
  p.mf_rhs1.write_to_file("static_contact_gears1.mfd", true);
  p.mf_rhs2.write_to_file("static_contact_gears2.mfd", true);
  gmm::vecsave("static_contact_gears1.U", U1);
  gmm::vecsave("static_contact_gears2.U", U2);
  gmm::vecsave("static_contact_gears.RHS", RHS);
  getfem::vtk_export exp1("static_contact_gears1.vtk", true);
  exp1.exporting(p.mf_u1);
  exp1.write_point_data(p.mf_u1, U1, "elastostatic_displacement_1");
  exp1.write_point_data(p.mf_u1, Forces1, "forces_1");
  exp1.write_point_data(p.mf_u1, CForces1, "contact_forces_1");
  getfem::vtk_export exp2("static_contact_gears2.vtk", true);
  exp2.exporting(p.mf_u2);
  exp2.write_point_data(p.mf_u2, U2, "elastostatic_displacement_2");
  exp2.write_point_data(p.mf_u2, Forces2, "forces_2");
  exp2.write_point_data(p.mf_u2, CForces2, "contact_forces_2");

  return 0;
}

