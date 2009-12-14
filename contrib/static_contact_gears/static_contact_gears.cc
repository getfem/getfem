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
         CONTACT_BOUNDARY_1 = 3, CONTACT_BOUNDARY_2 = 4,
         BLOCK_1 = 5, BLOCK_2 = 6};
  getfem::mesh mesh;         /* the mesh */
  getfem::mesh_im mim;       /* the integration methods */
  getfem::mesh_fem mf_u;     /* main mesh_fem, for the elastostatic solution */
  getfem::mesh_fem mf_rhs;   /* mesh_fem for the right hand side (f(x),..)   */
  scalar_type lambda, mu;    /* elastic coefficients.                        */

  scalar_type residual;      /* max residual for the iterative solvers       */
  scalar_type rot_angle;     /* rotation angle of the pinion gear            */
  scalar_type threshold;     /* threshold distance for contact finding       */

  size_type N;               /* dimension of the problem                     */

  bool solve(plain_vector &, plain_vector &, plain_vector &, plain_vector &);
  void init(void);
  elastostatic_problem(void) : mim(mesh), mf_u(mesh), mf_rhs(mesh) {}

};


/* Define the problem parameters, import the mesh, set finite element
   and integration methods and detect the boundaries.
*/
void elastostatic_problem::init(void) {

  std::string FEM_TYPE  = "FEM_QK(3,1)";
  std::string INTEGRATION = "IM_HEXAHEDRON(5)";

  /* First step : import the mesh */
  for (size_type ii = 0; ii <= 1; ii++) {
    getfem::mesh tmpmesh;
    getfem::mesh_region block, contact_boundary, dirichlet_boundary;
    if (ii == 0) {
      block = mesh.region(BLOCK_1);
      contact_boundary = mesh.region(CONTACT_BOUNDARY_1);
      dirichlet_boundary = mesh.region(DIRICHLET_BOUNDARY_1);
      getfem::import_mesh("gmsh:./gear1.msh",tmpmesh);
    } else {
      block = mesh.region(BLOCK_2);
      contact_boundary = mesh.region(CONTACT_BOUNDARY_2);
      dirichlet_boundary = mesh.region(DIRICHLET_BOUNDARY_2);
      getfem::import_mesh("gmsh:./gear2.msh",tmpmesh);
    }
    for (dal::bv_visitor cv(tmpmesh.convex_index()); !cv.finished(); ++cv) {
      size_type newcv =
        mesh.add_convex_by_points(tmpmesh.trans_of_convex(cv),
                                  tmpmesh.points_of_convex(cv).begin());
      block.add(newcv);
      for (size_type f = 0; f < tmpmesh.structure_of_convex(cv)->nb_faces(); f++) {
        if (tmpmesh.region(113).is_in(cv,f)) {
          contact_boundary.add(newcv,f);
        } else if (tmpmesh.region(133).is_in(cv,f) ||
                   tmpmesh.region(142).is_in(cv,f) ||
                   tmpmesh.region(143).is_in(cv,f) ||
                   tmpmesh.region(173).is_in(cv,f) ||
                   tmpmesh.region(182).is_in(cv,f) ||
                   tmpmesh.region(183).is_in(cv,f)  ) {
          dirichlet_boundary.add(newcv,f);
//          mesh.region(DIRICHLET_BOUNDARY).add(newcv,f);
        }
      }
    }
  }

  N = mesh.dim();

  residual = 1e-6;
  rot_angle = -1.5e-2;
  threshold = 10.;

  lambda = 1.18e+5;
  mu = 0.83e+5;

  mf_u.set_qdim(dim_type(N));

  /* set the finite element on the mf_u */
  getfem::pfem pf_u = getfem::fem_descriptor(FEM_TYPE);
  mf_u.set_finite_element(pf_u);

  getfem::pintegration_method ppi = getfem::int_method_descriptor(INTEGRATION);
  mim.set_integration_method(ppi);

  /* set the finite element on mf_rhs */
  std::string data_fem_name = "FEM_QK(3,1)";
  mf_rhs.set_finite_element(mesh.convex_index(), 
                            getfem::fem_descriptor(data_fem_name));

}

/*  Construction and solution of the Model.
*/
bool elastostatic_problem::solve(plain_vector &U, plain_vector &RHS,
                                 plain_vector &Forces, plain_vector &CForces) {
  size_type nb_dof_rhs = mf_rhs.nb_dof();

  cout << "mf_rhs.nb_dof() :" << mf_rhs.nb_dof() << endl;
  cout << "mf_u.nb_dof()   :" << mf_u.nb_dof() << endl;

  getfem::model md;
  md.add_fem_variable("u", mf_u);

  // Linearized elasticity brick.
  md.add_initialized_scalar_data("lambda", lambda);
  md.add_initialized_scalar_data("mu", mu);
  getfem::add_isotropic_linearized_elasticity_brick(md, mim, "u", "lambda", "mu");

  // Nonlinear elasticity brick.
//  base_vector p(2); p[0] = lambda; p[1] = mu;
//  getfem::SaintVenant_Kirchhoff_hyperelastic_law pl;
//  getfem::mdbrick_nonlinear_elasticity<>  ELAS(pl, mim, mf_u, p);
  
  // Defining the contact condition.
  std::string varname_u="u";
  std::string dataname_r="r";
  std::string dataname_alpha;
  md.add_initialized_scalar_data(dataname_r, 1.);
  getfem::add_frictionless_contact_brick
    (md, varname_u, CONTACT_BOUNDARY_1, CONTACT_BOUNDARY_2,
     dataname_r, dataname_alpha, false);

//cout << "no_cn: " << no_cn << endl;

  // Defining the DIRICHLET condition.
  plain_vector F(nb_dof_rhs * N);
  dal::bit_vector
    cn = mf_rhs.basic_dof_on_region(mesh.region(DIRICHLET_BOUNDARY_1));
  for (dal::bv_visitor i(cn); !i.finished(); ++i) {
    base_node node = mf_rhs.point_of_basic_dof(i);
    F[i*N] = -node[1] * rot_angle;
    F[i*N+1] = node[0] * rot_angle;
  }
  md.add_initialized_fem_data("DirichletData1", mf_rhs, F);
  getfem::add_Dirichlet_condition_with_multipliers
    (md, mim, "u", mf_u, DIRICHLET_BOUNDARY_1, "DirichletData1");

  gmm::clear(F);
  md.add_initialized_fem_data("DirichletData2", mf_rhs, F);
  getfem::add_Dirichlet_condition_with_multipliers
    (md, mim, "u", mf_u, DIRICHLET_BOUNDARY_2, "DirichletData2");

  // Defining the surface pressure term for the NEUMANN boundary.
//  base_vector f(N);
//  for (size_type i = 0; i < nb_dof_rhs; ++i) {
//    gmm::copy(f, gmm::sub_vector(F, gmm::sub_interval(i*N, N)));
//  }
//  getfem::mdbrick_source_term<> NEUMANN(FRICTION, mf_rhs, F, NEUMANN_BOUNDARY_NUM);

  gmm::iteration iter(residual, 1, 40000);

  gmm::default_newton_line_search ls;
  getfem::standard_solve(md, iter, getfem::rselect_linear_solver(md,"superlu"), ls);
                      
  gmm::resize(U, mf_u.nb_dof());
  gmm::resize(RHS, md.nb_dof());
  gmm::resize(CForces, mf_u.nb_dof());

  gmm::copy(md.real_variable("u"), U);
  gmm::copy(md.real_rhs(), RHS);

  gmm::copy(gmm::sub_vector(RHS, md.interval_of_variable("u")), Forces);
  gmm::scale(Forces, -1.0);

  gmm::mult(gmm::sub_matrix(md.real_tangent_matrix(),
                            md.interval_of_variable("u"),
                            md.interval_of_variable("contact_multiplier") ),
            md.real_variable("contact_multiplier"), 
            CForces);
  gmm::scale(CForces, -1.0);

  return (iter.converged());
}

/**************************************************************************/
/*  main program.                                                         */
/**************************************************************************/

int main(int argc, char *argv[]) {

  elastostatic_problem p;
  p.init();

  plain_vector U(p.mf_u.nb_dof());
  plain_vector RHS(p.mf_u.nb_dof());
  plain_vector Forces(p.mf_u.nb_dof());
  plain_vector CForces(p.mf_u.nb_dof());
  if (!p.solve(U, RHS, Forces, CForces)) cout << "Solve has failed\n";

  p.mesh.write_to_file("static_contact_gears.mesh");
  p.mf_u.write_to_file("static_contact_gears.mf", true);
  p.mf_rhs.write_to_file("static_contact_gears.mfd", true);
  gmm::vecsave("static_contact_gears.U", U);
  gmm::vecsave("static_contact_gears.RHS", RHS);
  getfem::vtk_export exp("static_contact_gears.vtk", true);
  exp.exporting(p.mf_u);
  exp.write_point_data(p.mf_u, U, "elastostatic_displacement");
  exp.write_point_data(p.mf_u, Forces, "forces");
  exp.write_point_data(p.mf_u, CForces, "contact_forces");

  return 0; 
}

