/*===========================================================================

 Copyright (C) 2009-2020 Konstantinos Poulios.

 This file is a part of GetFEM

 GetFEM  is  free software;  you  can  redistribute  it  and/or modify it
 under  the  terms  of the  GNU  Lesser General Public License as published
 by  the  Free Software Foundation;  either version 3 of the License,  or
 (at your option) any later version along with the GCC Runtime Library
 Exception either version 3.1 or (at your option) any later version.
 This program  is  distributed  in  the  hope  that it will be useful,  but
 WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 or  FITNESS  FOR  A PARTICULAR PURPOSE.  See the GNU Lesser General Public
 License and GCC Runtime Library Exception for more details.
 You  should  have received a copy of the GNU Lesser General Public License
 along  with  this program;  if not, write to the Free Software Foundation,
 Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301, USA.

===========================================================================*/

#include "getfem/getfem_model_solvers.h"
#include "getfem/getfem_import.h"   /* import functions (load a mesh from file)    */
#include "getfem/getfem_export.h"   /* export functions (save solution in a file)  */
#include "getfem/getfem_contact_and_friction_nodal.h"
#include "getfem/getfem_contact_and_friction_integral.h"
#include "getfem/getfem_contact_and_friction_large_sliding.h"
#include "gmm/gmm.h"

using std::endl; using std::cout; using std::cerr;
using std::ends; using std::cin;

/* some GetFEM types that we will be using */
using bgeot::dim_type;
using bgeot::size_type;   /* = unsigned long */
using bgeot::short_type;
using bgeot::scalar_type; /* = double */
using bgeot::base_node;   /* geometrical nodes (derived from base_small_vector)*/
using bgeot::base_matrix; /* small dense matrix. */

/* definition of some matrix/vector types. These ones are built
 * using the predefined types in Gmm++
 */
typedef getfem::model_real_plain_vector  plain_vector;

/*
  structure for the elastostatic problem
*/
struct elastostatic_contact_problem {

  enum { DIRICHLET_BOUNDARY_1 = 0, DIRICHLET_BOUNDARY_2 = 1,
         CONTACT_BOUNDARY_1 = 1001, CONTACT_BOUNDARY_2 = 1002 };
  getfem::mesh mesh;           /* the mesh */
  getfem::mesh_im mim;         /* the integration methods */
  getfem::mesh_fem mf_u;       /* main mesh_fem, for the elastostatic solution */
  getfem::mesh_fem mf_rhs;     /* mesh_fem for the right hand side             */
  getfem::mesh_fem mf_mult;    /* mesh_fem for the multipliers                 */
  scalar_type lambda, mu;      /* elastic coefficients.                        */

  scalar_type residual;        /* max residual for the iterative solvers       */
  scalar_type rot_angle;       /* rotation angle of the pinion gear            */
  scalar_type frict_coeff;     /* friction coefficient                         */
//scalar_type threshold;       /* threshold distance for contact finding       */

  size_type N;                 /* dimension of the problem                     */

  bool frictionless;           /* flag for frictionless model                  */
  int contact_algo;            /* contact algorithm (0:nodal, 1-4: integral
                                  >=5:integral large sliding)                  */

  // Vectors holding the ids of mesh region pairs expected to come in contact
  // with each other
  std::vector<size_type> cb_rgs1, cb_rgs2;
  size_type max_rg;            /* maximum id used to define a region in the mesh */

  std::string datafilename;
  bgeot::md_param PARAM;

  bool solve(void);
  void init(void);
  elastostatic_contact_problem(void)
    : mim(mesh), mf_u(mesh), mf_rhs(mesh), mf_mult(mesh), max_rg(1) {}

};


/* Define the problem parameters, import the mesh, set finite element
   and integration methods and detect the boundaries.
*/
void elastostatic_contact_problem::init(void) {

  std::string FEM_TYPE  = PARAM.string_value("FEM_TYPE","FEM name");
  std::string INTEGRATION = PARAM.string_value("INTEGRATION",
                                               "Name of integration method");
  /* First step : import the mesh */
  std::string meshname_1 = PARAM.string_value("MESHNAME_GEAR1",
                                              "Mesh filename for the 1st gear");
  std::string meshname_2 = PARAM.string_value("MESHNAME_GEAR2",
                                              "Mesh filename for the 2nd gear");
  size_type nb_cf_pairs;
  for (size_type swap = 0; swap <= 1; swap++) {
    std::vector<size_type> &cb_rgs = swap ? cb_rgs2 : cb_rgs1;
    // Contact faces
    const std::vector<bgeot::md_param::param_value> &cf_params
      = PARAM.array_value( swap ? "CONTACT_FACES_2" : "CONTACT_FACES_1" );
    if (swap == 0) nb_cf_pairs = cf_params.size();
    else GMM_ASSERT1(nb_cf_pairs == cf_params.size(),
                     "Vectors CONTACT_FACES_1 and CONTACT_FACES_2 " <<
                     "should have the same size.");
    std::vector<size_type> cfs(nb_cf_pairs);
    for (size_type i = 0; i < cfs.size(); ++i) {
      GMM_ASSERT1(cf_params[i].type_of_param() == bgeot::md_param::REAL_VALUE,
                  (swap ? "CONTACT_FACES_2" :  "CONTACT_FACES_1")
                   << " should be an integer array.");
      cfs[i] = size_type(cf_params[i].real()+0.5);
      cb_rgs.push_back(++max_rg);
    }
    // Dirichlet faces
    const std::vector<bgeot::md_param::param_value> &df_params
      = PARAM.array_value(swap ? "DIRICHLET_FACES_2" : "DIRICHLET_FACES_1");
    std::vector<size_type> dfs(df_params.size());
    for (size_type i = 0; i < dfs.size(); ++i) {
      GMM_ASSERT1(df_params[i].type_of_param() == bgeot::md_param::REAL_VALUE,
                  (swap ? "DIRICHLET_FACES_2" : "DIRICHLET_FACES_1")
                   << " should be an integer array.");
      dfs[i] = size_type(df_params[i].real()+0.5);
    }

    getfem::mesh tmpmesh;
    getfem::import_mesh(swap ? meshname_2 : meshname_1, tmpmesh);

    N = tmpmesh.dim();
    base_node Pmin(N), Pmax(N);
    tmpmesh.bounding_box(Pmin, Pmax);
    bool reduce_dim = false;
    if (gmm::abs(Pmax[N-1] - Pmin[N-1]) < 1.E-10) reduce_dim = true;

    getfem::mesh_region dirichlet_boundary =
      mesh.region(swap ? DIRICHLET_BOUNDARY_2 : DIRICHLET_BOUNDARY_1);

    // Add convexes to the main mesh
    for (dal::bv_visitor cv(tmpmesh.convex_index()); !cv.finished(); ++cv) {
      size_type newcv;
      if (reduce_dim) {
        std::vector<base_node> pt_tab;
        for (size_type i=0; i < tmpmesh.nb_points_of_convex(cv); ++i) {
          base_node PP = tmpmesh.points_of_convex(cv)[i];
          base_node P(N-1);
          for (size_type j=0; j < N-1; ++j) P[j] = PP[j];
          pt_tab.push_back(P);
        }
        newcv =
          mesh.add_convex_by_points(tmpmesh.trans_of_convex(cv), pt_tab.begin());
      } else {
        newcv =
          mesh.add_convex_by_points(tmpmesh.trans_of_convex(cv),
                                    tmpmesh.points_of_convex(cv).begin());
      }
      for (short_type f = 0; f < tmpmesh.structure_of_convex(cv)->nb_faces(); f++) {
        for (size_type i = 0; i < dfs.size(); ++i)
          if (tmpmesh.region(dfs[i]).is_in(cv,f)) {
            dirichlet_boundary.add(newcv,f);
            break;
          }
        for (size_type i = 0; i < cfs.size(); ++i)
          if (tmpmesh.region(cfs[i]).is_in(cv,f)) {
            GMM_ASSERT1(!dirichlet_boundary.is_in(newcv,f),
                        "Overlapping contact and dirichlet regions");
            mesh.region(cb_rgs[i]).add(newcv,f);
          }
      }
    }
  }

  N = mesh.dim();

  residual = PARAM.real_value("RESIDUAL", "Residual for Newton solver");
  if (residual == 0.) residual = 1e-10;

  rot_angle = PARAM.real_value("ROT_ANGLE", "Rotation angle of the first gear");
  frict_coeff = PARAM.real_value("FRICTION_COEFFICIENT", "Friction coefficient");
  if (frict_coeff == 0.) frictionless = true;
  else frictionless = false;
//  threshold = 10.;  //FIXME

  mu = PARAM.real_value("MU", "Lamé coefficient mu");
  lambda = PARAM.real_value("LAMBDA", "Lamé coefficient lambda");

  mf_u.set_qdim(dim_type(N));

  /* set the finite element on the mf_u */
  getfem::pfem pf_u = getfem::fem_descriptor(FEM_TYPE);
  mf_u.set_finite_element(pf_u);
  GMM_ASSERT1(pf_u->is_lagrange(), "You are using a non-lagrange FEM. "
                                << "For this problem you need a lagrange FEM.");

  getfem::pintegration_method ppi = getfem::int_method_descriptor(INTEGRATION);
  mim.set_integration_method(ppi);

  /* set the finite element on mf_rhs */
  mf_rhs.set_finite_element(getfem::fem_descriptor(FEM_TYPE));

  datafilename = PARAM.string_value("ROOTFILENAME","Base name of data files.");

  contact_algo = int(PARAM.int_value("CONTACT_ALGO","Algorithm for imposing the contact condition."));

  if (contact_algo != 0) { // integral contact
    std::string MULT_FEM_TYPE  = PARAM.string_value("MULT_FEM_TYPE","FEM name for the multipliers");
    if (frictionless && contact_algo >= 1 && contact_algo <= 4)
      mf_mult.set_qdim(dim_type(1));
    else
      mf_mult.set_qdim(dim_type(N));
    getfem::pfem pf_mult = getfem::fem_descriptor(MULT_FEM_TYPE);
    mf_mult.set_finite_element(pf_mult);

    getfem::mesh_region &mr1 = mesh.region(CONTACT_BOUNDARY_1);
    getfem::mesh_region &mr2 = mesh.region(CONTACT_BOUNDARY_2);
    for (std::vector<size_type>::const_iterator rg_it=cb_rgs1.begin();
         rg_it != cb_rgs1.end(); rg_it++)
      mr1 = getfem::mesh_region::merge(mr1, mesh.region(*rg_it));
    for (std::vector<size_type>::const_iterator rg_it=cb_rgs2.begin();
         rg_it != cb_rgs2.end(); rg_it++)
      mr2 = getfem::mesh_region::merge(mr2, mesh.region(*rg_it));

    dal::bit_vector dol = mf_mult.basic_dof_on_region(CONTACT_BOUNDARY_1);
    mf_mult.reduce_to_basic_dof(dol);
  }

}

/*  Construction and solution of the Model.
*/
bool elastostatic_contact_problem::solve() {
  size_type nb_dof_rhs = mf_rhs.nb_dof();

  cout << "mf_u.nb_dof()   :" << mf_u.nb_dof() << endl;

  getfem::model md;
  md.add_fem_variable("u", mf_u);

  // Linearized elasticity brick.
  md.add_initialized_scalar_data("lambda", lambda);
  md.add_initialized_scalar_data("mu", mu);
  getfem::add_isotropic_linearized_elasticity_brick(md, mim, "u", "lambda", "mu");

  // Defining the contact condition.
  md.add_initialized_scalar_data
    ("r", mu * (3*lambda + 2*mu) / (lambda + mu) );  // r ~= Young modulus
  std::string multname_n, multname_t;
  if (contact_algo == 0) {
    if (frictionless) {
      getfem::add_nodal_contact_between_nonmatching_meshes_brick
        (md, mim, "u", multname_n, "r", cb_rgs1, cb_rgs2);
    } else {
      std::string dataname_frict_coeff="friction_coefficient";
      md.add_initialized_scalar_data(dataname_frict_coeff, frict_coeff);
      getfem::add_nodal_contact_between_nonmatching_meshes_brick
        (md, mim, "u", multname_n, multname_t,
         "r", dataname_frict_coeff, cb_rgs1, cb_rgs2);
    }
  } else {
    md.add_fem_variable("mult", mf_mult);
    if (contact_algo >= 1 && contact_algo <= 4) { // integral contact
      if (frictionless)
        getfem::add_integral_contact_between_nonmatching_meshes_brick
          (md, mim, "u", "u", "mult", "r",
           CONTACT_BOUNDARY_1, CONTACT_BOUNDARY_2, contact_algo);
      else {
        md.add_initialized_scalar_data("f_coeff", frict_coeff);
        getfem::add_integral_contact_between_nonmatching_meshes_brick
          (md, mim, "u", "u", "mult", "r", "f_coeff",
           CONTACT_BOUNDARY_1, CONTACT_BOUNDARY_2, contact_algo);
      }
    }
    else { // large sliding is for the moment always frictionless
      std::string u0_str("");
//      if (frict_coeff > scalar_type(0)) {
//        u0_str = "u0";
//        md.add_fem_variable(u0_str, mf_u);
//      }
      md.add_initialized_scalar_data("f_coeff", frict_coeff);
      size_type indb =
      getfem::add_integral_large_sliding_contact_brick_raytracing
      (md, "r", 20., "f_coeff", "1", false, false);
      getfem::add_contact_boundary_to_large_sliding_contact_brick
      (md, indb, mim, CONTACT_BOUNDARY_1, false, true, "u", "mult", u0_str);
      getfem::add_contact_boundary_to_large_sliding_contact_brick
      (md, indb, mim, CONTACT_BOUNDARY_2, true, false, "u", "", u0_str);
    }
  }

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

  cout << "md.nb_dof()   :" << md.nb_dof() << endl;

  gmm::iteration iter(residual, 1, 40000);

//  getfem::default_newton_line_search ls;
  getfem::simplest_newton_line_search ls(50, 5., 5., 0.6, 1e-1);

  #if defined(GMM_USES_MUMPS)
    getfem::standard_solve(md, iter, getfem::rselect_linear_solver(md,"mumps"), ls);
  #else
    getfem::standard_solve(md, iter, getfem::default_linear_solver<getfem::model_real_sparse_matrix, getfem::model_real_plain_vector>(md), ls);
  #endif

  if (!iter.converged()) return false; // Solution has not converged

  plain_vector VM(mf_rhs.nb_dof());
  getfem::compute_isotropic_linearized_Von_Mises_or_Tresca
      (md, "u", "lambda", "mu", mf_rhs, VM, false);

  if (!getfem::MPI_IS_MASTER())
    return true;

  // Prepare results
  plain_vector U(mf_u.nb_dof());
  plain_vector RHS(md.nb_dof());
  plain_vector Forces(mf_u.nb_dof());
  plain_vector NCForces(mf_u.nb_dof());
  plain_vector TCForces(0);

  gmm::copy(md.real_variable("u"), U);
  gmm::copy(md.real_rhs(), RHS);

  gmm::copy(gmm::sub_vector(RHS, md.interval_of_variable("u")), Forces);
  gmm::scale(Forces, -1.0);

  // Export results
  mesh.write_to_file(datafilename + ".mesh");
  mf_u.write_to_file(datafilename + ".mf", true);
  mf_rhs.write_to_file(datafilename + ".mfd", true);
  gmm::vecsave(datafilename + ".U", U);
  gmm::vecsave(datafilename + ".VM", VM);
  gmm::vecsave(datafilename + ".RHS", RHS);
  getfem::vtk_export exp(datafilename + ".vtk", true);
  exp.exporting(mesh);
  exp.write_point_data(mf_u, U, "elastostatic_displacement");
  exp.write_point_data(mf_u, Forces, "forces");
  exp.write_point_data(mf_u, NCForces, "normal_contact_forces");
  exp.write_point_data(mf_rhs, VM, "von_mises_stresses");

  if (contact_algo == 0) {
    gmm::mult(gmm::sub_matrix(md.real_tangent_matrix(),
                              md.interval_of_variable("u"),
                              md.interval_of_variable(multname_n) ),
              md.real_variable(multname_n),
              NCForces);
    gmm::scale(NCForces, -1.0);

    if (!frictionless) {
      gmm::resize(TCForces, mf_u.nb_dof());
      gmm::mult(gmm::sub_matrix(md.real_tangent_matrix(),
                                md.interval_of_variable("u"),
                                md.interval_of_variable(multname_t) ),
                md.real_variable(multname_t),
                TCForces);
      gmm::scale(TCForces, -1.0);
      exp.write_point_data(mf_u, TCForces, "tangential_contact_forces");
    }
  }

  return true; // Solution has converged
}

/**************************************************************************/
/*  main program.                                                         */
/**************************************************************************/

int main(int argc, char *argv[]) {

  GETFEM_MPI_INIT(argc, argv); // For parallelized version

  elastostatic_contact_problem p;
  p.PARAM.read_command_line(argc, argv);
  p.init();
  if (!p.solve()) cout << "Solve has failed\n";

  GETFEM_MPI_FINALIZE;

  return 0;
}
