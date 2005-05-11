/**
 * Linear Elastostatic plate problem.
 *
 * This program is used to check that getfem++ is working. This is also 
 * a good example of use of Getfem++.
*/

#include <getfem_assembling.h> /* import assembly methods (and norms comp.) */
#include <getfem_linearized_plates.h>
#include <getfem_export.h>   /* export functions (save solution in a file)  */
#include <getfem_regular_meshes.h>
#include <getfem_modeling.h>
#include <gmm.h>
#include <getfem_derivatives.h>
#include <getfem_mesh_im_level_set.h>
#include <getfem_mesh_fem_level_set.h>
#include <getfem_mesh_fem_product.h>
#include <getfem_mesh_fem_global_function.h>
#include <getfem_mesh_fem_sum.h>

/* try to enable the SIGFPE if something evaluates to a Not-a-number
 * of infinity during computations
 */
#ifdef GETFEM_HAVE_FEENABLEEXCEPT
#  include <fenv.h>
#endif

/* some Getfem++ types that we will be using */
using bgeot::base_small_vector; /* special class for small (dim<16) vectors */
using bgeot::base_node;  /* geometrical nodes(derived from base_small_vector)*/
using bgeot::scalar_type; /* = double */
using bgeot::size_type;   /* = unsigned long */
using bgeot::base_matrix; /* small dense matrix. */

/* definition of some matrix/vector types. These ones are built
 * using the predefined types in Gmm++
 */
typedef getfem::modeling_standard_sparse_vector sparse_vector;
typedef getfem::modeling_standard_sparse_matrix sparse_matrix;
typedef getfem::modeling_standard_plain_vector  plain_vector; 

struct crackPlate_problem{

  enum { DIRICHLET_BOUNDARY_NUM = 0, NEUMANN_BOUNDARY_NUM = 1};
  getfem::getfem_mesh mesh;  /* the mesh */
  getfem::mesh_level_set mls;       /* the integration methods.              */
  getfem::mesh_im_level_set mim, mim_subint;    /* the integration methods.              */
  getfem::mesh_fem mf_pre_ut, mf_pre_u3, mf_pre_theta ; 
  getfem::mesh_fem_level_set mfls_ut, mfls_u3, mfls_theta ; 
  getfem::mesh_fem_global_function mf_sing_ut, mf_sing_u3, mf_sing_theta ;
  getfem::mesh_fem mf_partition_of_unity;
  getfem::mesh_fem_product mf_ut_product, mf_u3_product, mf_theta_product;
  getfem::mesh_fem_sum mf_ut_sum, mf_u3_sum, mf_theta_sum;

  getfem::mesh_fem& mf_ut() { return mf_ut_sum; }
  getfem::mesh_fem& mf_u3() { return mf_u3_sum; }
  getfem::mesh_fem& mf_theta() { return mf_theta_sum; }

  getfem::mesh_fem mf_rhs;   /* mesh_fem for the right hand side (f(x),..)   */
  getfem::mesh_fem mf_coef;  /* mesh_fem used to represent pde coefficients  */
  getfem::level_set ls;      /* The two level sets defining the crack.       */
 
  scalar_type residu;       /* max residu for the iterative solvers        */
  scalar_type cutoff_radius, enr_area_radius;
  int enrichment_option;
  std::string datafilename;
  ftool::md_param PARAM;
  size_type mitc ;

  scalar_type lambda, mu;    /* Lamé coefficients.                           */
  scalar_type epsilon, pressure ;
  // methods
  bool solve(plain_vector &U);
  void init(void);
  crackPlate_problem(void) : mls(mesh), mim(mls), mim_subint(mls),
  			mf_pre_ut(mesh), mf_pre_u3(mesh), mf_pre_theta(mesh), 
			mfls_ut(mls, mf_pre_ut), mfls_u3(mls, mf_pre_u3), mfls_theta(mls, mf_pre_theta),
			mf_sing_ut(mesh),mf_sing_u3(mesh),mf_sing_theta(mesh),
			mf_partition_of_unity(mesh),
			mf_ut_product(mf_partition_of_unity, mf_sing_ut),
			mf_u3_product(mf_partition_of_unity, mf_sing_u3),
			mf_theta_product(mf_partition_of_unity, mf_sing_theta),
			mf_ut_sum(mesh), mf_u3_sum(mesh), mf_theta_sum(mesh),
			mf_rhs(mesh),  mf_coef(mesh),
			ls(mesh, 1, true) {} 

};

void crackPlate_problem::init(void) {
  
  const char *MESH_TYPE = PARAM.string_value("MESH_TYPE","Mesh type ");
  const char *FEM_TYPE_UT  = PARAM.string_value("FEM_TYPE_UT","FEM name");
  const char *FEM_TYPE_U3  = PARAM.string_value("FEM_TYPE_U3","FEM name");
  const char *FEM_TYPE_THETA = PARAM.string_value("FEM_TYPE_THETA","FEM name");
  const char *INTEGRATION = PARAM.string_value("INTEGRATION",
					       "Name of integration method");
  const char *INTEGRATION_CT = PARAM.string_value("INTEGRATION_CT",
					       "Name of integration method");
  const char *SIMPLEX_INTEGRATION = PARAM.string_value("SIMPLEX_INTEGRATION",
					 "Name of simplex integration method");
  enrichment_option = PARAM.int_value("ENRICHMENT_OPTION",
				      "Enrichment option");

  cout << "MESH_TYPE=" << MESH_TYPE << "\n";
  cout << "FEM_TYPE_UT="  << FEM_TYPE_UT << "\n";
  cout << "FEM_TYPE_U3="  << FEM_TYPE_U3 << "\n";
  cout << "FEM_TYPE_THETA="  << FEM_TYPE_THETA << "\n";
  cout << "INTEGRATION=" << INTEGRATION << "\n";
  cout << "INTEGRATION_CT=" << INTEGRATION_CT << "\n";
  
  // build the mesh :
  bgeot::pgeometric_trans pgt = 
  bgeot::geometric_trans_descriptor(MESH_TYPE);
  size_type N = pgt->dim();
  if (N != 2)
    DAL_THROW(getfem::failure_error, "For a plate problem, N should be 2");
  std::vector<size_type> nsubdiv(N);
  std::fill(nsubdiv.begin(),nsubdiv.end(),
	    PARAM.int_value("NX", "Number of space steps "));
  getfem::regular_unit_mesh(mesh, nsubdiv, pgt,
			    PARAM.int_value("MESH_NOISED") != 0);

  base_small_vector tt(N); tt[1] = -0.5;
  mesh.translation(tt); 
  
  // setting parameters :
  epsilon = PARAM.real_value("EPSILON", "thickness") ; 
  mitc = (PARAM.int_value("MITC") != 0);
  pressure = PARAM.real_value("PRESSURE") ;
  datafilename = PARAM.string_value("ROOTFILENAME","Base name of data files.");
  residu = PARAM.real_value("RESIDU"); if (residu == 0.) residu = 1e-10;
  enr_area_radius = PARAM.real_value("RADIUS_ENR_AREA",
				     "radius of the enrichment area");
  mu = PARAM.real_value("MU", "Lamé coefficient mu");
  lambda = PARAM.real_value("LAMBDA", "Lamé coefficient lambda");
  cutoff_radius = PARAM.real_value("CUTOFF", "Cutoff");
  
  mf_ut().set_qdim(N);
  mf_theta().set_qdim(N);
  
  // set the finite element method on ut, u3, theta
  getfem::pfem pf_ut = getfem::fem_descriptor(FEM_TYPE_UT);
  getfem::pfem pf_u3 = getfem::fem_descriptor(FEM_TYPE_U3);
  getfem::pfem pf_theta = getfem::fem_descriptor(FEM_TYPE_THETA);
  getfem::pintegration_method ppi = 
    getfem::int_method_descriptor(INTEGRATION);
  getfem::pintegration_method sppi = 
    getfem::int_method_descriptor(SIMPLEX_INTEGRATION);
  getfem::pintegration_method ppi_ct = 
    getfem::int_method_descriptor(INTEGRATION_CT) ;
    
  mim.set_integration_method(mesh.convex_index(), ppi);
  mim_subint.set_integration_method(mesh.convex_index(), ppi_ct) ;
  mls.add_level_set(ls);
  mim.set_simplex_im(sppi);
  mf_pre_ut.set_finite_element(mesh.convex_index(), pf_ut);
  mf_pre_u3.set_finite_element(mesh.convex_index(), pf_u3);
  mf_pre_theta.set_finite_element(mesh.convex_index(), pf_theta);
  mf_partition_of_unity.set_classical_finite_element(1);

  /* set the finite element on mf_rhs (same as mf_u is DATA_FEM_TYPE is
     not used in the .param file */
  const char *data_fem_name = PARAM.string_value("DATA_FEM_TYPE");
  if (data_fem_name == 0) {
    if (!pf_ut->is_lagrange()) {
      DAL_THROW(dal::failure_error, "You are using a non-lagrange FEM. "
		<< "In that case you need to set "
		<< "DATA_FEM_TYPE in the .param file");
    }
    mf_rhs.set_finite_element(mesh.convex_index(), pf_ut);
  } else {
    mf_rhs.set_finite_element(mesh.convex_index(), 
			      getfem::fem_descriptor(data_fem_name));
  }
  
  /* set the finite element on mf_coef. Here we use a very simple element
   *  since the only function that need to be interpolated on the mesh_fem 
   * is f(x)=1 ... */
  mf_coef.set_finite_element(mesh.convex_index(),
			     getfem::classical_fem(pgt,0));

  /* set boundary conditions : Dirichlet on the right face */
  cout << "Selecting Neumann and Dirichlet boundaries\n";
  getfem::convex_face_ct border_faces;
  getfem::outer_faces_of_mesh(mesh, border_faces);
  for (getfem::convex_face_ct::const_iterator it = border_faces.begin();
       it != border_faces.end(); ++it) {
    assert(it->f != size_type(-1));
    base_node un = mesh.normal_of_face_of_convex(it->cv, it->f);
    un /= gmm::vect_norm2(un);
    if ( un[0] >= 1. - 1.0E-7) { // new Dirichlet
      mesh.add_face_to_set(DIRICHLET_BOUNDARY_NUM, it->cv, it->f);
      } else {
      mesh.add_face_to_set(NEUMANN_BOUNDARY_NUM, it->cv, it->f);
      }
  }

}

bool crackPlate_problem::solve(plain_vector &U) {

  cout << "solving\n" ;
  size_type nb_dof_rhs = mf_rhs.nb_dof();
  size_type N = mesh.dim();

  ls.reinit();  
  for (size_type d = 0; d < ls.get_mesh_fem().nb_dof(); ++d) {
    ls.values(0)[d] = (ls.get_mesh_fem().point_of_dof(d))[1];
    ls.values(1)[d] = -0.5 + (ls.get_mesh_fem().point_of_dof(d))[0];
  }
  ls.touch();
  
  mls.adapt();
  mim.adapt();
  mim_subint.adapt();
  mfls_ut.adapt();
  mfls_u3.adapt();
  mfls_theta.adapt();
  
  cout << "setting singularties... \n" ;
  std::vector<getfem::pglobal_function> utfunc(4) ;
  //std::vector<getfem::pglobal_function> u3func(5) ;
  std::vector<getfem::pglobal_function> u3func(4) ;
  std::vector<getfem::pglobal_function> theta_func(4) ;
  
  // A mettre a jour !!!
  /////////////////////////////////////////////////////////////////////
  for (size_type i = 0 ; i < 4 ; ++i){
        utfunc[i] = isotropic_crack_singular_2D(i, ls,
			 (enrichment_option == 2) ? 0.0 : cutoff_radius);
	u3func[i] = isotropic_crack_singular_2D(i, ls,
			 (enrichment_option == 2) ? 0.0 : cutoff_radius);
	theta_func[i] = isotropic_crack_singular_2D(i, ls,
			 (enrichment_option == 2) ? 0.0 : cutoff_radius);
  }
  /////////////////////////////////////////////////////////////////////:
  
  mf_sing_ut.set_functions(utfunc);
  mf_sing_u3.set_functions(u3func);
  mf_sing_theta.set_functions(theta_func);
  
  if (enrichment_option == 2) {
    dal::bit_vector enriched_dofs;
    plain_vector X(mf_partition_of_unity.nb_dof());
    plain_vector Y(mf_partition_of_unity.nb_dof());
    getfem::interpolation(ls.get_mesh_fem(), mf_partition_of_unity,
			  ls.values(1), X);
    getfem::interpolation(ls.get_mesh_fem(), mf_partition_of_unity,
			  ls.values(0), Y);
    for (size_type j = 0; j < mf_partition_of_unity.nb_dof(); ++j) {
      if (gmm::sqr(X[j]) + gmm::sqr(Y[j]) <= gmm::sqr(enr_area_radius))
	enriched_dofs.add(j);
    }
    if (enriched_dofs.card() < 3)
      DAL_WARNING(0, "There is " << enriched_dofs.card() <<
		  " enriched dofs for the crack tip");
    mf_ut_product.set_enrichment(enriched_dofs);
    mf_u3_product.set_enrichment(enriched_dofs);
    mf_theta_product.set_enrichment(enriched_dofs);
    mf_ut_sum.set_mesh_fems(mf_ut_product, mfls_ut);
    mf_u3_sum.set_mesh_fems(mf_u3_product, mfls_u3);
    mf_theta_sum.set_mesh_fems(mf_theta_product, mfls_theta);
  }
  else { 
    mf_ut_sum.set_mesh_fems(mf_sing_ut, mfls_ut);
    mf_u3_sum.set_mesh_fems(mf_sing_u3, mfls_u3);
    mf_theta_sum.set_mesh_fems(mf_sing_theta, mfls_theta);
    }

  U.resize(mf_ut().nb_dof() + mf_u3().nb_dof() + mf_theta().nb_dof());
  
getfem::mdbrick_abstract<> *ELAS, *SIMPLE;

  // Linearized plate brick.
  getfem::mdbrick_isotropic_linearized_plate<>
    ELAS1(mim, mim_subint, mf_ut(), mf_u3(), mf_theta(), mf_coef, lambda,
	  mu, epsilon);
  if (mitc) ELAS1.set_mitc();  
  ELAS = &ELAS1;

  // Defining the surface source term.
  plain_vector F(nb_dof_rhs * 3); 
  plain_vector M(nb_dof_rhs * 2);
  for (size_type i = 0; i < nb_dof_rhs; ++i) F[3*i+2] = pressure;
  
  getfem::mdbrick_plate_source_term<> VOL_F(*ELAS, mf_rhs, F, M);
  
  getfem::mdbrick_plate_clamped_support<> SIMPLE1
    (VOL_F, mf_rhs, DIRICHLET_BOUNDARY_NUM, 0, 1);
       
  SIMPLE = &SIMPLE1 ;
  
  getfem::mdbrick_plate_closing<> final_model(*SIMPLE, 0, 1);
  
  // Generic solve.
  cout << "Total number of variables : " << final_model.nb_dof() << endl;
  getfem::standard_model_state MS(final_model);
  gmm::iteration iter(residu, 1, 40000);
  getfem::standard_solve(MS, final_model, iter);
  
  return (iter.converged());
}

/************************************************************
 * main program
 ************************************************************/
 
 
int main(int argc, char *argv[]) {
  dal::exception_callback_debug cb;
  dal::exception_callback::set_exception_callback(&cb); // to debug ...

#ifdef GETFEM_HAVE_FEENABLEEXCEPT /* trap SIGFPE */
  feenableexcept(FE_DIVBYZERO | FE_INVALID);
#endif

  // getfem::getfem_mesh_level_set_noisy();


  try {
    crackPlate_problem p;
    p.PARAM.read_command_line(argc, argv);
    p.init();
    p.mesh.write_to_file(p.datafilename + ".mesh") ;
    plain_vector U(p.mf_ut().nb_dof() + p.mf_u3().nb_dof() 
      + p.mf_theta().nb_dof());
    if (!p.solve(U)) DAL_THROW(dal::failure_error, "Solve has failed");
  }
  
    DAL_STANDARD_CATCH_ERROR;

  return 0; 
}



















