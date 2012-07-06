/*===========================================================================
 
 Copyright (C) 2002-2012 Yves Renard, Julien Pommier.
 
 This file is a part of GETFEM++
 
 Getfem++  is  free software;  you  can  redistribute  it  and/or modify it
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
  
/**
 * Linear Elastostatic problem with a crack.
 *
 * This program is used to check that getfem++ is working. This is also 
 * a good example of use of Getfem++.
 */

#include "getfem/getfem_assembling.h" /* import assembly methods (and norms comp.) */
#include "getfem/getfem_export.h"   /* export functions (save solution in a file)  */
#include "getfem/getfem_derivatives.h"
#include "getfem/getfem_regular_meshes.h"
#include "getfem/getfem_model_solvers.h"
#include "getfem/getfem_mesh_im_level_set.h"
#include "getfem/getfem_mesh_fem_level_set.h"
#include "getfem/getfem_mesh_fem_product.h"
#include "getfem/getfem_mesh_fem_global_function.h"
#include "getfem/getfem_spider_fem.h"
#include "getfem/getfem_mesh_fem_sum.h"
#include "gmm/gmm.h"
#include "gmm/gmm_inoutput.h"

#include "crack_exact_solution.h"

/* definition of some matrix/vector types. These ones are built
 * using the predefined types in Gmm++
 */
typedef getfem::modeling_standard_sparse_vector sparse_vector;
typedef getfem::modeling_standard_sparse_matrix sparse_matrix;
typedef getfem::modeling_standard_plain_vector  plain_vector;

/**************************************************************************/
/*  Structure for the crack problem.                                      */
/**************************************************************************/

struct crack_problem {

  enum { DIRICHLET_BOUNDARY_NUM = 0, NEUMANN_BOUNDARY_NUM = 1, NEUMANN_BOUNDARY_NUM1=2, NEUMANN_HOMOGENE_BOUNDARY_NUM=3, MORTAR_BOUNDARY_IN=42, MORTAR_BOUNDARY_OUT=43};
  getfem::mesh mesh;  /* the mesh */
  getfem::mesh_level_set mls;       /* the integration methods.              */
  getfem::mesh_im_level_set mim;    /* the integration methods.              */
  getfem::mesh_fem mf_pre_u, mf_pre_mortar;
  getfem::mesh_fem mf_mult;
  getfem::mesh_fem_level_set mfls_u, mfls_mortar;
  getfem::mesh_fem_global_function mf_sing_u;
  getfem::mesh_fem mf_partition_of_unity;
  getfem::mesh_fem_product mf_product;
  getfem::mesh_fem_sum mf_u_sum;
  
  base_small_vector translation;

  struct spider_param {
    getfem::spider_fem *fem;
    scalar_type theta0;
    scalar_type radius;
    unsigned Nr;
    unsigned Ntheta;
    int K;
    int bimat_enrichment;
    scalar_type epsilon;
  };

  spider_param spider;

  getfem::mesh_fem mf_us;

  getfem::mesh_fem& mf_u() { return mf_u_sum; }
  // getfem::mesh_fem& mf_u() { return mf_us; }
  
  scalar_type lambda, mu;    /* Lame coefficients.                */
  getfem::mesh_fem mf_rhs;   /* mesh_fem for the right hand side (f(x),..)   */
  getfem::mesh_fem mf_p;     /* mesh_fem for the pressure for mixed form     */
#ifdef VALIDATE_XFEM
  crack_exact_solution exact_sol;
#endif
  
  
  int bimaterial;           /* For bimaterial interface fracture */
  double lambda_up, lambda_down, mu_up, mu_down;  /*Lame coeff for bimaterial case*/
  getfem::level_set ls;      /* The two level sets defining the crack.       */
  getfem::level_set ls2, ls3; /* The two level-sets defining the add. cracks.*/
 
  scalar_type residual;      /* max residual for the iterative solvers      */
  bool mixed_pressure, add_crack;
  unsigned dir_with_mult;
  scalar_type cutoff_radius, cutoff_radius1, cutoff_radius0, enr_area_radius;
  
  size_type cutoff_func;

  typedef enum { NO_ENRICHMENT=0, 
		 FIXED_ZONE=1, 
		 GLOBAL_WITH_MORTAR=2,
		 GLOBAL_WITH_CUTOFF=3,
		 SPIDER_FEM_ALONE=4,
		 SPIDER_FEM_ENRICHMENT=5 } enrichment_option_enum;
  enrichment_option_enum enrichment_option;

  std::string datafilename;
  
  int reference_test;
  std::string GLOBAL_FUNCTION_MF, GLOBAL_FUNCTION_U;

  bgeot::md_param PARAM;

  bool solve(plain_vector &U);
  void init(void);
  crack_problem(void) : mls(mesh), mim(mls), 
			mf_pre_u(mesh), mf_pre_mortar(mesh), mf_mult(mesh),
			mfls_u(mls, mf_pre_u), mfls_mortar(mls, mf_pre_mortar), 
			mf_sing_u(mesh),
			mf_partition_of_unity(mesh),
			mf_product(mf_partition_of_unity, mf_sing_u),

			mf_u_sum(mesh), mf_us(mesh), mf_rhs(mesh), mf_p(mesh),
#ifdef VALIDATE_XFEM
			exact_sol(mesh), 
#endif
			ls(mesh, 1, true), ls2(mesh, 1, true),
			ls3(mesh, 1, true) {}

};

std::string name_of_dof(getfem::pdof_description dof) {
  char s[200];
  sprintf(s, "UnknownDof[%p]", (void*)dof);
  for (dim_type d = 0; d < 4; ++d) {
    if (dof == getfem::lagrange_dof(d)) {
      sprintf(s, "Lagrange[%d]", d); goto found;
    }
    if (dof == getfem::normal_derivative_dof(d)) {
      sprintf(s, "D_n[%d]", d); goto found;
    }
    if (dof == getfem::global_dof(d)) {
      sprintf(s, "GlobalDof[%d]", d);
    }
    if (dof == getfem::mean_value_dof(d)) {
      sprintf(s, "MeanValue[%d]", d);
    }
    if (getfem::dof_xfem_index(dof) != 0) {
      sprintf(s, "Xfem[idx:%d]", int(dof_xfem_index(dof)));
    }
    
    for (dim_type r = 0; r < d; ++r) {
      if (dof == getfem::derivative_dof(d, r)) {
	sprintf(s, "D_%c[%d]", "xyzuvw"[r], d); goto found;
      }
      for (dim_type t = 0; t < d; ++t) {
	if (dof == getfem::second_derivative_dof(d, r, t)) {
	  sprintf(s, "D2%c%c[%d]", "xyzuvw"[r], "xyzuvw"[t], d); 
	  goto found;
	}
      }
    }
  }
 found:
  return s;
}


/* Read parameters from the .param file, build the mesh, set finite element
 * and integration methods and selects the boundaries.
 */
void crack_problem::init(void) {
  std::string MESH_TYPE = PARAM.string_value("MESH_TYPE","Mesh type ");
  std::string FEM_TYPE  = PARAM.string_value("FEM_TYPE","FEM name");
  std::string INTEGRATION = PARAM.string_value("INTEGRATION",
					       "Name of integration method");
  std::string SIMPLEX_INTEGRATION = PARAM.string_value("SIMPLEX_INTEGRATION",
					 "Name of simplex integration method");
  std::string SINGULAR_INTEGRATION = PARAM.string_value("SINGULAR_INTEGRATION");

  add_crack = (PARAM.int_value("ADDITIONAL_CRACK", "An additional crack ?") != 0);
  enrichment_option = enrichment_option_enum(PARAM.int_value("ENRICHMENT_OPTION",
							     "Enrichment option"));
  cout << "MESH_TYPE=" << MESH_TYPE << "\n";
  cout << "FEM_TYPE="  << FEM_TYPE << "\n";
  cout << "INTEGRATION=" << INTEGRATION << "\n";

  bimaterial = int(PARAM.int_value("BIMATERIAL", "bimaterial interface crack"));
  
  if (bimaterial == 1){
    mu = PARAM.real_value("MU", "Lame coefficient mu"); 
    mu_up = PARAM.real_value("MU_UP", "Lame coefficient mu"); 
    mu_down = PARAM.real_value("MU_DOWN", "Lame coefficient mu"); 
    lambda_up = PARAM.real_value("LAMBDA_UP", "Lame Coef");
    lambda_down = PARAM.real_value("LAMBDA_DOWN", "Lame Coef");
    lambda = PARAM.real_value("LAMBDA", "Lame coefficient lambda");
  }
  else{

    mu = PARAM.real_value("MU", "Lame coefficient mu");
    lambda = PARAM.real_value("LAMBDA", "Lame coefficient lambda");
  }
  
  
  spider.radius = PARAM.real_value("SPIDER_RADIUS","spider_radius");
  spider.Nr = unsigned(PARAM.int_value("SPIDER_NR","Spider_Nr "));
  spider.Ntheta = unsigned(PARAM.int_value("SPIDER_NTHETA","Ntheta "));
  spider.K = int(PARAM.int_value("SPIDER_K","K "));
  spider.theta0 = 0;
  spider.bimat_enrichment = int(PARAM.int_value("SPIDER_BIMAT_ENRICHMENT","spider_bimat_enrichment"));

  /* The following constants are taken from the Chang and Xu paper of 2007 
     titled The singular stress field and stress intensity factors of a crack 
     terminating at a bimaterial interface published in IJMECSCI. epsilon 
     being the constant present in the asymptotic displacement analytic solution: 
     r^{1/2} * cos (\epsilon Ln(r)) * f(\theta) and  r^{1/2} * sin (\epsilon Ln(r)) * f(\theta).
  */
  
  if (spider.bimat_enrichment == 1){
    scalar_type nu1 = (lambda_up) / (2.*lambda_up + mu_up);
    scalar_type nu2 = (lambda_down) / (2.*lambda_down + mu_down);
    scalar_type kappa1 = 3. - 4. * nu1;
    scalar_type kappa2 = 3. - 4. * nu2;
    if (lambda_up == lambda_down && mu_up == mu_down)
      cout << "ERROR... Connot use the spider bimaterial enrichment with an isotropic homogenuous material (beta = 0/0!!!)... You should either use a bimaterial or disable the spider bimaterial enrichment" << endl;
    
    scalar_type beta = (mu_up*(kappa2-1.) - mu_down*(kappa1-1.)) / (mu_up*(kappa2+1.) - mu_down*(kappa1+1.));
      spider.epsilon = 1./(2.*M_PI) * log( (1.-beta) / (1.+beta) );
  
    
    //spider.epsilon = PARAM.real_value("SPIDER_BIMAT_ENRICHMENT","spider_bimat_enrichment");
  }
  
  translation.resize(2); 
  translation[0] =0.5;
  translation[1] =0.;

  /* First step : build the mesh */
  bgeot::pgeometric_trans pgt = 
    bgeot::geometric_trans_descriptor(MESH_TYPE);
  size_type N = pgt->dim();
  std::vector<size_type> nsubdiv(N);
  std::fill(nsubdiv.begin(),nsubdiv.end(),
	    PARAM.int_value("NX", "Nomber of space steps "));
  getfem::regular_unit_mesh(mesh, nsubdiv, pgt,
			    PARAM.int_value("MESH_NOISED") != 0);
  base_small_vector tt(N); tt[1] = -0.5;
  mesh.translation(tt); 
  

  scalar_type refinement_radius;
  refinement_radius = PARAM.real_value("REFINEMENT_RADIUS","Refinement Radius");
  size_type refinement_process;
  refinement_process = PARAM.int_value("REFINEMENT_PROCESS","Refinement process");
  dal::bit_vector conv_to_refine;
  size_type ref = 1;  
  if (refinement_radius > 0){
    while(ref <= refinement_process){
      conv_to_refine.clear();
      for(dal::bv_visitor i(mesh.convex_index()); !i.finished(); ++i){
      for(size_type j=0; j < 3; ++j)
	if(fabs(mesh.points()[mesh.ind_points_of_convex(i)[j]][0])<refinement_radius 
	   && fabs(mesh.points()[mesh.ind_points_of_convex(i)[j]][1])<refinement_radius){
	  conv_to_refine.add(i);
	}
      }
      mesh.Bank_refine(conv_to_refine);
      
      ref = ref + 1; 
      refinement_radius = refinement_radius/3.;
      if(refinement_radius > 1e-16 )
	cout<<"refining process step " << ref << "... refining "<< conv_to_refine.size() <<" convexes..." << endl ; 
      
    }
  cout<<"refining process complete." << endl ;
  }



  datafilename = PARAM.string_value("ROOTFILENAME","Base name of data files.");
  residual = PARAM.real_value("RESIDUAL");
  if (residual == 0.) residual = 1e-10;
  enr_area_radius = PARAM.real_value("RADIUS_ENR_AREA",
				     "radius of the enrichment area");
  
  GLOBAL_FUNCTION_MF = PARAM.string_value("GLOBAL_FUNCTION_MF");
  GLOBAL_FUNCTION_U = PARAM.string_value("GLOBAL_FUNCTION_U");

  reference_test = int(PARAM.int_value("REFERENCE_TEST", "Reference test")); 

  unsigned EXACT_SOL_NUM = unsigned(PARAM.int_value("EXACT_SOL_NUM", "Exact solution function number")); 
  

  cutoff_func = PARAM.int_value("CUTOFF_FUNC", "cutoff function");

  cutoff_radius = PARAM.real_value("CUTOFF", "Cutoff");
  cutoff_radius1 = PARAM.real_value("CUTOFF1", "Cutoff1");
  cutoff_radius0 = PARAM.real_value("CUTOFF0", "Cutoff0");
  mf_u().set_qdim(dim_type(N));

  /* set the finite element on the mf_u */
  getfem::pfem pf_u = 
    getfem::fem_descriptor(FEM_TYPE);
  getfem::pintegration_method ppi = 
    getfem::int_method_descriptor(INTEGRATION);
  getfem::pintegration_method simp_ppi = 
    getfem::int_method_descriptor(SIMPLEX_INTEGRATION);
  getfem::pintegration_method sing_ppi = (SINGULAR_INTEGRATION.size() ?
		getfem::int_method_descriptor(SINGULAR_INTEGRATION) : 0);
  
  mim.set_integration_method(mesh.convex_index(), ppi);
  mls.add_level_set(ls);
  if (add_crack) { mls.add_level_set(ls2); mls.add_level_set(ls3); }

  mim.set_simplex_im(simp_ppi, sing_ppi);
  mf_pre_u.set_finite_element(mesh.convex_index(), pf_u);
  mf_pre_mortar.set_finite_element(mesh.convex_index(), 
				   getfem::fem_descriptor(PARAM.string_value("MORTAR_FEM_TYPE")));
  mf_mult.set_finite_element(mesh.convex_index(), pf_u);
  mf_mult.set_qdim(dim_type(N));
  mf_partition_of_unity.set_classical_finite_element(1);
  
//   if (enrichment_option == 3 || enrichment_option == 4) {
//     spider = new getfem::spider_fem(spider_radius, mim, spider_Nr,
// 				    spider_Ntheta, spider_K, translation,
// 				    theta0);
//     mf_us.set_finite_element(mesh.convex_index(),spider->get_pfem());
//     for (dal::bv_visitor_c i(mf_us.convex_index()); !i.finished(); ++i) {
//       if (mf_us.fem_of_element(i)->nb_dof(i) == 0) {
// 	mf_us.set_finite_element(i,0);
//       }
//     }
//   }

  mixed_pressure =
    (PARAM.int_value("MIXED_PRESSURE","Mixed version or not.") != 0);
  dir_with_mult = unsigned(PARAM.int_value("DIRICHLET_VERSINO"));
  if (mixed_pressure) {
    std::string FEM_TYPE_P  = PARAM.string_value("FEM_TYPE_P","FEM name P");
    mf_p.set_finite_element(mesh.convex_index(),
			    getfem::fem_descriptor(FEM_TYPE_P));
  }

  /* set the finite element on mf_rhs (same as mf_u is DATA_FEM_TYPE is
     not used in the .param file */
  std::string data_fem_name = PARAM.string_value("DATA_FEM_TYPE");
  if (data_fem_name.size() == 0) {
    GMM_ASSERT1(pf_u->is_lagrange(), "You are using a non-lagrange FEM. "
		<< "In that case you need to set "
		<< "DATA_FEM_TYPE in the .param file");
    mf_rhs.set_finite_element(mesh.convex_index(), pf_u);
  } else {
    mf_rhs.set_finite_element(mesh.convex_index(), 
			      getfem::fem_descriptor(data_fem_name));
  }
  
  /* set boundary conditions
   * (Neuman on the upper face, Dirichlet elsewhere) */
  cout << "Selecting Neumann and Dirichlet boundaries\n";
  getfem::mesh_region border_faces;
  getfem::outer_faces_of_mesh(mesh, border_faces);
  for (getfem::mr_visitor i(border_faces); !i.finished(); ++i) {
    
    base_node un = mesh.normal_of_face_of_convex(i.cv(), i.f());
    un /= gmm::vect_norm2(un);
    if(bimaterial == 1) {
      
      if (un[0]  > 1.0E-7 ) { // new Neumann face
	mesh.region(DIRICHLET_BOUNDARY_NUM).add(i.cv(), i.f());
      } else {
	if (un[1]  > 1.0E-7 ) {
	  cout << "normal = " << un << endl;
	  mesh.region(NEUMANN_BOUNDARY_NUM1).add(i.cv(), i.f());
	}
	else {
	  if (un[1]  < -1.0E-7 ) {
	    cout << "normal = " << un << endl;
	    mesh.region(NEUMANN_BOUNDARY_NUM).add(i.cv(), i.f());
	  }
	  else {
	    if (un[0]  < -1.0E-7 ) {
	      cout << "normal = " << un << endl;
	      mesh.region(NEUMANN_HOMOGENE_BOUNDARY_NUM).add(i.cv(), i.f());
	    }
	  }
	}
      }
    }
    else {
      
#ifdef VALIDATE_XFEM
      mesh.region(DIRICHLET_BOUNDARY_NUM).add(i.cv(), i.f());
#else
      base_node un = mesh.normal_of_face_of_convex(i.cv(), i.f());
      un /= gmm::vect_norm2(un);
      if (un[0] - 1.0 < -1.0E-7) { // new Neumann face
	mesh.region(NEUMANN_BOUNDARY_NUM).add(i.cv(), i.f());
      } else {
	cout << "normal = " << un << endl;
	mesh.region(DIRICHLET_BOUNDARY_NUM).add(i.cv(), i.f());
      }
#endif
    }
  }
  
  
  
#ifdef VALIDATE_XFEM
  exact_sol.init(EXACT_SOL_NUM, lambda, mu, ls);
#endif
}


base_small_vector ls_function(const base_node P, int num = 0) {
  scalar_type x = P[0], y = P[1];
  base_small_vector res(2);
  switch (num) {
    case 0: {
      res[0] = y;
      res[1] = -.5 + x;
    } break;
    case 1: {
      res[0] = gmm::vect_dist2(P, base_node(0.5, 0.)) - .25;
      res[1] = gmm::vect_dist2(P, base_node(0.25, 0.0)) - 0.27;
    } break;
    case 2: {
      res[0] = x - 0.25;
      res[1] = gmm::vect_dist2(P, base_node(0.25, 0.0)) - 0.35;
    } break;
    default: assert(0);
  }
  return res;
}

bool crack_problem::solve(plain_vector &U) {
  size_type nb_dof_rhs = mf_rhs.nb_dof();
  size_type N = mesh.dim();
  ls.reinit();  
  cout << "ls.get_mesh_fem().nb_dof() = " << ls.get_mesh_fem().nb_dof() << "\n";
  for (size_type d = 0; d < ls.get_mesh_fem().nb_dof(); ++d) {
    ls.values(0)[d] = ls_function(ls.get_mesh_fem().point_of_basic_dof(d), 0)[0];
    ls.values(1)[d] = ls_function(ls.get_mesh_fem().point_of_basic_dof(d), 0)[1];
  }
  ls.touch();

  if (add_crack) {
    ls2.reinit();
    for (size_type d = 0; d < ls2.get_mesh_fem().nb_dof(); ++d) {
      ls2.values(0)[d] = ls_function(ls2.get_mesh_fem().point_of_basic_dof(d), 1)[0];
      ls2.values(1)[d] = ls_function(ls2.get_mesh_fem().point_of_basic_dof(d), 1)[1];
    }
    ls2.touch();
    
    ls3.reinit();
    for (size_type d = 0; d < ls3.get_mesh_fem().nb_dof(); ++d) {
      ls3.values(0)[d] = ls_function(ls2.get_mesh_fem().point_of_basic_dof(d), 2)[0];
      ls3.values(1)[d] = ls_function(ls2.get_mesh_fem().point_of_basic_dof(d), 2)[1];
    }
    ls3.touch(); 
  }

  mls.adapt();
  mim.adapt();
  mfls_u.adapt();
  mfls_mortar.adapt(); mfls_mortar.set_qdim(2);

  bool load_global_fun = GLOBAL_FUNCTION_MF.size() != 0;
 
  cout << "Setting up the singular functions for the enrichment\n";

  std::vector<getfem::pglobal_function> vfunc(4);
  if (!load_global_fun) {
    cout << "Using default singular functions\n";
    for (unsigned i = 0; i < vfunc.size(); ++i){
      /* use the singularity */
      getfem::abstract_xy_function *s = 
	new getfem::crack_singular_xy_function(i);
      if (enrichment_option != FIXED_ZONE && 
	  enrichment_option != GLOBAL_WITH_MORTAR) {
	/* use the product of the singularity function
	   with a cutoff */
	getfem::abstract_xy_function *c = 
	  new getfem::cutoff_xy_function(int(cutoff_func),
					 cutoff_radius, 
					 cutoff_radius1,cutoff_radius0);
	s = new getfem::product_of_xy_functions(*s, *c);
      }
      vfunc[i] = getfem::global_function_on_level_set(ls, *s);
    }
  } else {
    cout << "Load singular functions from " << GLOBAL_FUNCTION_MF << " and " << GLOBAL_FUNCTION_U << "\n";
    getfem::mesh *m = new getfem::mesh(); 
    m->read_from_file(GLOBAL_FUNCTION_MF);
    getfem::mesh_fem *mf_c = new getfem::mesh_fem(*m); 
    mf_c->read_from_file(GLOBAL_FUNCTION_MF);
    std::fstream f(GLOBAL_FUNCTION_U.c_str(), std::ios::in);
    plain_vector W(mf_c->nb_dof());


  
    for (unsigned i=0; i < mf_c->nb_dof(); ++i) {
      f >> W[i]; GMM_ASSERT1(f.good(), "problem while reading " << GLOBAL_FUNCTION_U);
      
      //cout << "The precalculated dof " << i << " of coordinates " << mf_c->point_of_dof(i) << " is "<< W[i] <<endl; 
      /*scalar_type x = pow(mf_c->point_of_dof(i)[0],2); scalar_type y = pow(mf_c->point_of_dof(i)[1],2);
	scalar_type r = std::sqrt(pow(x,2) + pow(y,2));
	scalar_type sgny = (y < 0 ? -1.0 : 1.0);
	scalar_type sin2 = sqrt(gmm::abs(.5-x/(2*r))) * sgny;
	scalar_type cos2 = sqrt(gmm::abs(.5+x/(2*r)));
	W[i] = std::sqrt(r) * sin2;
      */
    }
    unsigned nb_func = mf_c->get_qdim();
    cout << "read " << nb_func << " global functions OK.\n";
    vfunc.resize(nb_func);
    getfem::interpolator_on_mesh_fem *global_interp = 
      new getfem::interpolator_on_mesh_fem(*mf_c, W);
    for (size_type i=0; i < nb_func; ++i) {
      /* use the precalculated function for the enrichment*/
      //getfem::abstract_xy_function *s = new getfem::crack_singular_xy_function(i);
      getfem::abstract_xy_function *s = new getfem::interpolated_xy_function(*global_interp,i);

      if (enrichment_option != FIXED_ZONE && 
	  enrichment_option != GLOBAL_WITH_MORTAR) {

	/* use the product of the enrichment function
	   with a cutoff */
	getfem::abstract_xy_function *c = 
	  new getfem::cutoff_xy_function(int(cutoff_func),
					 cutoff_radius, 
					 cutoff_radius1,cutoff_radius0);
	s = new getfem::product_of_xy_functions(*s, *c);
      }    
      vfunc[i] = getfem::global_function_on_level_set(ls, *s);
    }    
  }
  
  
  mf_sing_u.set_functions(vfunc);
  

  if (enrichment_option == SPIDER_FEM_ALONE || 
      enrichment_option == SPIDER_FEM_ENRICHMENT) {
    spider.fem = new getfem::spider_fem(spider.radius, mim, spider.Nr,
					spider.Ntheta, spider.K, translation,
					spider.theta0, spider.bimat_enrichment,
					spider.epsilon);
    mf_us.set_finite_element(mesh.convex_index(),spider.fem->get_pfem());
    for (dal::bv_visitor_c i(mf_us.convex_index()); !i.finished(); ++i) {
      if (mf_us.fem_of_element(i)->nb_dof(i) == 0) {
	mf_us.set_finite_element(i,0);
      }
    }
    spider.fem->check();
  }

  switch (enrichment_option) {

  case FIXED_ZONE :
    {
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
	GMM_WARNING0("There is " << enriched_dofs.card() <<
		     " enriched dofs for the crack tip");
      mf_product.set_enrichment(enriched_dofs);
      mf_u_sum.set_mesh_fems(mf_product, mfls_u);
    }
    break;
  

    case GLOBAL_WITH_MORTAR: {
      // Selecting the element in the enriched domain

      dal::bit_vector cvlist_in_area;
      dal::bit_vector cvlist_out_area;
      for (dal::bv_visitor cv(mesh.convex_index()); 
	   !cv.finished(); ++cv) {
	bool in_area = true;
	/* For each element, we test all of its nodes. 
	   If all the nodes are inside the enrichment area,
	   then the element is completly inside the area too */ 
	for (unsigned j=0; j < mesh.nb_points_of_convex(cv); ++j) {
	  if (gmm::sqr(mesh.points_of_convex(cv)[j][0] - translation[0]) + 
	      gmm::sqr(mesh.points_of_convex(cv)[j][1] - translation[1]) > 
	      gmm::sqr(enr_area_radius)) {
	    in_area = false; break;
	  }
	}

	/* "remove" the global function on convexes outside the enrichment
	   area */
	if (!in_area) {
	  cvlist_out_area.add(cv);
	  mf_sing_u.set_finite_element(cv, 0);
	  mf_u().set_dof_partition(cv, 1);
	} else cvlist_in_area.add(cv);
      }

      /* extract the boundary of the enrichment area, from the
	 "inside" point-of-view, and from the "outside"
	 point-of-view */
      getfem::mesh_region r_border, r_enr_out;
      getfem::outer_faces_of_mesh(mesh, r_border);

      getfem::outer_faces_of_mesh(mesh, cvlist_in_area, 
				  mesh.region(MORTAR_BOUNDARY_IN));
      getfem::outer_faces_of_mesh(mesh, cvlist_out_area, 
				  mesh.region(MORTAR_BOUNDARY_OUT));
      for (getfem::mr_visitor v(r_border); !v.finished(); ++v) {
	mesh.region(MORTAR_BOUNDARY_OUT).sup(v.cv(), v.f());
      }
      mf_u_sum.set_mesh_fems(mf_sing_u, mfls_u);
    } break;

  case GLOBAL_WITH_CUTOFF :{
    if(cutoff_func == 0)
      cout<<"Using exponential Cutoff..."<<endl;
    else
      cout<<"Using Polynomial Cutoff..."<<endl;
    mf_u_sum.set_mesh_fems(mf_sing_u, mfls_u);
  } break;

  case SPIDER_FEM_ALONE : { 
    mf_u_sum.set_mesh_fems(mf_us); 
  } break;
    
  case SPIDER_FEM_ENRICHMENT : {
    mf_u_sum.set_mesh_fems(mf_us, mfls_u); 
  } break;
    
  case NO_ENRICHMENT: {
    mf_u_sum.set_mesh_fems(mfls_u);
  } break;
  
  }
  

  U.resize(mf_u().nb_dof());


  if (mixed_pressure) cout << "Number of dof for P: " << mf_p.nb_dof() << endl;
  cout << "Number of dof for u: " << mf_u().nb_dof() << endl;

  // Linearized elasticity brick.
  getfem::mdbrick_isotropic_linearized_elasticity<>
    ELAS(mim, mf_u(), mixed_pressure ? 0.0 : lambda, mu);

  
  if(bimaterial == 1){
    cout<<"______________________________________________________________________________"<<endl;
    cout<<"CASE OF BIMATERIAL CRACK  with lambda_up = "<<lambda_up<<" and lambda_down = "<<lambda_down<<endl;
    cout<<"______________________________________________________________________________"<<endl;
    std::vector<double> bi_lambda(ELAS.lambda().mf().nb_dof());
    std::vector<double> bi_mu(ELAS.lambda().mf().nb_dof());
    
    cout<<"ELAS.lambda().mf().nb_dof()==="<<ELAS.lambda().mf().nb_dof()<<endl;
    GMM_ASSERT1(!ELAS.lambda().mf().is_reduced(), "To be adapted");

    for (size_type ite = 0; ite < ELAS.lambda().mf().nb_dof();ite++) {
      if (ELAS.lambda().mf().point_of_basic_dof(ite)[1] > 0){
	bi_lambda[ite] = lambda_up;
	bi_mu[ite] = mu_up;
      }
	else{
	  bi_lambda[ite] = lambda_down;
	  bi_mu[ite] = mu_down;
	}
    }
    
    //cout<<"bi_lambda.size() = "<<bi_lambda.size()<<endl;
    // cout<<"ELAS.lambda().mf().nb_dof()==="<<ELAS.lambda().mf().nb_dof()<<endl;
    
    ELAS.lambda().set(bi_lambda);
    ELAS.mu().set(bi_mu);
  }
  

  getfem::mdbrick_abstract<> *pINCOMP;
  if (mixed_pressure) {
    getfem::mdbrick_linear_incomp<> *incomp
      = new getfem::mdbrick_linear_incomp<>(ELAS, mf_p);
    incomp->penalization_coeff().set(1.0/lambda);
    pINCOMP = incomp;
  } else pINCOMP = &ELAS;

  // Defining the volumic source term.
  plain_vector F(nb_dof_rhs * N);
  getfem::interpolation_function(mf_rhs, F, sol_f);
  
  // Volumic source term brick.
  getfem::mdbrick_source_term<> VOL_F(*pINCOMP, mf_rhs, F);

  // Defining the Neumann condition right hand side.
  gmm::clear(F);
  
  // Neumann condition brick.
  
  getfem::mdbrick_abstract<> *pNEUMANN;
  
  
  if(bimaterial ==  1){
    //down side
    for(size_type i = 1; i<F.size(); i=i+2) 
      F[i] = -0.4;
    for(size_type i = 0; i<F.size(); i=i+2) 
      F[i] = -0.2;
  }
  
  getfem::mdbrick_source_term<>  NEUMANN(VOL_F, mf_rhs, F,NEUMANN_BOUNDARY_NUM);   
  //left side (crack opening side)
  gmm::clear(F);
  for(size_type i = 1; i<F.size(); i=i+2) 
    F[i] = 0.;
  for(size_type i = 0; i<F.size(); i=i+2) 
    F[i] = -0.;
  getfem::mdbrick_source_term<> NEUMANN_HOM(NEUMANN, mf_rhs, F,NEUMANN_HOMOGENE_BOUNDARY_NUM);
   
    //upper side
  gmm::clear(F);
  for(size_type i = 1; i<F.size(); i=i+2) 
    F[i] = 0.4;
  for(size_type i = 0; i<F.size(); i=i+2) 
    F[i] = 0.;
  getfem::mdbrick_source_term<> NEUMANN1(NEUMANN_HOM, mf_rhs, F,NEUMANN_BOUNDARY_NUM1);
  
  if (bimaterial == 1)
    pNEUMANN = & NEUMANN1;
  else
    pNEUMANN = & NEUMANN;
  
  
  
  //toto_solution toto(mf_rhs.linked_mesh()); toto.init();
  //assert(toto.mf.nb_dof() == 1);
  
  // Dirichlet condition brick.
  getfem::mdbrick_Dirichlet<> DIRICHLET(*pNEUMANN, DIRICHLET_BOUNDARY_NUM, mf_mult);
  
  if (bimaterial == 1)
    DIRICHLET.rhs().set(exact_sol.mf,0);
  else {
#ifdef VALIDATE_XFEM
    DIRICHLET.rhs().set(exact_sol.mf,exact_sol.U);
#endif
  }
  DIRICHLET.set_constraints_type(getfem::constraints_type(dir_with_mult));

  getfem::mdbrick_abstract<> *final_model = &DIRICHLET;

 if (enrichment_option == GLOBAL_WITH_MORTAR) {
    /* add a constraint brick for the mortar junction between
       the enriched area and the rest of the mesh */
    /* we use mfls_u as the space of lagrange multipliers */
    getfem::mesh_fem &mf_mortar = mfls_mortar; 
    /* adjust its qdim.. this is just evil and dangerous
       since mf_u() is built upon mfls_u.. it would be better
       to use a copy. */
    mf_mortar.set_qdim(2); // EVIL 
    getfem::mdbrick_constraint<> &mortar = 
      *(new getfem::mdbrick_constraint<>(DIRICHLET,0));

    cout << "Handling mortar junction\n";

    /* list of dof of mf_mortar for the mortar condition */
    std::vector<size_type> ind_mortar;
    /* unfortunately , dof_on_region sometimes returns too much dof
       when mf_mortar is an enriched one so we have to filter them */
    GMM_ASSERT1(!mf_mortar.is_reduced(), "To be adapted");
    sparse_matrix M(mf_mortar.nb_dof(), mf_mortar.nb_dof());
    getfem::asm_mass_matrix(M, mim, mf_mortar, MORTAR_BOUNDARY_OUT);
    for (dal::bv_visitor_c d(mf_mortar.basic_dof_on_region(MORTAR_BOUNDARY_OUT)); 
	 !d.finished(); ++d) {
      if (M(d,d) > 1e-8) ind_mortar.push_back(d);
      else cout << "  removing non mortar dof" << d << "\n";
    }

    cout << ind_mortar.size() << " dof for the lagrange multiplier)\n";

    sparse_matrix H0(mf_mortar.nb_dof(), mf_u().nb_dof()), 
      H(ind_mortar.size(), mf_u().nb_dof());


    gmm::sub_index sub_i(ind_mortar);
    gmm::sub_interval sub_j(0, mf_u().nb_dof());
    
    /* build the mortar constraint matrix -- note that the integration
       method is conformal to the crack
     */
    getfem::asm_mass_matrix(H0, mim, mf_mortar, mf_u(), 
			    MORTAR_BOUNDARY_OUT);
    gmm::copy(gmm::sub_matrix(H0, sub_i, sub_j), H);

    gmm::clear(H0);
    getfem::asm_mass_matrix(H0, mim, mf_mortar, mf_u(), 
			    MORTAR_BOUNDARY_IN);
    gmm::add(gmm::scaled(gmm::sub_matrix(H0, sub_i, sub_j), -1), H);


    /* because of the discontinuous partition of mf_u(), some levelset
       enriched functions do not contribute any more to the
       mass-matrix (the ones which are null on one side of the
       levelset, when split in two by the mortar partition, may create
       a "null" dof whose base function is all zero..
    */
    sparse_matrix M2(mf_u().nb_dof(), mf_u().nb_dof());
    getfem::asm_mass_matrix(M2, mim, mf_u(), mf_u());

    for (size_type d = 0; d < mf_u().nb_dof(); ++d) {
      if (M2(d,d) < 1e-10) {
	cout << "  removing null mf_u() dof " << d << "\n";
	size_type n = gmm::mat_nrows(H);
	gmm::resize(H, n+1, gmm::mat_ncols(H));
	H(n, d) = 1;
      }
    }

    
    
    getfem::base_vector R(gmm::mat_nrows(H));
    mortar.set_constraints(H,R);

    final_model = &mortar;
  }

  // Generic solve.
  cout << "Total number of variables : " << final_model->nb_dof() << endl;
  getfem::standard_model_state MS(*final_model);
  gmm::iteration iter(residual, 1, 40000);
  cout << "Solving..." << endl;
  getfem::standard_solve(MS, *final_model, iter);
  cout << "Solving... done" << endl;
  // Solution extraction
  gmm::copy(ELAS.get_solution(MS), U);
  
  if(reference_test)
    {
      cout << "Exporting reference solution...";
      GMM_ASSERT1(!mf_u().is_reduced(), "To be adapted");
      dal::bit_vector blocked_dof = mf_u().basic_dof_on_region(5);
      getfem::mesh_fem mf_refined(mesh, dim_type(N));
      std::string FEM_DISC = PARAM.string_value("FEM_DISC","fem disc ");
      mf_refined.set_finite_element(mesh.convex_index(),
				    getfem::fem_descriptor(FEM_DISC));
      
      plain_vector W(mf_refined.nb_dof());
      getfem::interpolation(mf_u(), mf_refined, U, W);
      
      
      mf_refined.write_to_file(datafilename + "_refined_test.meshfem", true);
      gmm::vecsave(datafilename + "_refined_test.U", W);
      cout << "done" << endl;
    }
  
  
  return (iter.converged());
}

void export_interpolated_on_line(const getfem::mesh_fem &mf,
				 const getfem::base_vector &U,
				 const base_node &x0,
				 const base_small_vector &dir,
				 const int nb_points, 
				 const std::string &filename) {
  getfem::mesh_trans_inv mti(mf.linked_mesh());
  scalar_type h = 1.0/(2*nb_points);
  for (int i=-nb_points; i <= nb_points; ++i) {
    mti.add_point(x0 + 2*(i*h)*dir);
  }

   getfem::base_vector V(mti.nb_points() * mf.get_qdim());
  getfem::base_matrix M;
  getfem::interpolation(mf, mti, U, V, M, 0, false);
  
  std::ofstream f(filename.c_str()); f.precision(16);
  
  for (size_type i=0; i < mti.nb_points(); ++i) {
    for (unsigned q=0; q < mf.get_qdim(); ++q) {
      f << V[i*mf.get_qdim()+q] << " ";
    }
    f << "\n";
  }
}
				 

/**************************************************************************/
/*  main program.                                                         */
/**************************************************************************/

int main(int argc, char *argv[]) {
  
  GMM_SET_EXCEPTION_DEBUG; // Exceptions make a memory fault, to debug.
  FE_ENABLE_EXCEPT;        // Enable floating point exception for Nan.

  //getfem::getfem_mesh_level_set_noisy();


  try {
    crack_problem p;
    p.PARAM.read_command_line(argc, argv);
    p.init();
    p.mesh.write_to_file(p.datafilename + ".mesh");

    plain_vector U(p.mf_u().nb_dof());
    if (!p.solve(U)) GMM_ASSERT1(false,"Solve has failed");
 
    //        for (size_type i = 4; i < U.size(); ++i)
    //U[i] = 0;
    //cout << "The solution" << U ;
     gmm::vecsave("crack.U", U);
     cout << "vecsave done"<<endl;
    { 
      getfem::mesh mcut;
      p.mls.global_cut_mesh(mcut);
      unsigned Q = p.mf_u().get_qdim();
      getfem::mesh_fem mf(mcut, dim_type(Q));
      mf.set_classical_discontinuous_finite_element(2, 0.001);
      // mf.set_finite_element
      //	(getfem::fem_descriptor("FEM_PK_DISCONTINUOUS(2, 2, 0.0001)"));
      plain_vector V(mf.nb_dof());

      getfem::interpolation(p.mf_u(), mf, U, V);

      getfem::stored_mesh_slice sl;
      getfem::mesh mcut_refined;

      unsigned NX = unsigned(p.PARAM.int_value("NX")), nn;
//  unsigned NX = (unsigned)sqrt(p.mesh.convex_index().card()); 
//    unsigned nn;
      if (NX < 6) nn = 24;
      else if (NX < 12) nn = 6;
      else if (NX < 30) nn = 3;
      else nn = 3;

      /* choose an adequate slice refinement based on the distance to the crack tip */
      std::vector<bgeot::short_type> nrefine(mcut.convex_index().last_true()+1);
      for (dal::bv_visitor cv(mcut.convex_index()); !cv.finished(); ++cv) {
	scalar_type dmin=0, d;
	base_node Pmin,P;
	for (unsigned i=0; i < mcut.nb_points_of_convex(cv); ++i) {
	  P = mcut.points_of_convex(cv)[i];
	  d = gmm::vect_norm2(ls_function(P));
	  if (d < dmin || i == 0) { dmin = d; Pmin = P; }
	}

	if (dmin < 1e-5)
	  nrefine[cv] = short_type(nn*8);
	else if (dmin < .1) 
	  nrefine[cv] = short_type(nn*2);
	else nrefine[cv] = short_type(nn);
	if (dmin < .01)
	  cout << "cv: "<< cv << ", dmin = " << dmin << "Pmin=" << Pmin << " " << nrefine[cv] << "\n";
      }

      {
	getfem::mesh_slicer slicer(mcut); 
	getfem::slicer_build_mesh bmesh(mcut_refined);
	slicer.push_back_action(bmesh);
	slicer.exec(nrefine, getfem::mesh_region::all_convexes());
      }
      /*
      sl.build(mcut, 
      getfem::slicer_build_mesh(mcut_refined), nrefine);*/

      getfem::mesh_im mim_refined(mcut_refined); 
      mim_refined.set_integration_method(getfem::int_method_descriptor
					 ("IM_TRIANGLE(6)"));

      getfem::mesh_fem mf_refined(mcut_refined, dim_type(Q));
      mf_refined.set_classical_discontinuous_finite_element(2, 0.001);
      plain_vector W(mf_refined.nb_dof());

      getfem::interpolation(p.mf_u(), mf_refined, U, W);


#ifdef VALIDATE_XFEM
      p.exact_sol.mf.set_qdim(dim_type(Q));
      assert(p.exact_sol.mf.nb_dof() == p.exact_sol.U.size());
      plain_vector EXACT(mf_refined.nb_dof());
      getfem::interpolation(p.exact_sol.mf, mf_refined, 
			    p.exact_sol.U, EXACT);

      plain_vector DIFF(EXACT); gmm::add(gmm::scaled(W,-1),DIFF);
#endif

      if (p.PARAM.int_value("VTK_EXPORT")) {
	getfem::mesh_fem mf_refined_vm(mcut_refined, 1);
	mf_refined_vm.set_classical_discontinuous_finite_element(1, 0.001);
	cerr << "mf_refined_vm.nb_dof=" << mf_refined_vm.nb_dof() << "\n";
	plain_vector VM(mf_refined_vm.nb_dof());

	cout << "computing von mises\n";
	getfem::interpolation_von_mises(mf_refined, mf_refined_vm, W, VM);

	plain_vector D(mf_refined_vm.nb_dof() * Q), 
	  DN(mf_refined_vm.nb_dof());
	
#ifdef VALIDATE_XFEM
	getfem::interpolation(mf_refined, mf_refined_vm, DIFF, D);
	for (unsigned i=0; i < DN.size(); ++i) {
	  DN[i] = gmm::vect_norm2(gmm::sub_vector(D, gmm::sub_interval(i*Q, Q)));
	}
#endif

	cout << "export to " << p.datafilename + ".vtk" << "..\n";
	getfem::vtk_export exp(p.datafilename + ".vtk",
			       p.PARAM.int_value("VTK_EXPORT")==1);

	exp.exporting(mf_refined); 
	//exp.write_point_data(mf_refined_vm, DN, "error");
	exp.write_point_data(mf_refined_vm, VM, "von mises stress");

	exp.write_point_data(mf_refined, W, "elastostatic_displacement");

	base_node line_x0(0.70001,0);
	base_small_vector line_dir(0, 0.5001);
	unsigned line_nb_points = 1000;
	export_interpolated_on_line(mf_refined_vm, VM, 
				    line_x0, line_dir, line_nb_points,
				    "von_mises_on_line.data");

#ifdef VALIDATE_XFEM

	plain_vector VM_EXACT(mf_refined_vm.nb_dof());


	/* getfem::mesh_fem_global_function mf(mcut_refined,Q);
	   std::vector<getfem::pglobal_function> cfun(4);
	   for (unsigned j=0; j < 4; ++j)
	   cfun[j] = getfem::isotropic_crack_singular_2D(j, p.ls);
	   mf.set_functions(cfun);
	   getfem::interpolation_von_mises(mf, mf_refined_vm, p.exact_sol.U,
	   VM_EXACT);
	*/


	getfem::interpolation_von_mises(mf_refined, mf_refined_vm, EXACT, VM_EXACT);
	getfem::vtk_export exp2("crack_exact.vtk");
	exp2.exporting(mf_refined);
	exp2.write_point_data(mf_refined_vm, VM_EXACT, "exact von mises stress");
	exp2.write_point_data(mf_refined, EXACT, "reference solution");
	
	export_interpolated_on_line(mf_refined_vm, VM_EXACT, 
				    line_x0, line_dir, line_nb_points,
				    "von_mises_on_line_exact.data");
#endif
	
	cout << "export done, you can view the data file with (for example)\n"
	  "mayavi -d " << p.datafilename << ".vtk -f  "
	  "WarpVector -m BandedSurfaceMap -m Outline\n";
      }

      //bool error_to_ref_sol = 0;
      
      if(p.PARAM.int_value("ERROR_TO_REF_SOL") == 1){
	cout << "Coputing error with respect to a reference solution..." << endl;
	std::string REFERENCE_MF = "crack_refined_test.meshfem";
	std::string REFERENCE_U = "crack_refined_test.U";
	
	cout << "Load reference meshfem and solution from " << REFERENCE_MF << " and " << REFERENCE_U << "\n";
	getfem::mesh ref_m; 
	ref_m.read_from_file(REFERENCE_MF);
	getfem::mesh_fem ref_mf(ref_m); 
	ref_mf.read_from_file(REFERENCE_MF);
	std::fstream f(REFERENCE_U.c_str(), std::ios::in);
	plain_vector ref_U(ref_mf.nb_dof());
      
	for (unsigned i=0; i < ref_mf.nb_dof(); ++i) {
	  f >> ref_U[i]; if (!f.good()) GMM_ASSERT1(f.good(), "problem while reading " << REFERENCE_U);
	}
	
	getfem::mesh_im ref_mim(ref_m);
	getfem::pintegration_method ppi = 
	  getfem::int_method_descriptor("IM_TRIANGLE(6)");
	ref_mim.set_integration_method(ref_m.convex_index(), ppi);
	plain_vector interp_U(ref_mf.nb_dof());
	getfem::interpolation(p.mf_u(), ref_mf, U, interp_U);
	
	//gmm::add(gmm::scaled(interp_U, -1), ref_U);
	  
	
	cout << "To ref L2 ERROR:" << getfem::asm_L2_dist(ref_mim, ref_mf, interp_U,
							  ref_mf, ref_U) << endl;
	
	cout << "To ref H1 ERROR:" << getfem::asm_H1_dist(ref_mim, ref_mf, interp_U,
							  ref_mf, ref_U) << endl;
	
	//cout << "L2 ERROR:"<< getfem::asm_L2_norm(ref_mim, ref_mf, ref_U)
	//     << endl << "H1 ERROR:"
	//     << getfem::asm_H1_norm(ref_mim, ref_mf, ref_U) << "\n";
	
      }
      
#ifdef VALIDATE_XFEM

      else {
	cout << "L2 ERROR:"<< getfem::asm_L2_dist(p.mim, p.mf_u(), U,
						  p.exact_sol.mf, p.exact_sol.U)
	     << endl << "H1 ERROR:"
	     << getfem::asm_H1_dist(p.mim, p.mf_u(), U,
				    p.exact_sol.mf, p.exact_sol.U) << "\n";
	
      }

      cout << "L2 norm of the solution:"  << getfem::asm_L2_norm(p.mim,p.mf_u(),U)<<endl;
      cout << "H1 norm of the solution:"  << getfem::asm_H1_norm(p.mim,p.mf_u(),U)<<endl;
      
      

      /* cout << "OLD ERROR L2:" 
	 << getfem::asm_L2_norm(mim_refined,mf_refined,DIFF) 
	 << " H1:" << getfem::asm_H1_dist(mim_refined,mf_refined,
	 EXACT,mf_refined,W)  << "\n";
	 
	 cout << "ex = " << p.exact_sol.U << "\n";
	 cout << "U  = " << gmm::sub_vector(U, gmm::sub_interval(0,8)) << "\n";
      */
#endif
    }

  }
  GMM_STANDARD_CATCH_ERROR;

  return 0; 
}
